package prog

import docking.dockscore._
import docking.docksearch._
import docking.docksearch.forcevector.{ForceVectorDocker, ForceVectorScore}
import io.Mol2Writer
import io.threadcso._
import jmolint.JmolCmds._
import jmolint.{JmolCmds, JmolFrame, JmolMoleculeReader, JmolPanel}
import model.Molecule
import opt.Action
import profiling.Profiler


object DockMain {
  val usage = "USAGE: scala DockMain -a (path to A, pdb format) " +
    "-b (path to B, pdb format) " +
    "-out (B's output path, xyz format) " +
    " [-dir d] [--randominit] [-ref a,b,c,...]" +
    "-docker (atompair|ehc|forcevector) [--consolelog] [--nogui] " +
    " [-initangle a] [initlevel i]" +
    " [-maxiters m]" +
    " [--ignorehydrogen] " +
    " [-surface s] " +
    " [-permeability p]" +
    " [-balance atomic,electric,hydrogenbond,bondstrength] " +
    " [-threshold t] "

  val frame = new JmolFrame(500, 500, false)
  val jmolPanel = frame.getPanel

  val viewInitCmds = getViewInitCmds

  def main(args: Array[String]): Unit = {
    parseArgs(args)
    val (docked, rmsd, _) = doMainDock(DockArgs.fullPathA, DockArgs.fullPathB, DockArgs.fullPathOut)

    new Mol2Writer(DockArgs.fullPathOut).write(docked)      // write docked b to file
    jmolPanel.openFiles(List(DockArgs.fullPathA, DockArgs.fullPathOut))  // show original a and modified b
    jmolPanel.execSeq(viewInitCmds)
    println(s"Finished with RMSD: $rmsd")
    Profiler.report
  }


  def doMainDock(pathA: String, pathB: String, pathOut: String) = {
    //val startTime = System.currentTimeMillis()

    jmolPanel.openFiles(List(pathA, pathB))
    jmolPanel.execSync(
      selectModel("2.1"),
      setLog(0),
      zoom(50),
      save
    )

    val molA: Molecule = JmolMoleculeReader.read(jmolPanel, 0)
    val molB: Molecule = JmolMoleculeReader.read(jmolPanel, 1)
//    println(s"Pairs: ${molA.atoms(DockArgs.ignoreAHydrogens).size * molB.surfaceAtoms.size}")
    val docker = getDocker

    val chan = OneOneBuf[Any](5)
    var dockResult = (null.asInstanceOf[Molecule], 0.0)
    (proc { dockResult = docker.dock(molA, molB.clone, chan); chan.close } ||
      showActions(chan, jmolPanel))()

    val docked = dockResult._1
    val score = dockResult._2
    (docked, getRMSD(docked, molB), score)
  }

  def showActions(chan: ?[Any], panel: JmolPanel) = proc {
    repeat {
      val msg = chan?();
      if (DockArgs.liveGui) {
        val s = msg match
        {
          case a: Action => val cmd = JmolCmds.cmd(a); panel.exec(cmd); cmd
          case "save" => panel.exec(JmolCmds.save); "save"
          case "reset" => panel.exec(JmolCmds.reset); panel.execSeq(viewInitCmds); "reset"
          case other => other.toString
        }
        if (DockArgs.consoleLog)
          println(s)
        ()
      }
    }
  }

  def getRMSD(docked: Molecule, molB: Molecule) = {
    if (DockArgs.fullPathsRef.isEmpty)
      docked.rmsd(molB)
    else {
      jmolPanel.openFiles(DockArgs.fullPathsRef)
      (for (i <- DockArgs.fullPathsRef.indices) yield {
        val refMol = JmolMoleculeReader.read(jmolPanel, i)
        docked.rmsd(refMol)
      }).min
    }
  }


  /**
    * Gets an instance of a docker given the name. In the future, this could be
    * enhanced to use reflection.
    */
  def getDocker = {
    val innerDocker = DockArgs.dockerName match {
      case "atompair" => new AtomPairDocker(new SurfaceDistanceScorer(DockArgs.surface))
      case "ehc" =>
        val scorer =
          if (DockArgs.scorerName == "ff")
            new ForceVectorScore(DockArgs.surface, DockArgs.ignoreAHydrogens, 1.25, DockArgs.geometricForceWeight, DockArgs.electricForceWeight, DockArgs.hydrogenBondsForceWeight, DockArgs.bondForceWeight)
          else
            new SurfaceDistanceScorer(DockArgs.surface)
        new EhcDocker(scorer, DockArgs.maxIters)

      case "forcevector" => new ForceVectorDocker(
        surface = DockArgs.surface,
        permeability = DockArgs.permeability,
        maxDecelerations = 10,
        ignoreAHydrogens = DockArgs.ignoreAHydrogens,
        threshold = DockArgs.threshold,
        geometricForceWeight = DockArgs.geometricForceWeight,
        electricForceWeight = DockArgs.electricForceWeight,
        hydrogenBondsForceWeight = DockArgs.hydrogenBondsForceWeight,
        bondForceWeight = DockArgs.bondForceWeight
      )

      case "forcevectorc" => new ForceVectorConcurrentDockerDocker(
        surface = DockArgs.surface,
        permeability = DockArgs.permeability,
        maxDecelerations = 10,
        ignoreAHydrogens = DockArgs.ignoreAHydrogens,
        threshold = DockArgs.threshold,
        geometricForceWeight = DockArgs.geometricForceWeight,
        electricForceWeight = DockArgs.electricForceWeight,
        hydrogenBondsForceWeight = DockArgs.hydrogenBondsForceWeight,
        bondForceWeight = DockArgs.bondForceWeight
      )

      case _ => sys.error(usage)
    }
    new MultipleInitialsDocker(innerDocker, DockArgs.initAngle, DockArgs.initConfigLevel, DockArgs.randomInit)
  }

  /**
    * Gets a list of commands from file viewinit.txt.
    * These commands are then run each time the JMOL view is reset.
    * If the file does not exist or cannot be read, returns an empty list.
    */
  def getViewInitCmds = {
    try {
      scala.io.Source.fromFile("viewinit.txt").getLines.toSeq
    } catch {
      case e: Exception => List[String]()
    }
  }

  def parseArgs(args: Array[String]) = {
    var i = 0

    try {
      while (i < args.length) {
        args(i).trim match {
          case "" => i += 1
          case "-a" => DockArgs.pathA = args(i + 1);i += 2
          case "-b" => DockArgs.pathB = args(i + 1);i += 2
          case "-out" => DockArgs.pathOut = args(i + 1); i += 2
          case "-dir" => DockArgs.dir = args(i+1); i += 2
          case "-ref" => DockArgs.pathRefs = args(i+1); i += 2
          case "--randominit" => DockArgs.randomInit = true; i += 1
          case "-initangle" => DockArgs.initAngle = Math.toRadians(args(i + 1).toDouble); i += 2
          case "-initlevel" => DockArgs.initConfigLevel = args(i + 1).toInt; i += 2
          case "-docker" => DockArgs.dockerName = args(i + 1); i += 2
          case "--consolelog" => DockArgs.consoleLog = true; i += 1
          case "-surface" => DockArgs.surface = args(i + 1).toDouble ; i += 2

          // Hill Climbing:
          case "-maxiters" => DockArgs.maxIters = args(i+1).toInt; i+=2
          case "-scorer" => DockArgs.scorerName = args(i + 1); i += 2

          // FF:
          case "-threshold" => DockArgs.threshold = args(i + 1).toDouble; i += 2
          case "-permeability" => DockArgs.permeability = args(i + 1).toDouble ; i += 2
          case "--ignoreAhydrogens" => DockArgs.ignoreAHydrogens = true; i += 1
          case "--nogui" => DockArgs.liveGui = false; i += 1
          case "-balance" =>
            val balanceStrs = args(i+1).split(',')
            if (balanceStrs.length != 4)
              sys.error(usage)

            val atomicForceWeight = balanceStrs(0).toDouble
            val electricForceWeight = balanceStrs(1).toDouble
            val hydrogenBondsForceWeight = balanceStrs(2).toDouble
            val bondForceWeight = balanceStrs(3).toDouble
            val sum = atomicForceWeight + electricForceWeight + hydrogenBondsForceWeight + bondForceWeight

            if (atomicForceWeight < 0 || electricForceWeight < 0 ||
              hydrogenBondsForceWeight < 0 || bondForceWeight < 0 || sum==0.0)
              sys.error(usage)

            DockArgs.geometricForceWeight = atomicForceWeight / sum
            DockArgs.electricForceWeight = electricForceWeight / sum
            DockArgs.hydrogenBondsForceWeight = hydrogenBondsForceWeight / sum
            DockArgs.bondForceWeight = bondForceWeight / sum

            i += 2




          case _ => sys.error(usage)
        }
      }
    }
    catch { case _:Exception => sys.error(usage) }

    if (!DockArgs.valid)
      sys.error(usage)
  }

  object DockArgs {
    var dir = ""
    var pathA = ""
    var pathB = ""
    var pathOut = ""
    var pathRefs = ""

    var dockerName = ""
    var scorerName = ""
    var consoleLog = false
    var liveGui = true
    var initAngle = Math.toRadians(90)
    var initConfigLevel = 0
    var randomInit = false

    // Hill Climbing specifict
    var maxIters = Int.MaxValue

    // Force vector docker specifics:
    var surface = 1.4
    var permeability = 0.5
    var ignoreAHydrogens = false
    var threshold = 1.0e-5
      // Force vector force balance: either all 3 set, or all 3 with default values
    var geometricForceWeight = 0.25
    var electricForceWeight = 0.25
    var hydrogenBondsForceWeight = 0.25
    var bondForceWeight = 0.25

    def valid = pathA != "" && pathB != "" && pathOut != "" && dockerName != "" &&
      Math.abs(geometricForceWeight + electricForceWeight + hydrogenBondsForceWeight + bondForceWeight - 1.0) < 1.0e-5 &&
      geometricForceWeight >= 0 && electricForceWeight >= 0 && hydrogenBondsForceWeight >= 0 && bondForceWeight >= 0 &&
      surface >= 0 &&
      permeability >= 0 && permeability <= 1

    def fullPathA = dir + pathA
    def fullPathB = dir + pathB
    def fullPathOut = dir + pathOut
    def fullPathsRef : Seq[String] =
      if (pathRefs == "")
        List[String]()
      else
        for (p <- pathRefs.split(",")) yield dir + p.trim

  }

}
