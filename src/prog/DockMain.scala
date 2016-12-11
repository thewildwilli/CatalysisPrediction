package prog

import docking.dockscore._
import docking.docksearch._
import docking.docksearch.forcevector.{DockingParamsHeuristic, ForceVectorDocker, ForceVectorScore}
import io.Mol2Writer
import io.threadcso._
import jmolint.JmolCmds._
import jmolint.{JmolCmds, JmolFrame, JmolMoleculeReader, JmolPanel}
import model.Molecule
import opt.Action
import profiling.Profiler


object DockMain {
  val usage = "USAGE: scala DockMain -a (path to A) " +
    "-b (path to B) " +
    "-out (B's output path, mol2 format) " +
    " [-dir d] [--randominit] [-workers w] [-ref a,b,c,...]" +
    "-docker (ehc|forcevector) [--consolelog] [--nogui] " +
    " [-initangle a] [initlevel i]" +
    " [-maxiters m]" +
    " [--ignoreAhydrogens] " +
    " [-surface s] " +
    " [-permeability p]" +
    " [-balance atomic,electric,hydrogenbond,bondstrength] " +
    " [-threshold t] "

  val frame = new JmolFrame(500, 500, false)
  val jmolPanel = frame.getPanel

  def main(args: Array[String]): Unit = {
    parseArgs(args)
    val (docked, rmsd, _) = doMainDock(DockArgs.fullPathA, DockArgs.fullPathB, DockArgs.fullPathOut)

    new Mol2Writer(DockArgs.fullPathOut).write(docked)      // write docked b to file
    jmolPanel.openFiles(List(DockArgs.fullPathA, DockArgs.fullPathOut))  // show original a and modified b
    jmolPanel.execSeq(DockArgs.viewInitCmds)
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
    val docker = getDocker(molA, molB)

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
      if (DockArgs.liveGui && DockArgs.workers <= 1) {
        val s = msg match
        {
          case a: Action => val cmd = JmolCmds.cmd(a); panel.exec(cmd); cmd
          case "save" => panel.exec(JmolCmds.save); "save"
          case "reset" => panel.exec(JmolCmds.reset); panel.execSeq(DockArgs.viewInitCmds); "reset"
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
  def getDocker(molA: Molecule, molB: Molecule) = {
    if (DockArgs.workers <= 1 )
        new MultipleInitialsDocker(getInnerDocker(molA, molB), DockArgs.initAngle,
          DockArgs.initConfigLevel, DockArgs.randomInit)
    else
      new MultipleInitialsConcurrentDocker(() => getInnerDocker(molA, molB), DockArgs.initAngle,
        DockArgs.initConfigLevel, DockArgs.randomInit, DockArgs.workers)
  }

  private def getInnerDocker(molA: Molecule, molB: Molecule) = {
    val ffparams = getFFParams(molA, molB)
    val scorer =
      if (DockArgs.scorerName == "ff")
        new ForceVectorScore(ffparams, 1.25, DockArgs.dockerName == "chain")
      else
        new SurfaceDistanceScorer(DockArgs.surface)

    DockArgs.dockerName match {
      case "atompair" => new AtomPairDocker(new SurfaceDistanceScorer(DockArgs.surface))
      case "ehc" => new EhcDocker(scorer, DockArgs.maxIters)
      case "forcevector" => new ForceVectorDocker(ffparams)
      case "forcevectorc" => new ForceVectorConcurrentDockerDocker(ffparams)
      case "chain" => new FFandEHC(ffparams, scorer, DockArgs.maxIters
      )
      case _ => sys.error(usage)
    }
  }

  def getFFParams(molA: Molecule, molB: Molecule) = {
    val params = DockingParamsHeuristic.estimate(molA, molB)
    if (DockArgs.balanceIsSet){
      params.geometricForceWeight = DockArgs.geometricForceWeight
      params.electricForceWeight = DockArgs.electricForceWeight
      params.hydrogenBondsForceWeight = DockArgs.hydrogenBondsForceWeight
      params.bondForceWeight = DockArgs.bondForceWeight
    }
    if (DockArgs.surfaceIsSet)
      params.surface = DockArgs.surface
    if (DockArgs.permeabilityIsSet)
      params.permeability = DockArgs.permeability
    if (DockArgs.ignoreHydrogensIsSet)
      params.ignoreAHydrogens = DockArgs.ignoreAHydrogens
    params
  }

  /**
    * Gets a list of commands from file viewinit.txt.
    * These commands are then run each time the JMOL view is reset.
    * If the file does not exist or cannot be read, returns an empty list.
    */
  def getViewInitCmds = {
    try {
      scala.io.Source.fromFile(DockArgs.dir + "viewinit.txt").getLines.toSeq
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
          case "-workers" => DockArgs.workers = args(i+1).toInt; i += 2; Profiler.setWorkers(DockArgs.workers + 8)
          case "-ref" => DockArgs.pathRefs = args(i+1); i += 2
          case "--randominit" => DockArgs.randomInit = true; i += 1
          case "-initangle" => DockArgs.initAngle = Math.toRadians(args(i + 1).toDouble); i += 2
          case "-initlevel" => DockArgs.initConfigLevel = args(i + 1).toInt; i += 2
          case "-docker" => DockArgs.dockerName = args(i + 1); i += 2
          case "--consolelog" => DockArgs.consoleLog = true; i += 1
          case "-surface" => DockArgs.surface = args(i + 1).toDouble; DockArgs.surfaceIsSet = true; i += 2

          // Hill Climbing:
          case "-maxiters" => DockArgs.maxIters = args(i+1).toInt; i+=2
          case "-scorer" => DockArgs.scorerName = args(i + 1); i += 2

          // FF:
          case "-threshold" => DockArgs.threshold = args(i + 1).toDouble; i += 2
          case "-permeability" => DockArgs.permeability = args(i + 1).toDouble; DockArgs.permeabilityIsSet = true; i += 2
          case "--ignoreAhydrogens" => DockArgs.ignoreAHydrogens = true; DockArgs.ignoreHydrogensIsSet = true; i += 1
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

            DockArgs.balanceIsSet = true
            i += 2
          case _ => sys.error(usage)
        }
      }
    }
    catch { case _:Exception => sys.error(usage) }

    if (!DockArgs.valid)
      sys.error(usage)

    DockArgs.viewInitCmds = getViewInitCmds
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
    var workers = 1

    var viewInitCmds = Seq[String]()

    // Hill Climbing specific
    var maxIters = Int.MaxValue - 1

    // Force vector docker specifics:
    var surface = 1.4
    var surfaceIsSet = false  // becomes true if overridden in program args

    var permeability = 0.5
    var permeabilityIsSet = false  // becomes true if overridden in program args

    var ignoreAHydrogens = false
    var ignoreHydrogensIsSet = false  // becomes true if overridden in program args

    var threshold = 1.0e-5
      // Force vector force balance: either all 3 set, or all 3 with default values
    var geometricForceWeight = 0.25
    var electricForceWeight = 0.25
    var hydrogenBondsForceWeight = 0.25
    var bondForceWeight = 0.25
    var balanceIsSet = false  // becomes true if overridden in program args

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
