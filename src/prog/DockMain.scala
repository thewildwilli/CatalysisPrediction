package prog

import docking.{Docker, DockLog}
import docking.dockscore._
import docking.docksearch._
import docking.docksearch.forcevector.{DockingParamsHeuristic, ForceVectorDocker, ForceVectorScore}
import docking.initials.{GlobeInitialsGenerator, RandomInitials}
import io.Mol2Writer
import io.threadcso._
import jmolint.JmolCmds._
import jmolint.{JmolCmds, JmolFrame, JmolMoleculeReader, JmolPanel}
import model.{Geometry, Molecule, Rotate}
import opt.Action
import profiling.Profiler


object DockMain {
  val usage = "USAGE: scala DockMain -a (path to A) " +
    "-b (path to B) " +
    "-out (B's output path, mol2 format) " +
    " [-dir d] [--randominit] [-workers w] [-ref a,b,c,...]" +
    "-docker (ehc|forcevector) [--consolelog] [--nogui] " +
    " [-initials (globe {angle} {level} ¦ (random {number})]"
    " [-maxiters m]" +
    " [--ignoreAhydrogens] " +
    " [-surface s] " +
    " [-permeability p]" +
    " [-balance atomic,electric,hydrogenbond,bondstrength] " +
    " [-threshold t] " +
    " [--pdbAddHydrogens] [--exit]"

  def main(args: Array[String]): Unit = {
    val dockArgs = parseArgs(args)
    val jmolPanel = new JmolFrame(500, 500, false, dockArgs.liveGuiEnabled).getPanel
    val (docked, (closestRef, rmsd), score) = doMainDock(jmolPanel, dockArgs)

    new Mol2Writer(dockArgs.fullPathOut).write(docked)      // write docked b to file
    jmolPanel.openFiles(List(dockArgs.fullPathA, dockArgs.fullPathOut) ++ dockArgs.fullPathsRef)
    jmolPanel.execSeq(dockArgs.viewInitCmds)
    println(s"Finished with RMSD: $rmsd ($closestRef), score: $score")
    Profiler.report
    if (dockArgs.exit)
      sys.exit(0)
  }

  def doMainDock(jmolPanel: JmolPanel, dockArgs: DockArgs) = {
    jmolPanel.openFiles(List(dockArgs.fullPathA, dockArgs.fullPathB), dockArgs.pdbAddHydrogens)
    jmolPanel.execSync(
      selectModel("2.1"),
      setLog(0),
      zoom(50),
      save
    )

    val molA: Molecule = JmolMoleculeReader.read(jmolPanel, 0)
    val molB: Molecule = JmolMoleculeReader.read(jmolPanel, 1)
    val docker = getDocker(molA, molB, dockArgs)
    val dockResult = getDockResult(jmolPanel, docker, molA, molB, dockArgs)
    val docked = dockResult._1
    val score = dockResult._2
    (docked, getRMSD(jmolPanel, docked, dockArgs), score)
  }

  def getDockResult(jmolPanel: JmolPanel, docker: Docker, molA: Molecule, molB: Molecule, dockArgs: DockArgs) = {
    if (dockArgs.liveGuiEnabled)
      getDockResultWithLiveUI(jmolPanel, docker, molA, molB, dockArgs)
    else {
      if (dockArgs.randomInit)
        doRandomRotation(molB)
      docker.dock(molA, molB.clone, DockLog.dummy)
    }
  }

  def getDockResultWithLiveUI(jmolPanel: JmolPanel, docker: Docker, molA: Molecule, molB: Molecule, dockArgs: DockArgs) = {
    val chan = OneOneBuf[Any](5000)
    val dockLog = new DockLog(chan, enabled = dockArgs.workers == 1)

    if (dockArgs.randomInit)
      doRandomRotation(molB, chan)

    var dockResult = (null.asInstanceOf[Molecule], 0.0)
    (proc { dockResult = docker.dock(molA, molB.clone, dockLog); chan.closeOut } ||
      showActions(chan, jmolPanel, dockArgs))()
    dockResult
  }

  def showActions(chan: ?[Any], panel: JmolPanel, dockArgs: DockArgs) = proc {
    repeat {
      val msg = chan?()
      val s = msg match
      {
        case a: Action => val cmds = JmolCmds.cmds(a); panel.exec(cmds:_*); cmds
        case "save" => panel.exec(JmolCmds.save); "save"
        case "reset" => panel.exec(JmolCmds.reset); panel.execSeq(dockArgs.viewInitCmds); "reset"
        case other => other.toString
      }
      if (dockArgs.consoleLog)
        println(s)
      sleepms(100)
    }
    chan.closeIn
  }

  def getRMSD(jmolPanel: JmolPanel, docked: Molecule, dockArgs: DockArgs) = {
    jmolPanel.openFiles(dockArgs.fullPathsRef, dockArgs.pdbAddHydrogens)
    (for (i <- dockArgs.fullPathsRef.indices) yield {
      val refMol = JmolMoleculeReader.read(jmolPanel, i)
      (dockArgs.fullPathsRef(i), docked.rmsd(refMol))
    }).minBy(_._2)
  }


  /**
    * Gets an instance of a docker given the name. In the future, this could be
    * enhanced to use reflection.
    */
  def getDocker(molA: Molecule, molB: Molecule, dockArgs: DockArgs) = {
      new MultipleInitialsConcurrentDocker(() => getInnerDocker(molA, molB, dockArgs),
        getInitialsGenerator(dockArgs), dockArgs.workers)
  }

  private def getInnerDocker(molA: Molecule, molB: Molecule, dockArgs: DockArgs) = {
    val ffparams = getFFParams(molA, molB, dockArgs)
    val scorer =
      if (dockArgs.scorerName == "ff")
        new ForceVectorScore(ffparams, 1.25, dockArgs.dockerName == "chain")
      else
        new SurfaceDistanceScorer(dockArgs.surface)

    dockArgs.dockerName match {
      case "atompair" => new AtomPairDocker(new SurfaceDistanceScorer(dockArgs.surface))
      case "ehc" => new EhcDocker(scorer, dockArgs.maxIters)
      case "forcevector" => new ForceVectorDocker(ffparams)
      case "chain" => new FFandEHC(ffparams, scorer, dockArgs.maxIters
      )
      case _ => sys.error(usage)
    }
  }

  private def getInitialsGenerator(dockArgs: DockArgs) =
    if (dockArgs.initials == "globe")
      new GlobeInitialsGenerator(dockArgs.initConfigLevel, dockArgs.initAngle)
    else if (dockArgs.initials == "random")
      new RandomInitials(dockArgs.initNumber)
    else ???

  def getFFParams(molA: Molecule, molB: Molecule, dockArgs: DockArgs) = {
    val params = DockingParamsHeuristic.estimate(molA, molB)
    if (dockArgs.balanceIsSet){
      params.geometricForceWeight = dockArgs.geometricForceWeight
      params.electricForceWeight = dockArgs.electricForceWeight
      params.hydrogenBondsForceWeight = dockArgs.hydrogenBondsForceWeight
      params.bondForceWeight = dockArgs.bondForceWeight
    }
    if (dockArgs.surfaceIsSet)
      params.surface = dockArgs.surface
    if (dockArgs.permeabilityIsSet)
      params.permeability = dockArgs.permeability
    if (dockArgs.ignoreHydrogensIsSet)
      params.ignoreAHydrogens = dockArgs.ignoreAHydrogens
    if (dockArgs.thresholdIsSet)
      params.threshold = dockArgs.threshold
    params
  }

  /**
    * Gets a list of commands from file viewinit.txt.
    * These commands are then run each time the JMOL view is reset.
    * If the file does not exist or cannot be read, returns an empty list.
    */
  def getViewInitCmds(dockArgs: DockArgs) = {
    try {
      scala.io.Source.fromFile(dockArgs.dir + "viewinit.txt").getLines.toSeq
    } catch {
      case e: Exception => List[String]()
    }
  }

  def parseArgs(args: Array[String]) = {
    val dockArgs = new DockArgs
    var i = 0

    try {
      while (i < args.length) {
        args(i).trim match {
          case "" => i += 1
          case "-a" => dockArgs.pathA = args(i + 1);i += 2
          case "-b" => dockArgs.pathB = args(i + 1);i += 2
          case "-out" => dockArgs.pathOut = args(i + 1); i += 2
          case "-dir" => dockArgs.dir = args(i+1); i += 2
          case "-workers" => dockArgs.workers = args(i+1).toInt; i += 2; Profiler.setWorkers(dockArgs.workers + 8)
          case "-ref" => dockArgs.pathRefs = args(i+1); i += 2
          case "--randominit" => dockArgs.randomInit = true; i += 1
          case "-docker" => dockArgs.dockerName = args(i + 1); i += 2
          case "--consolelog" => dockArgs.consoleLog = true; i += 1
          case "-surface" => dockArgs.surface = args(i + 1).toDouble; dockArgs.surfaceIsSet = true; i += 2
          case "--pdbAddHydrogens" => dockArgs.pdbAddHydrogens = true; i += 1

          // Initials:
          case "-initials" =>
            dockArgs.initials = args(i + 1)
            dockArgs.initials match {
              case "globe" =>
                dockArgs.initAngle = Math.toRadians(args(i + 2).toDouble)
                dockArgs.initConfigLevel = args(i + 3).toInt
                i += 4
              case "random" =>
                dockArgs.initNumber = args(i+2).toInt
                i += 3
              case unknown => sys.error(s"Unknown initials type: $unknown, $usage")
            }

          // Hill Climbing:
          case "-maxiters" => dockArgs.maxIters = args(i+1).toInt; i+=2
          case "-scorer" => dockArgs.scorerName = args(i + 1); i += 2

          // FF:
          case "-threshold" => dockArgs.threshold = args(i + 1).toDouble; dockArgs.thresholdIsSet = true; i += 2
          case "-permeability" => dockArgs.permeability = args(i + 1).toDouble; dockArgs.permeabilityIsSet = true; i += 2
          case "--ignoreAhydrogens" => dockArgs.ignoreAHydrogens = true; dockArgs.ignoreHydrogensIsSet = true; i += 1
          case "--nogui" => dockArgs.liveGuiSelected = false; i += 1
          case "--exit" => dockArgs.exitSelected = true; i += 1
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

            dockArgs.geometricForceWeight = atomicForceWeight / sum
            dockArgs.electricForceWeight = electricForceWeight / sum
            dockArgs.hydrogenBondsForceWeight = hydrogenBondsForceWeight / sum
            dockArgs.bondForceWeight = bondForceWeight / sum

            dockArgs.balanceIsSet = true
            i += 2
          case unknown => sys.error(s"Unknown argument: $unknown, $usage")
        }
      }
    }
    catch { case e:Exception => sys.error(s"${e.getMessage}, $usage") }

    if (!dockArgs.valid)
      sys.error(usage)

    dockArgs.viewInitCmds = getViewInitCmds(dockArgs)
    dockArgs
  }

  private def doRandomRotation(m: Molecule, log: ![Any] = null): Unit = {
    val (axis, angle) = Geometry.randomAxisAndAngle
    m.rotate(m.getGeometricCentre, axis, angle)
    if (log != null)
      log! new Rotate(m.getGeometricCentre, axis, angle)
  }

  class DockArgs {
    var dir = ""
    var pathA = ""
    var pathB = ""
    var pathOut = ""
    var pathRefs = ""

    var dockerName = ""
    var scorerName = ""
    var consoleLog = false
    var liveGuiSelected = true

    var initials = "globe"  // either globe or random
    var initAngle = Math.toRadians(90.0)  // only for globe
    var initConfigLevel = 0 // only for globe
    var initNumber = 10  // only for random


    var randomInit = false  // shuffle B molecule at the start
    var workers = 1

    var viewInitCmds = Seq[String]()
    var pdbAddHydrogens = false

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
    var thresholdIsSet = false

      // Force vector force balance: either all 3 set, or all 3 with default values
    var geometricForceWeight = 0.25
    var electricForceWeight = 0.25
    var hydrogenBondsForceWeight = 0.25
    var bondForceWeight = 0.25
    var balanceIsSet = false  // becomes true if overridden in program args

    var exitSelected = false // exit after dock

    def liveGuiEnabled = !this.exitSelected && this.liveGuiSelected && this.workers == 1
    def exit = !this.liveGuiEnabled

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
        List[String](fullPathB)
      else
        for (p <- pathRefs.split(",")) yield dir + p.trim
  }
}
