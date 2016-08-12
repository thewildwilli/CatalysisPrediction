package prog

import docking.dockscore._
import docking.docksearch._
import io.Mol2Writer
import io.threadcso._
import jmolint.JmolCmds._
import jmolint.{JmolCmds, JmolFrame, JmolMoleculeReader, JmolPanel}
import model.Molecule
import opt.Action


object DockMain {
  val usage = "USAGE: scala DockMain -a (path to A, pdb format) " +
    "-b (path to B, pdb format) " +
    "-out (B's output path, xyz format) " +
    "-docker (atompair|forcevector) [--consolelog] " +
    " [--ignorehydrogen] " +
    " [-balance atomic,electric,bond] "

  val frame = new JmolFrame(500, 500, false)
  val jmolPanel = frame.getPanel

  val viewInitCmds = getViewInitCmds

  def main(args: Array[String]): Unit = {
    val startTime = System.currentTimeMillis()
    parseArgs(args)

    jmolPanel.openAndColor((DockArgs.pathA, "gray"), (DockArgs.pathB, "red"))
    jmolPanel.exec(
      setLog(0),
      zoom(50),
      save
    )

    val molA: Molecule = JmolMoleculeReader.read(jmolPanel, 0)
    val molB: Molecule = JmolMoleculeReader.read(jmolPanel, 1)
    val docker = getDocker

    val chan = OneOneBuf[Any](5)
    var dockResult = (null.asInstanceOf[Molecule], 0.0)
    (proc { dockResult = docker.dock(molA, molB, chan); chan.close } ||
      showActions(chan, jmolPanel))()

    val docked = dockResult._1
    val score = dockResult._2

    new Mol2Writer(DockArgs.pathOut).write(docked)      // write docked b to file
    jmolPanel.openAndColor((DockArgs.pathA, "gray"), (DockArgs.pathOut, "red"))  // show original a and modified b
    jmolPanel.execSeq(viewInitCmds)

    println(s"Finished with score: $score, total time: ${System.currentTimeMillis()-startTime}ms")
  }

  def showActions(chan: ?[Any], panel: JmolPanel) = proc {
    repeat {
      val s = chan? match {
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

  /**
    * Gets an instance of a docker given the name. In the future, this could be
    * enhanced to use reflection.
    */
  def getDocker = DockArgs.dockerName match {
      case "atompair" => new AtomPairDocker(new SurfaceDistanceScorer(1.4))

      case "forcevector" => new ForceVectorDocker(
        surface = 1.4,
        maxDecelerations = 10,
        ignoreHydrogen = DockArgs.ignoreHydrogen,
        threshold = DockArgs.threshold,
        atomicForceWeight = DockArgs.atomicForceWeight,
        electricForceWeight = DockArgs.electricForceWeight,
        bondForceWeight = DockArgs.bondForceWeight
      )

      case _ => sys.error(usage)
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
        args(i) match {
          case "-a" => DockArgs.pathA = args(i + 1);i += 2
          case "-b" => DockArgs.pathB = args(i + 1);i += 2
          case "-out" => DockArgs.pathOut = args(i + 1); i += 2
          case "-docker" => DockArgs.dockerName = args(i + 1); i += 2
          case "-scorer" => DockArgs.scorerName = args(i + 1); i += 2
          case "--consolelog" => DockArgs.consoleLog = true; i += 1
          case "--ignorehydrogen" => DockArgs.ignoreHydrogen = true; i += 1
          case "-balance" =>
            val balanceStrs = args(i+1).split(',')
            if (balanceStrs.length != 3)
              sys.error(usage)

            val atomicForceWeight = balanceStrs(0).toDouble
            val electricForceWeight = balanceStrs(1).toDouble
            val bondForceWeight = balanceStrs(2).toDouble
            val sum = atomicForceWeight + electricForceWeight + bondForceWeight

            if (atomicForceWeight < 0 || electricForceWeight < 0 || bondForceWeight < 0 || sum==0.0)
              sys.error(usage)

            DockArgs.atomicForceWeight = atomicForceWeight / sum
            DockArgs.electricForceWeight = electricForceWeight / sum
            DockArgs.bondForceWeight = bondForceWeight / sum

            i += 2

          case _ => sys.error(usage)
        }
      }
    }
    catch { case _:Exception => sys.error(usage) }

    if (!DockArgs.valid)
      sys.error(usage)

    print(DockArgs)
  }

  object DockArgs {
    var pathA = ""
    var pathB = ""
    var pathOut = ""
    var dockerName = ""
    var scorerName = ""
    var consoleLog = false

    // Force vector docker specifics:
    var ignoreHydrogen = false
    var threshold = 1.0e-5
      // Force vector force balance: either all 3 set, or all 3 with default values
    var atomicForceWeight = 1.0
    var electricForceWeight = 1.0
    var bondForceWeight = 1.0

    def valid = pathA != "" && pathB != "" && pathOut != "" && dockerName != "" &&
      Math.abs(atomicForceWeight + electricForceWeight + bondForceWeight - 1.0) < 1.0e-5
  }

}
