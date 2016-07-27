package prog

import docking.DockingState
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
    " [-balance gravity,electric,bondenergy] "

  val frame = new JmolFrame(500, 500, true)
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
    val scorer: Scorer = new SurfaceDistanceScorer(1.4)
    val docker = getDocker

    val chan = OneOne[Any]
    var docked = null.asInstanceOf[DockingState]
    (proc { docked = docker.dock(molA, molB, scorer, chan); chan.close } ||
      showActions(chan, jmolPanel, scorer))()

    new Mol2Writer(DockArgs.pathOut).write(docked.b)      // write docked b to file
    jmolPanel.openAndColor((DockArgs.pathA, "gray"), (DockArgs.pathOut, "red"))  // show original a and modified b
    jmolPanel.execSeq(viewInitCmds)

    println(s"Finished with score: ${scorer.score(docked)}, total time: ${System.currentTimeMillis()-startTime}ms")
  }

  def showActions(chan: ?[Any], panel: JmolPanel, scorer: Scorer) = proc {
    repeat {
      val s = chan? match {
        case a: Action => val cmd = JmolCmds.cmd(a); panel.exec(cmd); cmd
        case d: DockingState => s"score: ${scorer.score(d)}"
        case "save" => panel.exec(JmolCmds.save); "save"
        case "reset" => panel.exec(JmolCmds.reset); panel.execSeq(viewInitCmds); "reset"
        case other => other.toString
      }
      if (DockArgs.consoleLog)
        println(s)
      ()
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
          case _ => sys.error(usage)
        }
      }
    }
    catch { case _:Exception => sys.error(usage) }

    if (!DockArgs.valid)
      sys.error(usage)

    print(DockArgs)
  }

  /**
    * Gets an instance of a docker given the name. In the future, this could be
    * enhanced to use reflection.
    */
  def getDocker = DockArgs.dockerName match {
      case "atompair" => AtomPairDocker
      case "forcevector" => new ForceVectorDocker(1.4, 10, DockArgs.ignoreHydrogen)
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

  object DockArgs {
    var pathA = ""
    var pathB = ""
    var pathOut = ""
    var dockerName = ""
    var scorerName = ""
    var consoleLog = false
    var ignoreHydrogen = false
    def valid = pathA != "" && pathB != "" && pathOut != "" && dockerName != ""
  }

}
