package prog

import docking.DockingState
import docking.dockscore._
import docking.docksearch._
import io.{PdbReader, XyzWriter}
import io.threadcso._
import jmolint.JmolCmds._
import jmolint.{JmolCmds, JmolFrame, JmolPanel}
import model.Molecule
import opt.Action


object DockMain {
  val usage = "USAGE: scala DockMain -a (path to A, pdb format) " +
    "-b (path to B, pdb format) " +
    "-out (B's output path, xyz format) " +
    "-docker (atompair|forcevector)"

  val frame = new JmolFrame(500, 500, true)
  val jmolPanel = frame.getPanel

  def showActions(chan: ?[Any], panel: JmolPanel, scorer: Scorer) = proc {
    repeat {
      chan? match {
        case a: Action => val cmd = JmolCmds.cmd(a); panel.execute(cmd); cmd
        case d: DockingState => s"score: ${scorer.score(d)}"
        case other => other.toString
      }
      ()
    }
  }

  def main(args: Array[String]): Unit = {
    parseArgs(args)
    val molA: Molecule = new PdbReader(DockArgs.pathA).read
    val molB: Molecule = new PdbReader(DockArgs.pathB).read
    val scorer: Scorer = new SurfaceDistanceScorer(1.4)
    val docker = getDocker

    jmolPanel.openAndColor((DockArgs.pathA, "gray"), (DockArgs.pathB, "red"))
    jmolPanel.execute(
      setLog(0),
      zoom(50),
      save
    )

    val chan = OneOne[Any]
    var docked = null.asInstanceOf[DockingState]
    (proc { docked = docker.dock(molA, molB, scorer, chan); chan.close } ||
      showActions(chan, jmolPanel, scorer))()

    new XyzWriter(DockArgs.pathOut).write(docked.b)      // write docked b to file
    jmolPanel.openAndColor((DockArgs.pathA, "gray"), (DockArgs.pathOut, "red"))  // show original a and modified b
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
          case _ => sys.error(usage)
        }
      }
    }
    catch { case _:Exception => sys.error(usage) }

    if (!DockArgs.valid)
      sys.error(usage)
  }

  /**
    * Gets an instance of a docker given the name. In the future, this could be
    * enhanced to use reflection.
    */
  def getDocker = DockArgs.dockerName match {
      case "atompair" => AtomPairDocker
      case "forcevector" => new ForceVectorDocker(1.4)
      case _ => sys.error(usage)
    }

  object DockArgs {
    var pathA = ""
    var pathB = ""
    var pathOut = ""
    var dockerName = ""
    var scorerName = ""
    def valid = pathA != "" && pathB != "" && pathOut != "" && dockerName != ""
  }

}
