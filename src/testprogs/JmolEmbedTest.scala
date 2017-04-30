package testprogs

import io.threadcso._
import docking._
import docking.dockscore._
import docking.docksearch._
import io.{PdbReader, XyzWriter}
import jmolint.{JmolCmds, JmolFrame, JmolMoleculeReader, JmolPanel}
import model.Molecule
import opt.Action
import JmolCmds._
import docking.docksearch.forcevector.{DockingParamsHeuristic, ForceVectorDocker}

object JmolEmbedTest {

  val frame = new JmolFrame(500, 500, true)
  val jmolPanel = frame.getPanel

  def showActions(chan: ?[Any], panel: JmolPanel, scorer: Scorer) = proc {
    repeat {
      val l = chan?
      val s = l match {
        case a: Action => val cmds = JmolCmds.cmds(a); panel.exec(cmds:_*); cmds
        case d: DockingState => s"score: ${scorer.score(d)}"
        case other => other.toString
      }
      ()
    }
  }

  def main(args: Array[String]) {
    System.out.println("loading")

    //val docker: Docker = AtomPairDocker

    jmolPanel.openAndColor((args(0), "gray"), (args(1), "red"))
    jmolPanel.exec(
      setLog(0),
      zoom(50),
      showAllModels,
      save
    )

    sleep(1*Sec)
    val molA: Molecule = JmolMoleculeReader.read(jmolPanel, 0)
    val molB: Molecule = new PdbReader(args(1)).read
    val scorer: Scorer = new SurfaceDistanceScorer(1.4)
    val docker: Docker = new ForceVectorDocker(DockingParamsHeuristic.estimate(molA, molB))

    val chan = OneOne[Any]
    var dockResult = (null.asInstanceOf[Molecule], 0.0)
    (proc { dockResult = docker.dock(molA, molB, new DockLog(chan)); chan.close } ||
      showActions(chan, jmolPanel, scorer))()
    val docked = dockResult._1
    val score = dockResult._2

    new XyzWriter(args(2)).write(docked)      // write docked b to file
    jmolPanel.openAndColor((args(0), "gray"), (args(2), "red"))  // show original a and modified b




  }
}
