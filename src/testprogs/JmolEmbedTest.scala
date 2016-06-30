package testprogs

import io.threadcso._

import docking._
import docking.dockscore._
import docking.docksearch._
import io.{PdbReader, XyzWriter}
import jmolint.{JmolCmds, JmolFrame, JmolPanel}
import model.Molecule
import opt.Action

import JmolCmds._

object JmolEmbedTest {

  val frame = new JmolFrame(500, 500, true)
  val jmolPanel = frame.getPanel

  def showActions(chan: ?[Any], panel: JmolPanel, scorer: Scorer) = proc {
    repeat {
      val l = chan?
      val s = l match {
        case a: Action => val cmd = JmolCmds.cmd(a); panel.execute(cmd); cmd
        case d: DockingState => s"score: ${scorer.score(d)}"
        case other => other.toString
      }

      //println(s)
      //sleep(Sec/10)
      //readLine()
      ()
    }
  }

  def main(args: Array[String]) {
    System.out.println("loading")
    val molA: Molecule = new PdbReader(args(0)).read
    val molB: Molecule = new PdbReader(args(1)).read
    val scorer: Scorer = new SurfaceDistanceScorer(1.4)
    val docker: Docker = new ForceVectorDocker(1.4)
    //val docker: Docker = AtomPairDocker

    jmolPanel.openAndColor((args(0), "gray"), (args(1), "red"))
    jmolPanel.execute(
      setLog(0),
      zoom(50),
      save
    )

    val chan = OneOne[Any]
    var docked = null.asInstanceOf[DockingState]
    (proc { docked = docker.dock(molA, molB, scorer, chan); chan.close } ||
      showActions(chan, jmolPanel, scorer))()

    new XyzWriter(args(2)).write(docked.b)      // write docked b to file
    jmolPanel.openAndColor((args(0), "gray"), (args(2), "red"))  // show original a and modified b



  }
}
