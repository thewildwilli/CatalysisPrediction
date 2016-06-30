package testprogs

import io.threadcso._
import java.awt.Dimension
import javax.swing.JFrame

import docking._
import docking.dockscore._
import docking.docksearch._
import io.{Pdb2DReader, XyzWriter}
import jmolint.JmolPanel
import model.Molecule
import opt.Action
import jmolint.JmolCommand

object JmolEmbedTest {

  val frame = new JFrame()
  frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE)
  val contentPane = frame.getContentPane
  val jmolPanel = new JmolPanel()

  jmolPanel.setPreferredSize(new Dimension(500,500))
  contentPane.add(jmolPanel)

  frame.pack()
  frame.setVisible(true)
  frame.setAlwaysOnTop(true)

  def printActions(chan: ?[Any], panel: JmolPanel, scorer: Scorer) = proc {
    var sbs = false
    repeat {
      val l = chan?
      val s = l match {
        case a: Action => val cmd = JmolCommand.cmd(a); panel.executeCmd(cmd); cmd
        case d: DockingState => s"score: ${scorer.score(d)}"
        case other => other.toString
      }

      println(s)
      //sleep(Sec/10)
      //readLine()
      ()
    }
  }

  def main(args: Array[String]) {
    System.out.println("loading")
    val molA: Molecule = new Pdb2DReader(args(0)).read
    val molB: Molecule = new Pdb2DReader(args(1)).read
    val scorer: Scorer = new SurfaceDistanceScorer(0)
    val docker: Docker = new ForceVectorDocker(1.4)
    //val docker: Docker = AtomPairDocker

    jmolPanel.getViewer.openFiles(Array(args(0), args(1)))
    jmolPanel.executeCmd(
      "set logLevel 0",

      "select model=1.1",
      "color gray",

      "select model=2.1",
      "color red",

      "model all",
      "select model=2.1",
      "zoom 50",
      "save state"
    )

    val chan = OneOne[Any]
    var docked = null.asInstanceOf[DockingState]
    (proc { docked = docker.dock(molA, molB, scorer, chan); chan.close } ||
      printActions(chan, jmolPanel, scorer))()

    val result: Molecule = molA.clone
    result.setElement("O")
    for (bAtom <- docked.b.Atoms)
      result.JAtoms.add(bAtom)

    new XyzWriter(args(2)).write(result)

    jmolPanel.getViewer.openFiles(Array(args(2)))

    println("done")

  }
}
