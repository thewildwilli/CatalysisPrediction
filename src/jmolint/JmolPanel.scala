package jmolint

import javax.swing._
import java.awt._

import org.jmol.adapter.smarter.SmarterJmolAdapter
import org.jmol.api._
import org.jmol.viewer.Viewer

// Created by Ernesto on 30/06/2016.
class JmolPanel extends JPanel {
  private val serialVersionUID: Long = -3661941083797644242L

  val adapter = new SmarterJmolAdapter
  val viewer = JmolViewer.allocateViewer(this, adapter).asInstanceOf[Viewer]
  private val currentSize: Dimension = new Dimension
  private val rectClip: Rectangle = new Rectangle

  override def paint(g: Graphics) {
    getSize(currentSize)
    g.getClipBounds(rectClip)
    viewer.renderScreenImage(g, currentSize.width, currentSize.height)
  }

  // CONVENIENCE METHODS:
  def exec(cmds: String*) {
    for (cmd <- cmds) viewer.evalString(cmd)
  }

  def execSeq(cmds: Seq[String]): Unit = {
    for (cmd <- cmds) viewer.evalString(cmd)
  }

  def execSync(cmds: String*) {
    for (cmd <- cmds) viewer.evalStringQuietSync(cmd, false, true)
    while (viewer.getScriptQueueInfo()) Thread.sleep(10)
  }

  def openFiles(paths: Seq[String]) {
    execSync(JmolCmds.loadFiles(paths))
    exec(JmolCmds.showAllModels)
  }


  def openAndColor(pathsAndColors: (String, String)*) {
    val paths = pathsAndColors.map(pair => pair._1)
    openFiles(paths)
    val colours = pathsAndColors.map(pair => pair._2)
    for (i <- colours.indices) {
      exec(JmolCmds.modelColor(s"${i+1}.1", colours(i)))

    }

  }


}
