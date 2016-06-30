package jmolint

import javax.swing._
import java.awt._

import org.jmol.adapter.smarter.SmarterJmolAdapter
import org.jmol.api.{JmolAdapter, JmolSimpleViewer}

// Created by Ernesto on 30/06/2016.
class JmolPanel extends JPanel {
  private val serialVersionUID: Long = -3661941083797644242L

  private val adapter = new SmarterJmolAdapter
  private val viewer = JmolSimpleViewer.allocateSimpleViewer(this, adapter)
  private val currentSize: Dimension = new Dimension
  private val rectClip: Rectangle = new Rectangle

  def getViewer = viewer

  override def paint(g: Graphics) {
    getSize(currentSize)
    g.getClipBounds(rectClip)
    viewer.renderScreenImage(g, currentSize.width, currentSize.height)
  }

  // CONVENIENCE METHODS:
  def execute(cmds: String*) {
    for (cmd <- cmds) viewer.evalString(cmd)
  }

  def openFiles(paths: Array[String]): Unit ={
    viewer.openFiles(paths)
    execute(JmolCmds.showAllModels)
  }

  def openFiles(paths: String*) {
    openFiles(paths.toArray)
  }

  def openAndColor(pathsAndColors: (String, String)*) {
    val paths = pathsAndColors.map(pair => pair._1)
    openFiles(paths.toArray)
    val colours = pathsAndColors.map(pair => pair._2)
    for (i <- colours.indices) {
      execute(JmolCmds.modelColor(s"${i+1}.1", colours(i)))

    }

  }


}
