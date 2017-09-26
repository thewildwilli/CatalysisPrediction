package jmolint

import java.awt.{Container, Dimension}
import javax.swing.JFrame

// Created by Ernesto on 30/06/2016.
class JmolFrame(width: Int, height: Int, alwaysOnTop: Boolean,
                visible: Boolean = true) extends JFrame {

  private val panel = new JmolPanel
  setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE)
  panel.setPreferredSize(new Dimension(width, height))
  getContentPane.add(panel)
  pack()
  setVisible(visible)
  setAlwaysOnTop(alwaysOnTop)
  setTitle("Real time docking progress")

  def getPanel: JmolPanel = panel
}
