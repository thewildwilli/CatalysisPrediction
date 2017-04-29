package jmolint

import model._
import opt.Action

object JmolCmds {
  def cmds(action: Action): Seq[String] = {
    action match {
      case t: Translate => List(s"translate {${t.v(0)} ${t.v(1)} ${t.v(2)} }")

      case r: Rotate =>
        // http://chemapps.stolaf.edu/jmol/docs/#rotate. {atom expression or point} {atom expression or point}
        // Example: rotate selected {0 100 0} {1 100 0} 10: rotate about point {0 100 0} with axis {1 0 0} 10 degrees
        val end = r.c + r.axis
        List(s"rotate selected {${r.c(0)} ${r.c(1)} ${r.c(2)}} {${end(0)} ${end(1)} ${end(2)}} ${Math.toDegrees(r.angRad)} ")

      case multi: MultiTransform =>
        (for (t <- multi.transforms) yield cmds(t)).flatten
    }
  }

  def setLog(level: Int) = s"set logLevel $level"
  def selectModel(model: String) = s"select model=$model"
  def setColor(c:String) = s"color $c"
  def modelColor(model:String, color:String) = selectModel(model) + "; " + setColor(color)
  def showAllModels = "model all"
  def zoom(percent: Int) = s"zoom $percent"
  def save = "save state"
  def reset = "restore state"
  def loadFiles(paths: Seq[String]) = "load files " + paths.map(p => "\"" + p.replace("\\", "/") + "\"").mkString(" ")


}