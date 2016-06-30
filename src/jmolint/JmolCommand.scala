package jmolint

import docking._
import opt.Action

object JmolCommand {
  def cmd(action: Action): String = {
    action match {
      case t: Translate => s"translate {${t.v(0)} ${t.v(1)} ${t.v(2)} }"

      case r: Rotate =>
        // http://chemapps.stolaf.edu/jmol/docs/#rotate. {atom expression or point} {atom expression or point}
        // Example: rotate selected {0 100 0} {1 100 0} 10: rotate about point {0 100 0} with axis {1 0 0} 10 degrees
        val end = r.c + r.axis
        s"rotate selected {${r.c(0)} ${r.c(1)} ${r.c(2)}} {${end(0)} ${end(1)} ${end(2)}} ${Math.toDegrees(r.angRad)} "

      case Reset => "restore state"
    }
  }
}