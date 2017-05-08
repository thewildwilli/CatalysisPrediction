package docking.docksearch.forcevector

import model.Atom
import HBond._

object Functions {
  def optimalDistance(a: Atom, b: Atom, surface: Double) = {
    if (canHBond(a, b) || canHBond(b, a))
      hBondDistance + 2 * surface
    else
      (a.radius + b.radius) * 0.9 + 2 * surface
    //(if (a.isH || b.isH) 0.0 else 2 * surface)
  }

  object expsquare {
    // Max of expsquare is reached at x=sqrt(2), y=0.13533528323661267. The root is 1.0.
    def apply(x: Double) = Math.exp(-Math.pow(x, 2)) * (Math.pow(x, 2) - 1)
    val maxX = Math.sqrt(2.0)
    val maxY = 0.13533528323661267
  }

  object explog {
    // Max of explog is reached at x=1.327864011995167, y=0.048630065613819703
    def apply(x: Double) = Math.exp(-Math.pow(x, 2)) * Math.log(x)
    val maxX = 1.327864011995167
    val maxY = 0.048630065613819703
  }

  def minusExpOverX(x: Double) = - Math.exp(-x) * 0.1 / x

  /** exp(-(x-1) ** 2)*log(x)
    * peaks at 1.6290055996317214 with value 0.3285225256677019
    * */
  object explog2 {
    def apply(x: Double) = Math.exp(-Math.pow(x - 1, 2)) * Math.log(x)
    val maxX = 1.6290055996317214
    val maxY = 0.3285225256677019
  }

  object xexp {
    def apply(x: Double) = x * Math.exp(-x)
    val maxX = 1.0
    val maxY = 0.367879441171442
  }

  def gaussian(x: Double) = Math.exp(-Math.pow(x, 2))

}
