package docking.docksearch.forcevector

import model.Atom
import HBond._

object Functions {
  def optimalDistance(a: Atom, b: Atom, surface: Double) = {
    if (canHBond(a, b) || canHBond(b, a))
      hBondDistance + 2 * surface
    else
      (a.radius + b.radius) * 0.9 + 2 * surface
    //(if (a.isElement("H") || b.isElement("H")) 0.0 else 2 * surface)
  }

  def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)

  def explog(x: Double) = Math.exp(-Math.pow(x,2))*Math.log(x)    // Max of explog is reached at 1.327864011995167

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
