//Created by Ernesto on 23/05/2016.
package model
import breeze.linalg._

/** Distances in Angstrongs */
class Atom (val id: Int, val element: String, initX: Double, initY: Double, initZ: Double,
            var partialCharge: Double = 0.0,
            val atomName: String = "",
            var substructureId: String = "",
            var substructureName: String = "",
            var bonds: List[Int] = List(),
            var isSurface: Boolean = false) {

  var coords = DenseVector(initX, initY, initZ)

  def isElement(e: String) = element.toUpperCase == e.toUpperCase
  def isOneOf(elems: String*) = elems.exists(e => isElement(e))
  lazy val radius = VanDerWaalsRadii(element)
  lazy val isH = isElement("H")  // cache

  def x = coords(0)
  def y = coords(1)
  def z = coords(2)

  def addBond(atomIndex: Int) = { bonds ::= atomIndex }

  def distTo(other: Atom): Double = distTo(other.coords)
  def distTo(point: DenseVector[Double]): Double =
    // this is twice as fast as norm(this.coords - point) :
    Math.sqrt( Math.pow(this.x - point(0), 2) + Math.pow(this.y - point(1), 2) + Math.pow(this.z - point(2),2))

  /** Gets the distance to other, a vector from this to other, and a vector in the
    * same direction with norm 1 */
  def distDifDir(other: Atom) = {
    val diffX = other.x - x
    val diffY = other.y - y
    val diffZ = other.z - z
    val dist = Math.sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ)
    (dist, DenseVector(diffX, diffY, diffZ), DenseVector(diffX/dist, diffY/dist, diffZ/dist))
  }

  def translate(v: DenseVector[Double]): Unit ={
    coords.+=(v)
  }

  def transform(m: DenseMatrix[Double]): Unit = {
    val updated = m * DenseVector(x, y, z, 1.0)
    coords = updated(0 to 2)
  }

  override def clone = {
    new Atom(this.id, element, this.x, this.y, this.z, this.partialCharge, this.atomName,
      this.substructureId, this.substructureName, this.bonds, this.isSurface)
  }

  override def toString: String = s"Atom $id $element at $x, $y, $z with charge $partialCharge"
}
