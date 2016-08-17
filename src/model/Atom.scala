//Created by Ernesto on 23/05/2016.
package model
import breeze.linalg.{DenseMatrix, DenseVector, norm}
import scala.collection.mutable.ArrayBuffer

/** Distances in Angstrongs */
class Atom (val id: Int, elem: String, initX: Double, initY: Double, initZ: Double,
            var partialCharge: Double = 0.0,
            val atomName: String = "",
            var substructureId: String = "",
            var substructureName: String = "",
            var bonds: List[Int] = List(),
            var isSurface: Boolean = false) {

  private var _element = "C"
  private var _radius = 1.4
  var coords = DenseVector(initX, initY, initZ)

  setElement(elem)

  def element = _element
  def setElement(e: String) = {_element = e; _radius = VanDerWaalsRadii(_element)}
  def isElement(e: String) = _element.toUpperCase == e.toUpperCase
  def isOneOf(elems: String*) = elems.exists(e => isElement(e))
  def radius = _radius

  def x = coords(0)
  def y = coords(1)
  def z = coords(2)

  def addBond(atomIndex: Int) = { bonds ::= atomIndex }

  def distTo(other: Atom): Double = distTo(other.coords)
  def distTo(point: DenseVector[Double]): Double = norm(this.coords - point) //  Math.sqrt( Math.pow(this.x - point(0), 2) + Math.pow(this.y - point(1), 2) + Math.pow(this.z - point(2),2))

  def translate(v: DenseVector[Double]): Unit ={
    coords.+=(v)
  }

  def transform(m: DenseMatrix[Double]): Unit = {
    val updated = m * DenseVector(x, y, z, 1.0)
    coords = updated(0 to 2)
  }

  override def clone = {
    new Atom(this.id, elem, this.x, this.y, this.z, this.partialCharge, this.atomName,
      this.substructureId, this.substructureName, this.bonds, this.isSurface)
  }

  override def toString: String = s"Atom $id $element at $x, $y, $z with charge $partialCharge"
}
