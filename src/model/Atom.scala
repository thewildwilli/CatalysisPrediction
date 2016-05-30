//Created by Ernesto on 23/05/2016.
package model
import breeze.linalg.{DenseMatrix, DenseVector}

/** Distances in Angstrongs */
class Atom (elem: String, initX: Double, initY: Double, initZ: Double) {

  setElement(elem)
  var coords = DenseVector(initX, initY, initZ)
  var isSurface = false

  var _element = "C"
  var _radius = 1.4
  def element = _element
  def setElement(e: String) = {_element = e; _radius = VanDerWaalsRadii(e)}
  def radius = _radius

  def x = coords(0)
  def y = coords(1)
  def z = coords(2)

  def distTo(other: Atom): Double = distTo(other.coords)
  def distTo(point: DenseVector[Double]) = Math.sqrt( Math.pow(this.x - point(0), 2) + Math.pow(this.y - point(1), 2) + Math.pow(this.z - point(2),2))

  def translate(v: DenseVector[Double]): Unit ={
    coords.+=(v)
  }

  def transform(m: DenseMatrix[Double]): Unit = {
    coords = m * coords
  }

  override def clone = new Atom("C", this.x, this.y, this.z)
  override def toString: String = "Atom at " + x + ", " + y + ", " + z
}
