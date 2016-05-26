//Created by Ernesto on 23/05/2016.
package model
import breeze.linalg.{DenseMatrix, DenseVector}

/** Distances in Armstrongs */
class Atom (initX: Double, initY: Double, initZ: Double) {

  var coords = DenseVector(initX, initY, initZ)
  var element = 'C'

  def x = coords(0)
  def y = coords(1)
  def z = coords(2)

  def distTo(other: Atom) = Math.sqrt( Math.pow(this.x - other.x, 2) + Math.pow(this.y - other.y, 2) + Math.pow(this.z - other.z,2))
  def translate(v: DenseVector[Double]): Unit ={
    coords.+=(v)
  }

  def transform(m: DenseMatrix[Double]): Unit = {
    coords = m * coords
  }

  override def clone = new Atom(this.x, this.y, this.z)
  override def toString: String = "Atom at " + x + ", " + y + ", " + z
}
