// Created by Ernesto on 23/05/2016.

package model

import breeze.linalg.{DenseVector, DenseMatrix}

import collection.JavaConverters._

class Molecule(val Atoms: scala.collection.mutable.Buffer[Atom]) {

  def this(l: Seq[Atom]) {this (l.toBuffer[Atom])}
  def this(){ this (new scala.collection.mutable.ArrayBuffer[Atom]())}


  def apply(i: Int) = Atoms(i)

  def translate(v: DenseVector[Double]): Unit = {
    for (a <- this.Atoms)
      a.translate(v)
  }

  /** Rotates the molecule with respect to X axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateX(centre: DenseVector[Double], angRad: Double) {
    // Translate the centre to the origin, perform rotation around the origin, then translate back
    translate(centre * -1.0)
    transform (DenseMatrix(
      (1.0, 0.0,                0.0),
      (0.0, Math.cos(angRad), -Math.sin(angRad)),
      (0.0, Math.sin(angRad), Math.cos(angRad))))
    translate(centre)
  }

  /** Rotates the molecule with respect to Y axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateY(centre: DenseVector[Double], angRad: Double) {
    translate(centre * -1.0)
    transform (DenseMatrix(
      (Math.cos(angRad) , 0.0, Math.sin(angRad)),
      (0.0              , 1.0, 0.0),
      (-Math.sin(angRad), 0.0, Math.cos(angRad))))
    translate(centre)
  }

  /** Rotates the molecule with respect to Z axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateZ(centre: DenseVector[Double], angRad: Double) {
    translate(centre * -1.0)
    transform (DenseMatrix(
      (Math.cos(angRad) , -Math.sin(angRad), 0.0),
      (Math.sin(angRad) , Math.cos(angRad) , 0.0),
      (0.0              , 0.0              , 1.0)))
    translate(centre)
  }

  def transform(m: DenseMatrix[Double]) = {
    for (a <- Atoms)
      a.transform(m)
  }

  def getGeometricCentre = {
    // Add all coords vectors together and then divide by the number of atoms. Can be optimized caching.
    (Atoms.map(a => a.coords).reduce(_+_))./(Atoms.size + 0.0)
  }

  /** The molecule radius is the distance from the geometricCentre to its farthest atom */
  def getRadius = {
    val centre = getGeometricCentre
    Atoms.map(a => a.distTo(centre)).max
  }

  def computeSurfaceAtoms3D = computeSurfaceAtoms(6)
  def computeSurfaceAtoms2D = computeSurfaceAtoms(4)
  private def computeSurfaceAtoms(maxNeighbours: Int) = {
    for (a <- Atoms) {
      a.isSurface = Atoms.filter(b => (b ne a) && a.distTo(b) <= 2.16).size <= maxNeighbours //ne = reference inequality. 216 is C=C bond length.
    }
  }


  /** Deep copy */
  override def clone = new Molecule(this.Atoms.map(a => a.clone))

  def setElement(e: String) = {
    for (a <- this.Atoms)
      a.setElement(e)
  }

  def JAtoms = Atoms.asJava


}
