// Created by Ernesto on 23/05/2016.

package model

import breeze.linalg.{DenseVector, DenseMatrix}

import collection.JavaConverters._

class Molecule(val Atoms: scala.collection.mutable.Set[Atom]) {

  def this(){ this (scala.collection.mutable.Set[Atom]()) }

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

  /** Deep copy */
  override def clone = {
    val result = new Molecule()
    for (a <- this.Atoms)
      result.Atoms.add(a.clone())
    result
  }

  def setElement(e: Char) = {
    for (a <- this.Atoms)
      a.element = e
  }

  def JAtoms = Atoms.asJava


}
