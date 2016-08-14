// Created by Ernesto on 23/05/2016.

package model

import breeze.linalg.{DenseMatrix, DenseVector}

import collection.JavaConverters._

class Molecule(var atomMap: Map[Int, Atom]) {

  def this(){ this (Map[Int, Atom]()) }
  def this(l: Iterable[Atom]) { this(); for (a <- l) this.atomMap += (a.id -> a) }

  def apply(i: Int) = atomMap(i)
  def Atoms = atomMap.values

  def translate(v: DenseVector[Double]): Unit = {
    for (a <- this.Atoms)
      a.translate(v)
    if (_geometricCentre != None)
      _geometricCentre = Some(_geometricCentre.get + v)
  }

  def rotate(centre: DenseVector[Double], axis: DenseVector[Double],
             angRad: Double) = {
    translate(centre * -1.0)
    transform(Geometry.rotate(axis, angRad))
    translate(centre)
  }

  /** Rotates the molecule with respect to X axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates. Returns this same object.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateX(centre: DenseVector[Double], angRad: Double) = {
    // Translate the centre to the origin, perform rotation around the origin, then translate back
    translate(centre * -1.0)
    transform (DenseMatrix(
      (1.0, 0.0,                0.0),
      (0.0, Math.cos(angRad), -Math.sin(angRad)),
      (0.0, Math.sin(angRad), Math.cos(angRad))))
    translate(centre)
    this
  }

  /** Rotates the molecule with respect to Y axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates. Returns this same object.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateY(centre: DenseVector[Double], angRad: Double) = {
    translate(centre * -1.0)
    transform (DenseMatrix(
      (Math.cos(angRad) , 0.0, Math.sin(angRad)),
      (0.0              , 1.0, 0.0),
      (-Math.sin(angRad), 0.0, Math.cos(angRad))))
    translate(centre)
    this
  }

  /** Rotates the molecule with respect to Z axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates. Returns this same object.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateZ(centre: DenseVector[Double], angRad: Double) = {
    translate(centre * -1.0)
    transform (DenseMatrix(
      (Math.cos(angRad) , -Math.sin(angRad), 0.0),
      (Math.sin(angRad) , Math.cos(angRad) , 0.0),
      (0.0              , 0.0              , 1.0)))
    translate(centre)
    this
  }

  def transform(m: DenseMatrix[Double]) = {
    for (a <- Atoms)
      a.transform(m)
    _geometricCentre = None         // could update the centre too
  }

  var _geometricCentre: Option[DenseVector[Double]] = None
  def getGeometricCentre = {
    // Add all coords vectors together and then divide by the number of atoms.
    if (_geometricCentre == None)
      _geometricCentre = Some(Atoms.map(a => a.coords).reduce(_+_) ./ (Atoms.size + 0.0))
    _geometricCentre.get
  }

  /** The molecule radius is the distance from the geometricCentre to its farthest atom */
  var _radius: Option[Double] = None
  def getRadius = {
    val centre = getGeometricCentre
    if (_radius == None)
      _radius = Some(Atoms.map(a => a.distTo(centre)).max)
    _radius.get
  }

  /* This is now taken from JMOL
  def computeSurfaceAtoms() = {
    for (a <- Atoms) {
      val pointsInSurface = Geometry.sphereOrientations(a.radius, Math.toRadians(45))
      a.isSurface = pointsInSurface.exists(point =>
        !Atoms.exists(other => (other ne a) && other.distTo(point) < other.radius) )
    }

  }*/

  def computeSurfaceAtoms2D() = {
    for (a <- Atoms) {
      a.isSurface = Atoms.count(b => (b ne a) && a.distTo(b) <= 2.16) <= 4 //ne = reference inequality. 216 is C=C bond length.
    }
  }



  /** Deep copy */
  override def clone = {
    val m = new Molecule(this.Atoms.map(a => a.clone));
    m._geometricCentre = this._geometricCentre;
    m._radius = this._radius
    m
  }

  def setElement(e: String) = {
    for (a <- this.Atoms)
      a.setElement(e)
  }

  /** Adds clones of all atoms from b to this. Modifies this molecule and returns itself */
  def importM(b: Molecule) = {
    for (bAtom <- b.Atoms)
      this.atomMap += (bAtom.id -> bAtom.clone)
    _geometricCentre = None
    _radius = None
    this
  }

  def JAtoms = atomMap.asJava


}