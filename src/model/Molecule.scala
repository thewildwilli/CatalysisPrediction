// Created by Ernesto on 23/05/2016.

package model

import breeze.linalg.{DenseMatrix, DenseVector}

import collection.JavaConverters._

class Molecule(var atomMap: Map[Int, Atom]) {

  def this(){ this (Map[Int, Atom]()) }
  def this(l: Iterable[Atom]) { this(); for (a <- l) this.atomMap += (a.id -> a) }

  def apply(i: Int) = atomMap(i)
  def atoms: Iterable[Atom] = atomMap.values
  def atoms(ignoreHydrogen: Boolean) = atomMap.values.filter(a => !(ignoreHydrogen && a.isElement("H")))
  def surfaceAtoms: Iterable[Atom] = surfaceAtoms(false)
  def surfaceAtoms(ignoreHydrogen: Boolean): Iterable[Atom] = atoms.filter(a => a.isSurface && !(ignoreHydrogen && a.isElement("H")))
  def innerAtoms(ignoreHydrogen: Boolean): Iterable[Atom] = atoms.filter(a => !a.isSurface && !(ignoreHydrogen && a.isElement("H")))

  def atomsBoundTo(a: Atom) = a.bonds.map(i => atomMap(i))

  def translate(v: DenseVector[Double]): Unit = {
    for (a <- this.atoms)
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
    rotate(centre, centre + DenseVector(1.0, 0.0, 0.0), angRad); this
  }

  /** Rotates the molecule with respect to Y axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates. Returns this same object.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateY(centre: DenseVector[Double], angRad: Double) = {
    rotate(centre, centre + DenseVector(0.0, 1.0, 0.0), angRad); this
  }

  /** Rotates the molecule with respect to Z axis around a centre, angle in radians.
    * The centre must be given as [X, Y, Z] coordinates. Returns this same object.
    * See https://en.wikipedia.org/wiki/Rotation_matrix
    */
  def rotateZ(centre: DenseVector[Double], angRad: Double) = {
    rotate(centre, centre + DenseVector(0.0, 0.0, 1.0), angRad); this
  }

  def transform(m: DenseMatrix[Double]) = {
    for (a <- atoms)
      a.transform(m)
    _geometricCentre = None         // could update the centre too
  }

  var _geometricCentre: Option[DenseVector[Double]] = None
  def getGeometricCentre = {
    // Add all coords vectors together and then divide by the number of atoms.
    if (_geometricCentre == None)
      _geometricCentre = Some(atoms.map(a => a.coords).reduce(_+_) ./ (atoms.size + 0.0))
    _geometricCentre.get
  }

  /** The molecule radius is the distance from the geometricCentre to its farthest atom */
  var _radius: Option[Double] = None
  def getRadius = {
    val centre = getGeometricCentre
    if (_radius == None)
      _radius = Some(atoms.map(a => a.distTo(centre)).max)
    _radius.get
  }

  /**
    * Computes Root-mean-square deviation with respect to another molecule.
    * All atoms are taken into account.
    * https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
    */
  def rmsd(other: Molecule) = {
    val n = this.atoms.size
    if (n != other.atoms.size) throw new Exception(s"Cannot compute RMSD of molecules of different sizes. I have ${this.atoms.size}, other has ${other.atoms.size}")
    if (n == 0) throw new Exception ("Calculating RMSD of empty molecule")

    // Atoms in the two molecules may have different ids, so they are accessed by index:
    val myAtoms = this.atoms.toIndexedSeq
    val otherAtoms = other.atoms.toIndexedSeq

    val squaresum = (for (i <- myAtoms.indices) yield {
      val a = myAtoms(i)
      val b = otherAtoms(i)
      Math.pow(a.x-b.x, 2) + Math.pow(a.y-b.y, 2) + Math.pow(a.z-b.z, 2)
    }).sum
    Math.sqrt(squaresum / n)
  }


  def computeSurfaceAtoms2D() = {
    for (a <- atoms) {
      a.isSurface = atoms.count(b => (b ne a) && a.distTo(b) <= 2.16) <= 4 //ne = reference inequality. 216 is C=C bond length.
    }
  }



  /** Deep copy */
  override def clone = {
    val m = new Molecule(this.atoms.map(a => a.clone));
    m._geometricCentre = this._geometricCentre;
    m._radius = this._radius
    m
  }

  def setElement(e: String) = {
    for (a <- this.atoms)
      a.setElement(e)
  }

  /** Adds clones of all atoms from b to this. Modifies this molecule and returns itself */
  def importM(b: Molecule) = {
    for (bAtom <- b.atoms)
      this.atomMap += (bAtom.id -> bAtom.clone)
    _geometricCentre = None
    _radius = None
    this
  }

  def JAtoms = atomMap.asJava


  /* This is now taken from JMOL
def computeSurfaceAtoms() = {
  for (a <- Atoms) {
    val pointsInSurface = Geometry.sphereOrientations(a.radius, Math.toRadians(45))
    a.isSurface = pointsInSurface.exists(point =>
      !Atoms.exists(other => (other ne a) && other.distTo(point) < other.radius) )
  }

}*/

}