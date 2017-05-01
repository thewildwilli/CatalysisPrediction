// Created by Ernesto on 23/05/2016.

package model

import breeze.linalg._
import profiling.Profiler

import collection.JavaConverters._

class Molecule(atomSeq: Iterable[Atom]) {
  var atomMap = Map[Int, Atom]()
  for (a <- atomSeq) this.atomMap += (a.id -> a)

  // Atoms collections:
  def apply(i: Int) = atomMap(i)
  def atoms: Iterable[Atom] = atomMap.values
  private def getAtoms(ignoreHydrogen: Boolean) = atomMap.values.filter(a => !(ignoreHydrogen && a.isH))
  private lazy val nonHydrogenAtoms = getAtoms(true).toList // cache
  def atoms(ignoreHydrogen: Boolean): Iterable[Atom] = if (ignoreHydrogen) nonHydrogenAtoms else atoms

  // Surface and inner atoms:
  private def getSurfaceAtoms(ignoreHydrogen: Boolean): Iterable[Atom] = atoms.filter(a => a.isSurface && !(ignoreHydrogen && a.isH))
  private lazy val allSurfaceAtoms = getSurfaceAtoms(false).toList  // cache
  private lazy val nonHydrogenSurfaceAtoms = getSurfaceAtoms(true).toList  // cache
  def surfaceAtoms: Iterable[Atom] = surfaceAtoms(false)
  def surfaceAtoms(ignoreHydrogen: Boolean) = if (ignoreHydrogen) nonHydrogenSurfaceAtoms else allSurfaceAtoms

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
    if (_geometricCentre.isDefined) {
      // transform geometric centre as well:
      val updatedCentre = m * DenseVector.vertcat(_geometricCentre.get, DenseVector(1.0))
      _geometricCentre = Some(updatedCentre(0 to 2))
    }
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
    if (_radius.isEmpty)
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

    Math.max(this.rmsdPrime(other), other.rmsdPrime(this))
  }

  private def rmsdPrime(other: Molecule) = {
    val squaresum = (
      for {a <- this.atoms(ignoreHydrogen = this.atoms.size > 5) } yield {
        val minDist = (for {b <- other.atoms if b.element == a.element} yield a.distTo(b)).min
        minDist * minDist
      }
    ).sum
    Math.sqrt(squaresum / this.atoms.size)
  }

  private def naiveRmsd(other: Molecule) = {
    val myAtoms = this.atoms.toIndexedSeq
    val otherAtoms = other.atoms.toIndexedSeq
    val squaresum = (for (i <- myAtoms.indices) yield {
      val a = myAtoms(i)
      val b = otherAtoms(i)
      Math.pow(a.x-b.x, 2) + Math.pow(a.y-b.y, 2) + Math.pow(a.z-b.z, 2)
    }).sum
    Math.sqrt(squaresum / this.atoms.size)
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

  /** Adds clones of all atoms from b to this. Modifies this molecule and returns itself */
  def importM(b: Molecule) = {
    for (bAtom <- b.atoms)
      this.atomMap += (bAtom.id -> bAtom.clone)
    _geometricCentre = None
    _radius = None
    this
  }

  def JAtoms = atomMap.asJava
}
