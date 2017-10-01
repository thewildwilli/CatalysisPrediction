package docking.docksearch.forcevector

import breeze.linalg._
import Functions._
import model.{Atom, BondEnergy, Molecule}

/**
  * Each force for the ForceVectorDocker must provide, given a pair
  * of atoms, a force vector. On calling the endIter method, it must
  * return the total accumulated score for all the pairs seen since the
  * last call to endIter. The idea is that ForceVectorDocker needs to do
  * only one pass of iterating all the atom pairs (as opposed to doing
  * 2 passes: one to get force vectors and the other to get scores).
  */
trait Force {
  val weight: Double
  protected var score = 0.0

  /** The apply method is called to get the force vector from an
    * atom of B to an atom of A, or None if it does not apply.
    */
  def apply(molA: Molecule, molB: Molecule, atomA: Atom,
            atomB: Atom, actualDist: Double,
            dir: DenseVector[Double]
           ): Option[DenseVector[Double]]

  def totalScore = score
  def weightedScore = score * weight

  /** Signals that the current atom pair iteration has ended. The force
    * returns the total score for all the pairs seen since the last call
    * to endIter
    */
  def endIter(): Double = try { score } finally { score = 0.0 }
}

class GeometricForce(val weight: Double, surface: Double) extends Force {
  override def apply(molA: Molecule, molB: Molecule, atomA: Atom, atomB: Atom, actualDist: Double,
                     dir: DenseVector[Double]): Option[DenseVector[Double]] = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1
      val optimal = atomA.radius + atomB.radius + 2 * surface
      val normalized = actualDist / optimal
      score += expsquare(normalized * expsquare.maxX) / expsquare.maxY
      val force = expsquare(normalized)
      Some(dir * force)
  }
}

class ElectricForce(val weight: Double, surface: Double) extends Force {
  override def apply(molA: Molecule, molB: Molecule, atomA: Atom, atomB: Atom, actualDist: Double,
                     dir: DenseVector[Double]): Option[DenseVector[Double]] = {
    val chargeProduct = atomA.partialCharge * atomB.partialCharge
    val optimal = optimalDistance(atomA, atomB, surface)
    val normalized = actualDist / optimal

    if (chargeProduct < 0) {
      if (atomA.isSurface && atomB.isSurface)
        score += (explog2(normalized * explog2.maxX) / explog2.maxY) * chargeProduct * -1.0

      val forceNorm = - explog(normalized) * chargeProduct
      Some(dir * forceNorm)
    } else if (chargeProduct > 0 && normalized <= 2) {  // atoms farther than twice the optimal don't count
      if (atomA.isSurface && atomB.isSurface)
        score += minusExpOverX(normalized) * chargeProduct

      val forceNorm = minusExpOverX(normalized) * chargeProduct
      Some(dir * forceNorm)
    } else None
  }
}

class HydrogenBondForce(val weight: Double, ignoreAHydrogens: Boolean) extends Force {

  private def getForce(vectToHBondSpot: DenseVector[Double], atomA: Atom, atomB: Atom) =
    if (vectToHBondSpot != null) {
      val dist = norm(vectToHBondSpot)
      val force = - (xexp(dist) / xexp.maxY) * atomA.partialCharge * atomB.partialCharge
      val dir = vectToHBondSpot / dist
      Some(force * dir)
    } else
      None

  override def apply(molA: Molecule, molB: Molecule, atomA: Atom, atomB: Atom, actualDist: Double,
                     dir: DenseVector[Double]): Option[DenseVector[Double]] = {
    // if A hydrogens are ignored, at least get reverse forces
    if (ignoreAHydrogens) {
      val vectToSpot = HydrogenBondForce.getVectToHBondSpot(atomB, atomA, molB)
      if (vectToSpot != null && atomA.isSurface && atomB.isSurface)
        score += gaussian(norm(vectToSpot), 100.0)

      val aToBForceVector = getForce(vectToSpot, atomB, atomA)
      if (aToBForceVector.isDefined)
        Some(-1.0 * aToBForceVector.get)
      else
        None
    }else {
      val vectToSpot = HydrogenBondForce.getVectToHBondSpot(atomA, atomB, molA)
      if (vectToSpot != null && atomA.isSurface && atomB.isSurface)
        score += gaussian(norm(vectToSpot))

      getForce(vectToSpot, atomA, atomB)
    }
  }
}

object HydrogenBondForce {
  val hBondDistance = 1.8
  def getVectToHBondSpot(h: Atom, other: Atom, hmol: Molecule) = {
    if (canHBond(h, other)){
      val hNeighbour = hmol(h.bonds(0));
      if (hNeighbour.partialCharge < 0 && hNeighbour.isOneOf("F", "O", "N")) {
        val (_, _, dirToNeighbour) = h.distDifDir(hNeighbour)
        val toSpot = -dirToNeighbour * hBondDistance // 1.8 Angstroms in the opposite direction
        val spot = h.coords + toSpot
        spot - other.coords
      } else
        null
    }
    else
      null
  }

  def canHBond(h: Atom, other: Atom) =
    h.isH && h.partialCharge > 0 && other.partialCharge < 0 &&
      h.bonds.size == 1 && other.isOneOf("F", "O", "N")
}

class BondForce(val weight: Double, molA: Molecule, molB: Molecule, surface: Double) extends Force() {
  val avgBondEnergy = getAvgBondEnergy(molA, molB)

  def getAvgBondEnergy(molA: Molecule, molB: Molecule) =
    (for (a <- molA.atoms; b <- molB.atoms) yield BondEnergy(a.element, b.element)).sum /
      (molA.atoms.size * molB.atoms.size)

  override def apply(molA: Molecule, molB: Molecule, atomA: Atom, atomB: Atom, actualDist: Double,
                     dir: DenseVector[Double]): Option[DenseVector[Double]] = {
    val bondEnergy = BondEnergy(atomA.element, atomB.element)
    val normalized = actualDist/optimalDistance(atomA, atomB, surface)
    val frac = bondEnergy / avgBondEnergy
    if (frac > 1) {
      score += (explog2(normalized * explog2.maxX) / explog2.maxY) * frac
      Some(explog(normalized) * frac * dir)
    }else if (frac < 1 && normalized <= 2) {
      score += minusExpOverX(normalized) * (1 - frac)
      Some(minusExpOverX(normalized) * (1 - frac) * dir)
    }else
      None
  }
}
