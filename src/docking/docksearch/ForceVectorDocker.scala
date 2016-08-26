package docking.docksearch

import breeze.linalg
import breeze.linalg._
import docking.dockscore.Scorer
import docking._
import opt._
import io.threadcso._
import model._
import profiling.Profiler

class ForceVectorDocker(val surface: Double = 1.4,
                        val permeability: Double = 0.5,
                        val maxDecelerations: Int = 10,
                        val ignoreAHydrogens: Boolean = false,
                        val threshold: Double = 1.0e-5,
                        val geometricForceWeight: Double = .25,
                        val electricForceWeight: Double = .25,
                        val hydrogenBondsForceWeight: Double = .25,
                        val bondForceWeight: Double = .25) extends Docker {

  val hBondDistance = 1.8

  val initialDeltaAngle = Math.toRadians(5)
  val initialDeltaSpace = 0.5
  val minDeltaSpace = 1.0           // minimum to be used only when the molecules are too far apart
  val minCoverage = 1.25            // how many times the optimal distance must always be considered, at least
  var maxCoverage = Double.PositiveInfinity
  var atomWithinMinCover = false

  var avgBondEnergy = 0.0
  var startingPos = 0
  var softness = 1.0

  var maxAngle = initialDeltaAngle
  var maxTranslate = initialDeltaSpace

  var currWindowScore = Double.NegativeInfinity
  var lastWindowScore = Double.NegativeInfinity
  var bestWindowScore = Double.NegativeInfinity
  var currIterScore = 0.0
  var lastIterScore = Double.NegativeInfinity

  var done = false
  var decelerations = 0
  val window = 5                  // Assess score every window iterations
  var approachPhase = true
  var shortestDistance = Double.NegativeInfinity
  var forceShortestDistance = Double.PositiveInfinity

  /**  Docks b into a from b's initial position and orientation, using force vectors.
    *   Try to minimize score
    */
  override def dock(molA: Molecule, molB: Molecule, log: ![Any]) = {
    maxAngle = initialDeltaAngle
    maxTranslate = initialDeltaSpace
    currWindowScore = 0.0
    lastWindowScore = Double.NegativeInfinity
    bestWindowScore = Double.NegativeInfinity
    currIterScore = 0.0
    lastIterScore = Double.NegativeInfinity
    done = false
    decelerations = 0
    softness = 1.0
    approachPhase = true
    atomWithinMinCover = false
    shortestDistance = Double.NegativeInfinity
    forceShortestDistance = Double.PositiveInfinity
    avgBondEnergy =
      (for (a <- molA.atoms; b <- molB.atoms) yield BondEnergy(a.element, b.element)).sum /
        (molA.atoms.size * molB.atoms.size)

    //println(s"Docking from pos $startingPos")
    startingPos+=1

    // move the molecule to the starting position


    var i = 1
    while (!done){
      val forces = Profiler.time("getForces") { getForces(molA, molB) }             // it is a list of pairs (atomB, force)

      val (translation, _) = Profiler.time("getTranslation") { getTranslation(forces) }
      Profiler.time("translate") { molB.translate(translation) }
      log!new Translate(translation)

      val (axis, angle) = Profiler.time("getRotation") {  getRotation(molB.getGeometricCentre, forces) }
      Profiler.time("rotate") {  molB.rotate(molB.getGeometricCentre, axis, angle)}
      log!new Rotate(molB.getGeometricCentre, axis, angle)

      Profiler.time("updateState") {   updateState(i, molA, molB)       }  // updates done and all other algorithm control variables
      i+=1

      //println(s"$currIterScore")
    }
    (molB, currIterScore)
  }

  private def updateState(i: Integer, molA: Molecule, molB: Molecule)  = {
    lastIterScore = currIterScore
    currIterScore = getScore(molA, molB, !approachPhase)

    if (currIterScore < lastIterScore) {
      maxAngle = maxAngle * 0.95
      maxTranslate = maxTranslate * 0.95
      softness = softness * 0.5
    }

    currWindowScore += currIterScore
    if (i % window == 0) {

      if (currWindowScore - lastWindowScore < threshold) {
        if (decelerations >= maxDecelerations  && !approachPhase) {
          done = true
        } else {
          // decelerate
          maxAngle = maxAngle * 0.75
          maxTranslate = maxTranslate * 0.75
          softness = Math.max(softness - 0.2, 0)
          decelerations += 1

          if (approachPhase && softness <= 0 && atomWithinMinCover){
            // end approach phase
            approachPhase = false
            currWindowScore = 0.0
            bestWindowScore = Double.NegativeInfinity
            decelerations = 0
            currIterScore = 0
          }
        }
      } else if (currWindowScore > bestWindowScore + threshold) {
        decelerations = 0
        if (approachPhase) {
          maxAngle = initialDeltaAngle
          maxTranslate = initialDeltaSpace
        }
      }
      //println(s"$i, score: $currWindowScore, soft: $softness, decel $decelerations, approach: $approachPhase, maxT: $maxTranslate")
      lastWindowScore = currWindowScore
      bestWindowScore = Math.max(bestWindowScore, currWindowScore)
      currWindowScore = 0
    }

  }


  /** Returns a pair (vector, boolean). The vector has the translation,
    * the boolean indicates if minimum translation had to be applied
    */
  private def getTranslation(forces: Iterable[(Atom, DenseVector[Double])]) = {
    val forcesOnly = forces.map(pair => pair._2)
    val netForce = forcesOnly.foldLeft(DenseVector(0.0, 0.0, 0.0))((a, b) => a + b)
    val totalForceAmount = forcesOnly.map(f => norm(f)).sum
    val minApplied = (totalForceAmount < minDeltaSpace) && approachPhase
    val translateDistance =
      if (minApplied)
        minDeltaSpace
      else
        Math.min(norm(netForce), maxTranslate)
    val translation =
      if ( norm(netForce) > 0)
        (netForce * translateDistance) / norm(netForce)
      else
        DenseVector(0.0, 0.0, 0.0)
    (translation, minApplied)
  }

  private def getRotation(centre: DenseVector[Double], forces: Iterable[(Atom, DenseVector[Double])]) = {
    // torque in Euler axis/angle format is cross(r, force) - see http://web.mit.edu/8.01t/www/materials/modules/chapter21.pdf
    val torques = forces.map { case (atomB, force) =>
      val r = atomB.coords - centre // radius vector from centre to atom
      linalg.cross(r, force)
    }
    val netTorque = torques.foldLeft(DenseVector(0.0, 0.0, 0.0))((a, b) => a + b)     // add all torques together
    val netTorqueNorm = norm(netTorque)
    val axis = if (netTorqueNorm > 0) netTorque / netTorqueNorm else DenseVector(0.0, 0.0, 0.0)

    // calculate moment of inertia with all atoms assumed weight 0.1 otherwise angles are too small:
    var moi = 0.0
    for( (atom, force ) <- forces) {
      val distToAxis = Geometry.distToLine(centre, centre + netTorque, atom.coords)
      moi += 0.1 * distToAxis*distToAxis
    }

    val angle = if (moi > 0) Math.min(netTorqueNorm / moi, maxAngle) else 0.0
    //println(Math.toDegrees(angle))
    (axis, angle)
  }

  /** Returns a list of pairs (Atom, DenseVector) containing the net force for
    * each atom in molB.     */
  private def getForces(molA: Molecule, molB: Molecule) = {
    val forces = for (atomB <- molB.surfaceAtoms) yield {
      var forceOnAtomB = DenseVector(0.0, 0.0, 0.0)
      for (atomA <- molA.surfaceAtoms(ignoreAHydrogens)) {
        val actualDistance = atomA.distTo(atomB)
        if (approachPhase || actualDistance <= optimalDistance(atomA, atomB) * minCoverage) {
          forceOnAtomB += atomToAtomForce(atomA, atomB, molA, molB)
          atomWithinMinCover |= actualDistance <= optimalDistance(atomA, atomB) * minCoverage;
        }
      }

      for (atomA <- molA.innerAtoms(ignoreAHydrogens))
        forceOnAtomB += (1-permeability) * getPenetrationPenaltyForce(atomA, atomB)

      (atomB, forceOnAtomB)
    }

    forces.filter(t => norm(t._2) > 0)
  }

  private def inLineOfSight(atomA: Atom, atomB: Atom, molA: Molecule, molB: Molecule) = {
    // can't check against all atoms - too expensive, and actually we don't really want that
    val allAtoms = molA.atomsBoundTo(atomA) ++ molB.atomsBoundTo(atomB)
    val actualDistance = atomA.distTo(atomB)
    !allAtoms.exists(a => (a ne atomA) && (a ne atomB)
      && !a.isElement("H")
      && a.distTo(atomA) < actualDistance
      && a.distTo(atomB) < actualDistance
      && Geometry.distToLine(atomA.coords, atomB.coords, a.coords) < a.radius * 0.8
    )
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
  private def molToAtomForce(molA: Molecule, atomB: Atom, molB: Molecule) = {
    molA.surfaceAtoms(ignoreAHydrogens)
      .map(atomA => atomToAtomForce(atomA, atomB, molA, molB))
      .reduce((a,b) => a+b)
  }

  /** calculates the force that atomA exerts on atomB
    * This force has 2 components:
    *   - Atomic attraction: atoms naturally attract each other so as to dock.
    *     They attract up to the optimal distance, and if close, they reject.
    *     This is a sort of simulated gravity.
    *   - Electric force: different charges attract, equal charges reject.
    *     Equal forces reject if the distance is closer than the optimal distance.
    *     Different forces attract if the distance is farther than the optimal distance.
    * */
  private def atomToAtomForce(atomA: Atom, atomB: Atom,
                              molA: Molecule, molB: Molecule): DenseVector[Double] = {
    val dif = atomA.coords - atomB.coords;                    // direction from b to a
    val actualDistance = norm(dif)
    val dir = dif / actualDistance;                                // normalized to length 1
    val optimal = optimalDistance(atomA, atomB)

    val geometricForce = dir * getGeometricForceNorm(atomA, atomB, actualDistance, optimal)
    val electricForce = dir * getElectricForceNorm(atomA, atomB, actualDistance, optimal)
    val hydrogenBondsForce = getHydrogenBondsForce(atomA, atomB, molA, molB, actualDistance)
    val bondForce = dir * getBondForceNorm(atomA, atomB, actualDistance, optimal)

    val force =     // weighted result:
      geometricForce * geometricForceWeight +
        electricForce * electricForceWeight +
        hydrogenBondsForce * hydrogenBondsForceWeight +
        bondForce * bondForceWeight

    force
  }


  /** Calculates the atomic force norm f such that the force that a exerts on b
    * is f * a normalized vector pointing from b to a.
    * Positive = attraction, negative = repulsion */
  private def getGeometricForceNorm(atomA: Atom, atomB: Atom, actualDistance: Double, optimal: Double) = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1
    if (atomA.isElement("H") || atomB.isElement("H"))
      0.0
    else {
      val normalized = actualDistance / optimal
      val force = explog(normalized / (softness * 3 + 1))
      force
    }
  }

  private def getElectricForceNorm(atomA: Atom, atomB: Atom, actualDistance: Double,
                                   optimal: Double) = {
    val chargeProduct = atomA.partialCharge*atomB.partialCharge
    //-chargeProduct / (actualDistance*actualDistance) // Coulomb

     if (chargeProduct < 0)
       - explog(actualDistance / optimal) * chargeProduct
     else if (chargeProduct > 0)
       minusExpOverX(actualDistance / optimal) * chargeProduct
     else
      0.0
  }

  private def getHydrogenBondsForce(atomA: Atom, atomB: Atom, molA: Molecule,
                                   molB: Molecule, actualDistance: Double) = {

    def getForce(vectToHBondSpot: DenseVector[Double]) =
      if (vectToHBondSpot != null) {
        val dist = norm(vectToHBondSpot)
        val force = - (xexp(dist) / xexp.maxY) * atomA.partialCharge * atomB.partialCharge
        val dir = vectToHBondSpot / dist
        force * dir
      } else
        DenseVector(0.0, 0.0, 0.0)

    val forceOnB = getForce(getVectToHBondSpot(atomA, atomB, molA))
    val forceOnA = getForce(getVectToHBondSpot(atomB, atomA, molB))
    forceOnB + (-1.0 * forceOnA)
  }

  private def getVectToHBondSpot(h: Atom, other: Atom, hmol: Molecule) = {
    if (canHBond(h, other)){
      val hNeighbour = hmol(h.bonds(0));
      if (hNeighbour.partialCharge < 0 && hNeighbour.isOneOf("F", "O", "N")) {
        val dirToNeighbour = hNeighbour.coords - h.coords
        val toSpot = -(dirToNeighbour / norm(dirToNeighbour)) * hBondDistance // 1.8 Angtrom in the opposite direction
        val spot = h.coords + toSpot
        spot - other.coords
      } else
        null
    }
    else
      null
  }

  private def canHBond(h: Atom, other: Atom) =
    h.isElement("H") && h.partialCharge > 0 && other.partialCharge < 0 &&
      other.isOneOf("F", "O", "N") && h.bonds.size == 1

  private def getBondForceNorm(atomA: Atom, atomB: Atom, actualDistance: Double, optimal: Double) = {
    val bondEnergy = BondEnergy(atomA.element, atomB.element)
    val normalized = actualDistance/optimal
    val frac = bondEnergy / avgBondEnergy
    if (frac > 1)
      explog(normalized) * frac
    else if (frac < 1)
      minusExpOverX(normalized) * (1-frac)
    else
      0.0
  }

  private def getPenetrationPenaltyForce(atomA: Atom, atomB: Atom) = {
    val actualDistance = atomA.distTo(atomB)
    val penalizationDistance = (atomA.radius + atomB.radius) * 0.9
    if ( (!atomA.isSurface || !atomB.isSurface) &&
      !atomA.isElement("H") && !atomB.isElement("H") &&
      actualDistance < penalizationDistance) {
      val dif = atomA.coords - atomB.coords;                         // direction from b to a
      val dir = dif / actualDistance;                                // normalized to length 1
      val normalized = actualDistance/penalizationDistance
      dir * Math.log(normalized)
    } else
      DenseVector(0.0, 0.0, 0.0)
  }

  /* --- Scoring --- */
  /**
    * Molecules assumed to be nonempty!
    *
    * @param molA
    * @param molB
    * @param onlyTargetRadius
    * @return
    */
  private def getScore(molA: Molecule, molB: Molecule, onlyTargetRadius: Boolean) = {
    var totalGeometricScore = 0.0
    var totalElectricScore = 0.0
    var totalHBondScore = 0.0
    var totalBondStrengthScore = 0.0

    var aCount = 0
    var bCount = 0
    for (a <- molA.surfaceAtoms(ignoreAHydrogens) ){
      aCount += 1
      for(b <- molB.surfaceAtoms) {
        bCount += 1
        val actualDistance = a.distTo(b)
        if (!onlyTargetRadius || actualDistance <= optimalDistance(a, b) * minCoverage) {
          val n = actualDistance / optimalDistance(a, b)
          if (!a.isElement("H") && !b.isElement("H"))
            totalGeometricScore += explog2(n * explog2.maxX) / explog2.maxY

          //TODO: I changed something here..
          val chargeProduct = a.partialCharge*b.partialCharge
          totalElectricScore += (
            if (chargeProduct < 0)
              (explog2(n * explog2.maxX) / explog2.maxY) * a.partialCharge * b.partialCharge * -1.0
            else if (chargeProduct > 0 && n <= 2)
              minusExpOverX(n) * chargeProduct
            else
              0.0
            )

          def hScore(vectToSpot: DenseVector[Double]) = {
            if (vectToSpot != null)
              gaussian(norm(vectToSpot)) //* a.partialCharge * b.partialCharge * -1.0
            else
              0.0
          }
          totalHBondScore += hScore(getVectToHBondSpot(a, b, molA)) + hScore(getVectToHBondSpot(b, a, molB))

          val bondEnergy = BondEnergy(a.element, b.element)
          val frac = bondEnergy / avgBondEnergy
          totalBondStrengthScore += (
            if (frac > 1)
              (explog2(n * explog2.maxX) / explog2.maxY) * frac
            else if (frac < 1)
              minusExpOverX(n) * (1-frac)
            else
              0.0
            )
        }
      }
    }
    val totalScore =
      totalGeometricScore * geometricForceWeight +
        totalElectricScore * electricForceWeight +
        totalHBondScore * hydrogenBondsForceWeight +
        totalBondStrengthScore * bondForceWeight
    //println(s"SCORES: geo: ${totalGeometricScore * geometricForceWeight}, electric: ${totalElectricScore * electricForceWeight}, hbond: ${totalHBondScore * hydrogenBondsForceWeight}, bondstrength: ${totalBondStrengthScore * bondForceWeight}")
    totalScore / (Math.min(aCount, bCount))
  }





  /* --- Distance Functions --- */

  private def optimalDistance(a: Atom, b: Atom) = {
    if (canHBond(a, b) || canHBond(b, a))
      hBondDistance + 2 * surface
    else
      (a.radius + b.radius) * 0.9 + 2 * surface
        //(if (a.isElement("H") || b.isElement("H")) 0.0 else 2 * surface)
  }

  private def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)

  private def explog(x: Double) = Math.exp(-Math.pow(x,2))*Math.log(x)    // Max of explog is reached at 1.327864011995167

  private def minusExpOverX(x: Double) = - Math.exp(-x) * 0.1 / x

  /** exp(-(x-1) ** 2)*log(x)
    * peaks at 1.6290055996317214 with value 0.3285225256677019
    * */
  private object explog2 {
    def apply(x: Double) = Math.exp(-Math.pow(x - 1, 2)) * Math.log(x)
    val maxX = 1.6290055996317214
    val maxY = 0.3285225256677019
  }

  private object xexp {
    def apply(x: Double) = x * Math.exp(-x)
    val maxX = 1.0
    val maxY = 0.367879441171442
  }

  private def gaussian(x: Double) = Math.exp(-Math.pow(x, 2))


}
