package docking.docksearch.forcevector

import breeze.linalg
import breeze.linalg._
import docking._
import model._
import profiling.Profiler
import HBond._
import Functions._

class ForceVectorDocker(val params: DockingParams) extends Docker {

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
  var bestScore = 0.0
  var bestDock: Molecule = null

  var done = false
  var decelerations = 0
  val window = 5                  // Assess score every window iterations
  var approachPhase = true
  var shortestDistance = Double.NegativeInfinity
  var forceShortestDistance = Double.PositiveInfinity

  val scorer = new ForceVectorScore(params, minCoverage)

  /**  Docks b into a from b's initial position and orientation, using force vectors.
    *   Try to minimize score
    */
  override def dock(molA: Molecule, molB: Molecule, log: DockLog) = {
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
    if (params.bondForceWeight > 0)
      avgBondEnergy = scorer.getAvgBondEnergy(molA, molB)

    //println(s"Docking from pos $startingPos")
    startingPos+=1
    bestDock = molB.clone

    var i = 1
    while (!done){
      val forces = Profiler.time("getForces") { getForces(molA, molB) }             // it is a list of pairs (atomB, force)

      val (translation, _) = getTranslation(forces)
       molB.translate(translation)
      log.action(new Translate(translation))

      val (axis, angle) = getRotation(molB.getGeometricCentre, forces)
      molB.rotate(molB.getGeometricCentre, axis, angle)
      log.action(new Rotate(molB.getGeometricCentre, axis, angle))

      Profiler.time("updateState") {   updateState(i, molA, molB)       }  // updates done and all other algorithm control variables
      i+=1
    }
    (bestDock, bestScore)
  }

  private def updateState(i: Integer, molA: Molecule, molB: Molecule)  = {
    lastIterScore = currIterScore
    currIterScore = scorer.getScore(molA, molB, !approachPhase)

    if (currIterScore < lastIterScore) {
      maxAngle = maxAngle * 0.95
      maxTranslate = maxTranslate * 0.95
      softness = softness * 0.5
    }

    currWindowScore += currIterScore
    if (i % window == 0) {

      if (currWindowScore - lastWindowScore < params.threshold) {
        if ((decelerations >= params.maxDecelerations  && !approachPhase) || decelerations > 100) {
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
      } else if (currWindowScore > bestWindowScore + params.threshold) {
        decelerations = 0
        if (approachPhase) {
          maxAngle = initialDeltaAngle
          maxTranslate = initialDeltaSpace
        }
      }
      if (currWindowScore > bestWindowScore){
        bestScore = currIterScore
        bestDock = molB.clone
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
      moi += 0.1 * distToAxis//*distToAxis
    }

    val angle = if (moi > 0) Math.min(netTorqueNorm / moi, maxAngle) else 0.0
    //println(Math.toDegrees(angle))
    (axis, angle)
  }

  /** Returns a list of pairs (Atom, DenseVector) containing the net force for
    * each atom in molB.     */
  private def getForces(molA: Molecule, molB: Molecule) = {
   // val (minDist, maxDist) = molA.atoms(ignoreAHydrogens).foldLeft(Double.PositiveInfinity, Double.NegativeInfinity)
   //   { case ((min, max), a) => {val d = a.distTo(molB.getGeometricCentre); (Math.min(min, d), Math.max(max, d)) }}
    val maxCover = Int.MaxValue// minDist + (maxDist - minDist) * 0.2

    val forces = for (atomB <- molB.surfaceAtoms) yield {
      var forceOnAtomB = DenseVector.zeros[Double](3)

      for (atomA <- molA.atoms(params.ignoreAHydrogens)){
        val opt = optimalDistance(atomA, atomB, params.surface)
        val cover =
          if (approachPhase)
            maxCover
          else
            opt * minCoverage

        // cache these values for reuse:
        val (actual, bToA, dir) = atomB.distDifDir(atomA)
        if (actual <= cover) {
          if (atomA.isSurface) {
            forceOnAtomB :+= atomToAtomForce(atomA, atomB, molA, molB, opt, actual, dir)
            atomWithinMinCover |= actual <= opt * minCoverage;
          } else {
            val penalisation = getPenetrationPenaltyForce(atomA, atomB, actual, bToA, dir)
            if (penalisation.isDefined)
              forceOnAtomB :+= (1 - params.permeability) * penalisation.get
          }
        }
      }

      (atomB, forceOnAtomB)
    }

    forces.filter(t => norm(t._2) > 0)
  }

  private def inLineOfSight(atomA: Atom, atomB: Atom, molA: Molecule, molB: Molecule) = {
    // can't check against all atoms - too expensive, and actually we don't really want that
    val allAtoms = molA.atomsBoundTo(atomA) ++ molB.atomsBoundTo(atomB)
    val actualDistance = atomA.distTo(atomB)
    !allAtoms.exists(a => (a ne atomA) && (a ne atomB)
      && !a.isH
      && a.distTo(atomA) < actualDistance
      && a.distTo(atomB) < actualDistance
      && Geometry.distToLine(atomA.coords, atomB.coords, a.coords) < a.radius * 0.8
    )
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
//  private def molToAtomForce(molA: Molecule, atomB: Atom, molB: Molecule) = {
//    molA.surfaceAtoms(params.ignoreAHydrogens)
//      .map(atomA => atomToAtomForce(atomA, atomB, molA, molB))
//      .reduce((a,b) => a+b)
//  }

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
                              molA: Molecule, molB: Molecule,
                              optimal: Double, actual: Double,
                              dir: DenseVector[Double]): DenseVector[Double] = {

    // staggered in-place updates proved to be fastest option by far
    val force = DenseVector.zeros[Double](3)
    if (params.geometricForceWeight > 0)
      force :+= dir * (getGeometricForceNorm(atomA, atomB, actual, optimal) * params.geometricForceWeight)
    if (params.electricForceWeight > 0)
      force :+= dir * (getElectricForceNorm(atomA, atomB, actual, optimal) * params.electricForceWeight)
    if (params.hydrogenBondsForceWeight > 0)
      force :+= getHydrogenBondsForce(atomA, atomB, molA, molB, actual) * params.hydrogenBondsForceWeight
    if (params.bondForceWeight > 0)
      force :+= dir * (getBondForceNorm(atomA, atomB, actual, optimal) * params.bondForceWeight)
    force
  }


  /** Calculates the atomic force norm f such that the force that a exerts on b
    * is f * a normalized vector pointing from b to a.
    * Positive = attraction, negative = repulsion */
  private def getGeometricForceNorm(atomA: Atom, atomB: Atom, actualDistance: Double, optimal: Double) = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1
    if (atomA.isH || atomB.isH)
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
     else if (chargeProduct > 0 && actualDistance / optimal <= 2) // atoms farther than twice the optimal dont count
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

    if (params.ignoreAHydrogens)  // if A hydrogens are ignored, at least get reverse forces
      -1.0 * getForce(getVectToHBondSpot(atomB, atomA, molB))
    else
      getForce(getVectToHBondSpot(atomA, atomB, molA))
  }



  private def getBondForceNorm(atomA: Atom, atomB: Atom, actualDistance: Double, optimal: Double) = {
    val bondEnergy = BondEnergy(atomA.element, atomB.element)
    val normalized = actualDistance/optimal
    val frac = bondEnergy / avgBondEnergy
    if (frac > 1)
      explog(normalized) * frac
    else if (frac < 1 && actualDistance / optimal <= 2)
      minusExpOverX(normalized) * (1-frac)
    else
      0.0
  }

  private def getPenetrationPenaltyForce(atomA: Atom, atomB: Atom, actual: Double,
                                         bToA: DenseVector[Double], dir: DenseVector[Double]) = {
    if (!atomA.isH && !atomB.isH && (!atomA.isSurface || !atomB.isSurface)) {
      val penalizationDistance = (atomA.radius + atomB.radius) * 0.9
      if (actual < penalizationDistance) {
        val normalized = actual / penalizationDistance
        val penalisationForce = dir * Math.log(normalized)
        Some(penalisationForce)
      } else None
    } else None
  }
}
