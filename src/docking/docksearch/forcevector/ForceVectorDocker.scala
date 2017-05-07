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

    startingPos+=1
    bestDock = molB.clone

    // initialise forces:
    var forces = List[Force]()
    if (params.geometricForceWeight > 0)
      forces ::= new GeometricForce(params.geometricForceWeight, params.surface)
    if (params.electricForceWeight > 0)
      forces ::= new ElectricForce(params.electricForceWeight, params.surface)
    if (params.hydrogenBondsForceWeight > 0)
      forces ::= new HydrogenBondForce(params.hydrogenBondsForceWeight, params.ignoreAHydrogens)

    var i = 1
    while (!done){
      val forceVectors = Profiler.time("getForces") { getForcesOnBAtoms(forces, molA, molB) }             // it is a list of pairs (atomB, force)

      val (translation, _) = getTranslation(forceVectors)
      molB.translate(translation)
      log.action(new Translate(translation))

      val (axis, angle) = getRotation(molB.getGeometricCentre, forceVectors)
      molB.rotate(molB.getGeometricCentre, axis, angle)
      log.action(new Rotate(molB.getGeometricCentre, axis, angle))

      updateState(i, molA, molB, forces) // updates done and all other algorithm control variables
      i+=1
    }
    (bestDock, bestScore)
  }

  private def updateState(i: Integer, molA: Molecule, molB: Molecule, forces: Seq[Force])  = {
    lastIterScore = currIterScore
    currIterScore = forces.map(f => f.weightedScore).sum
    for (f <- forces) f.endIter()  // reset forces score counters

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
  private def getForcesOnBAtoms(forces: Seq[Force], molA: Molecule, molB: Molecule) = {
   // val (minDist, maxDist) = molA.atoms(ignoreAHydrogens).foldLeft(Double.PositiveInfinity, Double.NegativeInfinity)
   //   { case ((min, max), a) => {val d = a.distTo(molB.getGeometricCentre); (Math.min(min, d), Math.max(max, d)) }}
    val maxCover = Int.MaxValue// minDist + (maxDist - minDist) * 0.2

    val forceVectors = for (atomB <- molB.surfaceAtoms) yield {
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
            forceOnAtomB :+= atomToAtomForce(forces, atomA, atomB, molA, molB, opt, actual, dir)
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

    forceVectors.filter(t => norm(t._2) > 0)
  }

  /** calculates the force that atomA exerts on atomB */
  private def atomToAtomForce(forces: Seq[Force], atomA: Atom, atomB: Atom,
                              molA: Molecule, molB: Molecule,
                              optimal: Double, actual: Double,
                              dir: DenseVector[Double]): DenseVector[Double] = {
    val totalForceVector = DenseVector.zeros[Double](3)
    for (force <- forces) {
      val forceVector = force(molA, molB, atomA, atomB, actual, dir)
      if (forceVector.isDefined)
        totalForceVector :+= forceVector.get
    }
    totalForceVector
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
