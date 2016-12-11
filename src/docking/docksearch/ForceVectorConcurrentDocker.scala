package docking.docksearch

import breeze.linalg
import breeze.linalg._
import docking.dockscore.Scorer
import docking._
import docking.docksearch.forcevector.DockingParams
import opt._
import io.threadcso._
import model._
import profiling.Profiler

class ForceVectorConcurrentDockerDocker(params: DockingParams) extends Docker {

  val workers = 1

  val hBondDistance = 1.8

  val initialDeltaAngle = Math.toRadians(5)
  val initialDeltaSpace = 0.5
  val minDeltaSpace = 1.0           // minimum to be used only when the molecules are too far apart
  val minCoverage = 1.25            // how many times the optimal distance must always be considered, at least
  @volatile var maxCoverage = Double.PositiveInfinity
  @volatile var atomWithinMinCover = false

  var avgBondEnergy = 0.0
  var startingPos = 0
  @volatile var softness = 1.0

  var maxAngle = initialDeltaAngle
  var maxTranslate = initialDeltaSpace

  var currScore = Double.NegativeInfinity
  var lastScore = Double.NegativeInfinity
  var bestScore = Double.NegativeInfinity
  var currIterScore = 0.0
  var lastIterScore = Double.NegativeInfinity

  var done = false
  var decelerations = 0
  val window = 5                  // Assess score every window iterations
  @volatile var approachPhase = true
  @volatile var forceLongestDistance = Double.NegativeInfinity
  @volatile var forceShortestDistance = Double.PositiveInfinity

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]) = {
    println(getScore(molA, molB, true))
    // First, make sure both molecules are centered -> move B into the centre of A
    val centreVect = molA.getGeometricCentre - molB.getGeometricCentre
    log!new Translate(centreVect)
    molB.translate(centreVect)

    molB.rotate(molB.getGeometricCentre, DenseVector(0.0, 1.0, 0.0), Math.toRadians(180))
    log!new Rotate(molB.getGeometricCentre, DenseVector(0.0, 1.0, 0.0), Math.toRadians(180))

    //molB.rotate(molB.getGeometricCentre, DenseVector(1.0, 0.0, 0.0), Math.toRadians(180))
    //log!new Rotate(molB.getGeometricCentre, DenseVector(1.0, 0.0, 0.0), Math.toRadians(180))

    log!"save"

    val radius = molA.getRadius + molB.getRadius
    val initialConfigs = Geometry.sphereOrientations(radius, Math.toRadians(90))

    avgBondEnergy =
      (for (a <- molA.atoms; b <- molB.atoms) yield BondEnergy(a.element, b.element)).sum /
        (molA.atoms.size * molB.atoms.size)

    val result =
      initialConfigs.map(pos => {
        Profiler.time("dock") {  dockFromPos(molA, molB, pos, log) }
      }).maxBy(p => p._2)

    Profiler.report
    result
  }

  /**  Docks b into a from b's initial position and orientation, using force vectors.
    *   Try to minimize score
    */
  private def dockFromPos(molA: Molecule, b: Molecule,
                          pos: DenseVector[Double], log: ![Any]) = {
    maxAngle = initialDeltaAngle
    maxTranslate = initialDeltaSpace
    currScore = 0.0
    lastScore = Double.NegativeInfinity
    bestScore = Double.NegativeInfinity
    currIterScore = 0.0
    lastIterScore = Double.NegativeInfinity
    done = false
    decelerations = 0
    softness = 1.0
    approachPhase = true
    maxCoverage = Double.PositiveInfinity
    atomWithinMinCover = false
    forceLongestDistance = Double.NegativeInfinity
    forceShortestDistance = Double.PositiveInfinity

    println(s"Docking from pos $startingPos")
    startingPos+=1

    // move the molecule to the starting position
    val molB = b.clone
    log!"reset"
    molB.translate(pos)
    log!new Translate(pos)

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

    currScore += currIterScore
    if (i % window == 0) {

      if (currScore - lastScore < params.threshold) {
        if (decelerations >= params.maxDecelerations  && !approachPhase) {
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
            currScore = 0.0
            bestScore = Double.NegativeInfinity
            decelerations = 0
            currIterScore = 0
          }
        }
      } else if (currScore > bestScore + params.threshold) {
        decelerations = 0
        if (approachPhase) {
          maxAngle = initialDeltaAngle
          maxTranslate = initialDeltaSpace
        }
      }
      println(s"$i, score: $currScore, soft: $softness, decel $decelerations, approach: $approachPhase, maxT: $maxTranslate, maxC: $maxCoverage")
      lastScore = currScore
      bestScore = Math.max(bestScore, currScore)
      currScore = 0
    }

  }


  /** Returns a pair (vector, boolean). The vector has the translation,
    * the boolean indicates if minimum translation had to be applied
    */
  private def getTranslation(forces: Iterable[(Atom, DenseVector[Double])]) = {
    val forcesOnly = forces.map(pair => pair._2)
    val netForce = forcesOnly.reduce((a, b) => a + b)
    val totalForceAmount = forcesOnly.map(f => norm(f)).sum
    val minApplied = (totalForceAmount < minDeltaSpace) && approachPhase
    val translateDistance =
      if (minApplied)
        minDeltaSpace
      else
        Math.min(norm(netForce), maxTranslate)
    val translation = (netForce * translateDistance) / norm(netForce)
    (translation, minApplied)
  }

  private def getRotation(centre: DenseVector[Double], forces: Iterable[(Atom, DenseVector[Double])]) = {
    // torque in Euler axis/angle format is cross(r, force) - see http://web.mit.edu/8.01t/www/materials/modules/chapter21.pdf
    val torques = forces.map { case (atomB, force) =>
      val r = atomB.coords - centre // radius vector from centre to atom
      linalg.cross(r, force)
    }
    val netTorque = torques.reduce((a, b) => a + b)     // add all torques together
    val netTorqueNorm = norm(netTorque)
    val axis = netTorque / netTorqueNorm // normalize

    // calculate moment of inertia with all atoms assumed weight 0.1 otherwise angles are too small:
    var moi = 0.0
    for( (atom, force ) <- forces) {
      val distToAxis = Geometry.distToLine(centre, centre + netTorque, atom.coords)
      moi += 0.1 * distToAxis*distToAxis
    }

    val angle = Math.min(netTorqueNorm / moi, maxAngle)
    //println(Math.toDegrees(angle))
    (axis, angle)
  }

  /** Returns a list of pairs (Atom, DenseVector) containing the net force for
    * each atom in molB.     */
  private def getForces(molA: Molecule, molB: Molecule) = {
    forceLongestDistance = Double.NegativeInfinity
    @volatile var shortest = Double.PositiveInfinity

    val forces = for (atomB <- molB.surfaceAtoms) yield {
      /*var geometricForceOnB = DenseVector(0.0, 0.0, 0.0)
      var electricForceOnB = DenseVector(0.0, 0.0, 0.0)
      var hBondForceOnB = DenseVector(0.0, 0.0, 0.0)
      var bondStrengthForceOnB = DenseVector(0.0, 0.0, 0.0)
      var penetrationPenalty = DenseVector(0.0, 0.0, 0.0)*/
      @volatile var forceOnAtomB = DenseVector(0.0, 0.0, 0.0)

      val toW = N2N[(Atom, Atom)](1, workers, "tasks to workers")
      val fromW = N2N[DenseVector[Double]](workers, 1, "results from workers")
      val workerProcs = || (for (w <- 0 until workers) yield worker(molA, molB, toW, fromW))
      (workerProcs ||
        proc { shortest = Math.min(forceShortestDistance, pushTasks(atomB, molA, toW)) } ||
        proc { forceOnAtomB = collectResults(fromW) } )()        // run concurrently

      (atomB, forceOnAtomB)
    }

    forceShortestDistance = shortest
    maxCoverage = forceLongestDistance + maxTranslate
    forces.filter(t => norm(t._2) > 0)
  }

  def pushTasks(atomB: Atom, molA: Molecule, out: ![(Atom, Atom)]) = {
    var shortest = Double.PositiveInfinity
    for (atomA <- molA.atoms(params.ignoreAHydrogens)) {
      if (atomA.isSurface) {
        val cover = if (approachPhase)
          Math.max(forceShortestDistance * minCoverage, maxCoverage * softness)
        else
          optimalDistance(atomA, atomB) * minCoverage

        val actualDistance = atomA.distTo(atomB)
        if (actualDistance <= cover) {
          out!(atomA, atomB)
          //forceOnAtomB += atomToAtomForce(atomA, atomB, molA, molB)

          forceLongestDistance = Math.max(forceLongestDistance, actualDistance) // should be only if nonzero
          shortest = Math.min(shortest, actualDistance)
          atomWithinMinCover |= actualDistance <= optimalDistance(atomA, atomB) * minCoverage;
        }
      } else
        out!(atomA, atomB)
    }
    out.closeOut
    shortest
  }

  def worker(molA: Molecule, molB: Molecule,
             in: ?[(Atom, Atom)], out: ![DenseVector[Double]]) = proc {
    repeat {
      val (a, b) = in.?()
      if (a.isSurface)
        out ! atomToAtomForce(a, b, molA, molB)
      else
        out ! getPenetrationPenaltyForce(a, b)
    }
    in.closeIn
    out.closeOut
  }

  def collectResults(in: ?[DenseVector[Double]]) = {
    var i = 0
    var force = DenseVector(0.0, 0.0, 0.0);
    repeat {
      force = force + in.?()
      i+=1
    }
    in.closeIn
    //println(s"got $i with norm ${norm(force)}")
    force
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
                              molA: Molecule, molB: Molecule)= {
    val dif = atomA.coords - atomB.coords;                    // direction from b to a
    val actualDistance = norm(dif)
    val dir = dif / actualDistance;                                // normalized to length 1
    val optimal = optimalDistance(atomA, atomB)

    val geometricForce = dir * getGeometricForceNorm(atomA, atomB, actualDistance, optimal)
    val electricForce = dir * getElectricForceNorm(atomA, atomB, actualDistance, optimal)
    val hydrogenBondsForce = getHydrogenBondsForce(atomA, atomB, molA, molB, actualDistance)
    val bondForce = dir * getBondForceNorm(atomA, atomB, actualDistance, optimal)

    val force =     // weighted result:
      geometricForce * params.geometricForceWeight +
        electricForce * params.electricForceWeight +
        hydrogenBondsForce * params.hydrogenBondsForceWeight +
        bondForce * params.bondForceWeight

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
      val force = dir * Math.log(normalized)
      (1 - params.permeability) * force
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
    for (a <- molA.surfaceAtoms(params.ignoreAHydrogens) ){
      aCount += 1
      for(b <- molB.surfaceAtoms) {
        bCount += 1
        val actualDistance = a.distTo(b)
        if (!onlyTargetRadius || actualDistance <= optimalDistance(a, b) * minCoverage) {
          val n = actualDistance / optimalDistance(a, b)
          if (!a.isElement("H") && !b.isElement("H"))
            totalGeometricScore += explog2(n * explog2.maxX) / explog2.maxY

          totalElectricScore += (explog2(n * explog2.maxX) / explog2.maxY) * a.partialCharge * b.partialCharge * -1.0

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
      totalGeometricScore * params.geometricForceWeight +
        totalElectricScore * params.electricForceWeight +
        totalHBondScore * params.hydrogenBondsForceWeight +
        totalBondStrengthScore * params.bondForceWeight
    //println(s"SCORES: geo: ${totalGeometricScore * geometricForceWeight}, electric: ${totalElectricScore * electricForceWeight}, hbond: ${totalHBondScore * hydrogenBondsForceWeight}, bondstrength: ${totalBondStrengthScore * bondForceWeight}")
    totalScore / (Math.min(aCount, bCount))
  }


  /* --- Distance Functions --- */

  private def optimalDistance(a: Atom, b: Atom) = {
    if (canHBond(a, b) || canHBond(b, a))
      hBondDistance + 2 * params.surface
    else
      (a.radius + b.radius) * 0.9 + 2 * params.surface
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
