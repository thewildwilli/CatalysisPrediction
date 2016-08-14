package docking.docksearch

import breeze.linalg
import breeze.linalg._
import docking.dockscore.Scorer
import docking._
import opt._
import io.threadcso._
import model._
import profiling.Profiler

class ForceVectorDocker(val surface: Double, val maxDecelerations: Int = 10,
                        val ignoreHydrogen: Boolean = false,
                        val threshold: Double = 1.0e-5,
                        val atomicForceWeight: Double = .34,
                        val electricForceWeight: Double = .33,
                        val bondForceWeight: Double = .33) extends Docker {


  val hBondDistance = 1.8

  val initialDeltaAngle = Math.toRadians(5)
  val initialDeltaSpace = 0.5
  val minDeltaSpace = 1.0           // minimum to be used only when the molecules are too far apart
  val minCoverage = 1.25            // how many times the optimal distance must always be considered, at least

  var avgBondEnergy = 0.0
  var startingPos = 0
  var softness = 1.0

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
  var approachPhase = true

  var hBondSpots: Map[Int, DenseVector[Double]] = null

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]) = {
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
      (for (a <- molA.Atoms; b <- molB.Atoms) yield BondEnergy(a.element, b.element)).sum /
        (molA.Atoms.size * molB.Atoms.size)

    hBondSpots = getHBondMap(molA)

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

    if (currIterScore - lastIterScore < 0) {
      maxAngle = maxAngle * 0.95
      maxTranslate = maxTranslate * 0.95
      softness = softness * 0.5
    }

    currScore += currIterScore
    if (i % window == 0) {

      if (currScore - lastScore < threshold) {
        if (decelerations >= 10  && !approachPhase) {
          done = true
        } else {
          // decelerate
          maxAngle = maxAngle * 0.9
          maxTranslate = maxTranslate * 0.9
          softness = Math.max(softness - 0.2, 0)
          decelerations += 1

          if (approachPhase && softness <= 0){
            // end approach phase
            approachPhase = false
            currScore = 0.0
            bestScore = Double.NegativeInfinity
            decelerations = 0
            currIterScore = 0
          }
        }
      } else if (currScore > bestScore + threshold) {
        decelerations = 0
        if (approachPhase) {
          maxAngle = initialDeltaAngle
          maxTranslate = initialDeltaSpace
        }
      }
      println(s"$i, score: $currScore, soft: $softness, decel $decelerations, approach: $approachPhase")
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
    val minApplied = totalForceAmount < minDeltaSpace
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
      //val distToAxis = norm(linalg.cross(atom.coords - centre, atom.coords - (centre+netTorque))) / norm(netTorque-centre)  //http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
      val distToAxis = Geometry.distToLine(centre, centre + netTorque, atom.coords)
      moi += 0.1 * distToAxis*distToAxis
    }

    val angle = Math.min(netTorqueNorm / moi, maxAngle)
    (axis, angle)
  }

  /** Returns a list of pairs (Atom, DenseVector) containing the net force for
    * each atom in molB.     */
  private def getForces(molA: Molecule, molB: Molecule) = {
    for (atomB <- molB.Atoms if !(ignoreHydrogen && atomB.isElement("H"))) yield {
      val distToACenter = atomB.distTo(molA.getGeometricCentre)
      var forceOnAtomB = DenseVector(0.0, 0.0, 0.0)
      for (atomA <- molA.Atoms if !(ignoreHydrogen && atomB.isElement("H"))){
        val cover = if (approachPhase) distToACenter else optimalDistance(atomA, atomB) * minCoverage;
        if (atomA.distTo(atomB) <= cover)
          forceOnAtomB += atomToAtomForce(atomA, atomB, molA, molB)
      }
      (atomB, forceOnAtomB)
    }

   /* molB.Atoms
      .filter(atomB => !(ignoreHydrogen && atomB.isElement("H")))
      .map(atomB => (atomB, molToAtomForce(molA, atomB)))*/
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
  private def molToAtomForce(molA: Molecule, atomB: Atom, molB: Molecule) = {
    molA.Atoms
      .filter(atomA => !(ignoreHydrogen && atomA.isElement("H")))
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

    // if someone in the way, then no force
    for (other <- molA.Atoms ++ molB.Atoms if (other ne atomA) && (other ne atomB)){
      val t = Geometry.project(atomB.coords, atomA.coords, other.coords)
      if (t > 0 && t < 1) {
        val d = Geometry.distToLine(atomB.coords, atomA.coords, other.coords)
        if (d < other.radius) {
          //println(s"Force from ${atomB.id} to ${atomA.id} is none because ${other.id} is in the way")
          return DenseVector(0.0, 0.0, 0.0)
        }
      }
    }


    // Atomic force:
    val atomicForceNorm = getAtomicForceNorm(actualDistance, optimalDistance(atomA, atomB))
    val atomicForce = dir * atomicForceNorm

    // electric force:
    val electricForce = getElectricForce(atomA, atomB, molA, molB, actualDistance)

    // bond force:
    val bondForce = dir * getBondForceNorm(atomA, atomB, actualDistance)

    // weighted result:
    val force =
      atomicForce * atomicForceWeight +
      electricForce * electricForceWeight +
      bondForce * bondForceWeight

    force
  }


  /** Calculates the atomic force norm f such that the force that a exerts on b
    * is f * a normalized vector pointing from b to a.
    * Positive = attraction, negative = repulsion */
  private def getAtomicForceNorm(actualDistance: Double, optimal: Double) = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1

    val normalized = actualDistance/optimal
    val force = explog(normalized / (softness * 3 + 1))
    force
  }

  private def getElectricForce(atomA: Atom, atomB: Atom, molA: Molecule,
                                   molB: Molecule, actualDistance: Double) = {

    if (canHydrogenBond(atomA, atomB, molA, molB)){
      val hBondSpot = getHydrogenBondSpot(atomA, molA)
      val bToSpot = hBondSpot - atomB.coords
      val dist = norm(bToSpot)
      val bToSpotDir = bToSpot / dist
      //val force = - (atomA.partialCharge * atomB.partialCharge) / Math.pow(dist, 1)
      val force = - dist * Math.exp(-Math.pow(dist, 1)) * atomA.partialCharge * atomB.partialCharge
      bToSpotDir * force
    } else
      DenseVector(0.0, 0.0, 0.0)

    //val chargeProduct = atomA.partialCharge*atomB.partialCharge
    //-chargeProduct / (actualDistance*actualDistance) // Coulomb

   /* (if (chargeProduct < 0)
      explog(actualDistance)
    else if (chargeProduct > 0)
      minusExpOverX(actualDistance)
    else
     0.0) * Math.abs(chargeProduct)
*/
  }

  private def canHydrogenBond(atomA: Atom, atomB: Atom, molA: Molecule, molB: Molecule) = {
    def canHBond(a: Atom, b: Atom, molA: Molecule, molB: Molecule) =
      (a.isElement("H") && a.partialCharge > 0 && b.partialCharge < 0
        && (b.isElement("F") || b.isElement("O") || b.isElement("N"))
        && a.bonds.size == 1
        && {
          val aNeighbour = molA(a.bonds(0));
          aNeighbour.partialCharge < 0 &&
            (aNeighbour.isElement("F") || aNeighbour.isElement("O") || aNeighbour.isElement("N"))
          }
      )
    canHBond(atomA, atomB, molA, molB) //|| canHBond(atomB, atomA, molB, molA)
  }

  // where atomA is a H and atom B is F, O, or N
  private def getHydrogenBondSpot(atomA: Atom, molA: Molecule) = {
    val aNeighbour = molA(atomA.bonds(0))
    val dirToNeighour = aNeighbour.coords - atomA.coords
    val toSpot = - (dirToNeighour / norm(dirToNeighour)) * 1.8 // 1.8 Angtrom in the opposite direction
    val spot = atomA.coords + toSpot
    spot
  }

  private def getBondForceNorm(atomA: Atom, atomB: Atom, actualDistance: Double) = {
    val bondEnergy = BondEnergy(atomA.element, atomB.element)
    val diff = bondEnergy - avgBondEnergy
    //TODO: take into account how far from average
    if (diff > 0)
      explog(actualDistance)
    else if (diff < 0)
      minusExpOverX(actualDistance)
    else
      0.0
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
    var totalScore = 0.0
    var aCount = 0
    var bCount = 0
    for (a <- molA.Atoms.filter(atomA => !(ignoreHydrogen && atomA.isElement("H")))){
      aCount += 1
      for(b <- molB.Atoms.filter(atomB => !(ignoreHydrogen && atomB.isElement("H")))) {
        bCount += 1
        val actualDistance = a.distTo(b)
        totalScore += (
          if (!onlyTargetRadius || actualDistance <= optimalDistance(a, b) * minCoverage) {
            val x = actualDistance / optimalDistance(a, b)
            val geometricScore = explog2(x * explog2.maxX) / explog2.maxY   // ignore H for geometric?

            // optimal distance for hydrogen bond is 1.97 Angstrom

            val electricScore = if (canHydrogenBond(a, b, molA, molB)) {
              val hBondSpot = getHydrogenBondSpot(a, molA)
              val bToSpot = hBondSpot - b.coords
              val bToSpotNorm = norm(bToSpot)
              -(a.partialCharge * b.partialCharge) / (bToSpotNorm*bToSpotNorm)


              //-(a.partialCharge * b.partialCharge) / (actualDistance * actualDistance)
              //val normalized = actualDistance / 1.8
              //- explog2(normalized * explog2.maxX) / explog2.maxY * a.partialCharge * b.partialCharge
            } else
              0.0

            geometricScore * atomicForceWeight + electricScore * electricForceWeight
          } else
            0.0
          )
      }
    }
    totalScore / (Math.min(aCount, bCount))       // over atom count of smallest molecule
  }





  /* --- Distance Functions --- */

  private def optimalDistance(a: Atom, b: Atom) = a.radius + b.radius + 2 * surface

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


  // ----
  def getHBondMap(mol: Molecule) = {
    var result = Map[Int, DenseVector[Double]]()
    for (atom <- mol.Atoms){
      if ( atom.isElement("H") && atom.partialCharge > 0
      && atom.bonds.size == 1
      && {
        val aNeighbour = mol(atom.bonds(0));
        aNeighbour.partialCharge < 0 &&
          (aNeighbour.isElement("F") || aNeighbour.isElement("O") || aNeighbour.isElement("N"))
      }) {
        val hBondSpot = getHydrogenBondSpot(atom, mol)
        val closestAtom = mol.Atoms.minBy(a => a.distTo(hBondSpot))
        if (closestAtom.distTo(hBondSpot) >= closestAtom.radius) {
          result += (atom.id -> hBondSpot)
        }
      }
    }
    result
  }
}
