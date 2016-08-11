package docking.docksearch

import breeze.linalg
import breeze.linalg._

import docking.dockscore.Scorer
import docking._
import opt._
import io.threadcso._
import model._


// FIX THE DESIGN A BIT!!!
class ForceVectorDocker(val surface: Double, val maxDecelerations: Int = 10,
                        val ignoreHydrogen: Boolean = false,
                        val atomicForceWeight: Double = .34,
                        val electricForceWeight: Double = .33,
                        val bondForceWeight: Double = .33) extends Docker {


  val initialDeltaAngle = Math.toRadians(5) //
  val initialDeltaSpace = 0.5
  val minDeltaSpace = 1.0       // minimum to be used only when the molecules are too far apart

  var avgBondEnergy = 0.0
  var startingPos =0

  var softness = 1.0
  var coverage = 1.0          // 0 = min, 1 = max

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]) = {
    // First, make sure both molecules are centered -> move B into the centre of A
    val centreVect = molA.getGeometricCentre - molB.getGeometricCentre
    log!new Translate(centreVect)
    molB.translate(centreVect)
    log!"save"

    val radius = molA.getRadius + molB.getRadius
    val initialConfigs = Geometry.sphereOrientations(radius, Math.toRadians(90))

    val avgBondEnergy =
      (for (a <- molA.Atoms; b <- molB.Atoms) yield BondEnergy(a.element, b.element)).sum /
        (molA.Atoms.size + molB.Atoms.size)

    initialConfigs.map(pos => {
      dockFromPos(molA, molB, pos, 1.0e-5, avgBondEnergy, log)
    }).maxBy(p => p._2)
  }

  /**  Docks b into a from b's initial position and orientation, using force vectors.
    *   Try to minimize score
    */
  private def dockFromPos(molA: Molecule, b: Molecule,
                           pos: DenseVector[Double],
                           threshold: Double,
                           avgBondEnergy: Double, log: ![Any]) = {

    println(s"Docking from pos $startingPos")
    startingPos+=1

    // move the molecule to the starting position
    val molB = b.clone
    log!"reset"
    molB.translate(pos)
    log!new Translate(pos)

    var maxAngle = initialDeltaAngle
    var maxTranslate = initialDeltaSpace

    var currScore = Double.NegativeInfinity
    var lastScore = Double.NegativeInfinity
    var bestScore = Double.NegativeInfinity
    var currIterScore = 0.0
    var lastIterScore = Double.NegativeInfinity

    var scoreChange = 0.0                             // how much has the score changed in the last iteration

    var done = false
    var decelerations = 0
    val span = 5

    softness = 1.0
    coverage = 1.0
    //coverage = norm(molA.getGeometricCentre - molB.getGeometricCentre)

    /*while (lastScore == Double.PositiveInfinity               // first and second iterations
      || scoreChange > threshold                             // good rate of score improvement
      || scoreChange > 0                                      // decelerated
      || minTranslationApplied  ) {                           // min translation had to be applied because the molecules were too far apart
*/
    var i = 1
    var approachPhase = true

    while (!done){
      val forces = getForces(molA, molB, avgBondEnergy)             // it is a list of pairs (atomB, force)

      var (translation, m) = getTranslation(forces, maxTranslate)
      //translation = reviewTranslation(translation, molA, molB)

      val (axis, angle) = getRotation(molB.getGeometricCentre, maxAngle, forces)

      molB.translate(translation)
      log!new Translate(translation)

      molB.rotate(molB.getGeometricCentre, axis, angle)
      log!new Rotate(molB.getGeometricCentre, axis, angle)

      lastIterScore = currIterScore
      currIterScore = getScore(molA, molB, !approachPhase)

      if (currIterScore - lastIterScore < 0) {
        maxAngle = maxAngle * 0.95
        maxTranslate = maxTranslate * 0.95
        softness = softness * 0.8
      }

      currScore += currIterScore
      if (i % span == 0) {
        scoreChange = currScore - lastScore

        if (scoreChange < threshold) {
          if (decelerations >= 3  && !approachPhase) {
            done = true
          } else {
            // decelerate
            maxAngle = maxAngle * 0.75
            maxTranslate = maxTranslate * 0.75

            softness = Math.max(softness - 0.2, 0)
            //coverage = Math.max(coverage - 0.1, 0)

            decelerations += 1

            if (approachPhase && softness <= 0){
              // end approach phase
              coverage = 0.0
              approachPhase = false
              currScore = Double.NegativeInfinity
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

        println(s"$i, score: $currScore, soft: $softness, cover: $coverage, decel $decelerations, approach: $approachPhase")

        lastScore = currScore
        bestScore = Math.max(bestScore, currScore)
        if (!done) currScore = 0
      }

      i+=1
    }
    (molB, getScore(molA, molB, true))
  }

  private def reviewTranslation(t: DenseVector[Double], molA: Molecule, molB: Molecule) = {
    for (a <- molA.Atoms.filter(atomA => !(ignoreHydrogen && atomA.isElement("H")))) {
      for (b <- molB.Atoms.filter(atomB => !(ignoreHydrogen && atomB.isElement("H")))) {
        if (a.distTo(b) < optimalDistance(a, b)) {
          val dif = a.coords - b.coords // direction from b to a
          val dir = dif / norm(dif)

          val component = t dot dir
          if (component > 0)
            t += dir * (-component)
        }
      }
    }
      t
  }

  /** Returns a pair (vector, boolean). The vector has the translation,
    * the boolean indicates if minimum translation had to be applied
    */
  private def getTranslation(forces: Seq[(Atom, DenseVector[Double])],
                             maxTranslate: Double) = {
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

  private def getRotation(centre: DenseVector[Double], maxAngle: Double,
                          forces: Seq[(Atom, DenseVector[Double])]) = {
    // torque in Euler axis/angle format is cross(r, force) - see http://web.mit.edu/8.01t/www/materials/modules/chapter21.pdf
    val torques = forces.map { case (atomB, force) =>
      val r = atomB.coords - centre // radius vector from centre to atom
      linalg.cross(r, force)
    }
    val netTorque = torques.reduce((a, b) => a + b)     // add all torques together
    val netTorqueNorm = norm(netTorque)
    val axis = netTorque / netTorqueNorm // normalize

    // calculate moment of inertia with all atoms assumed weight 1:
    var moi = 0.0
    for( (atom, force ) <- forces){
      val distToAxis = norm(linalg.cross(atom.coords - centre, atom.coords - (centre+netTorque))) / norm(netTorque-centre)  //http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
      moi += 1.0* distToAxis*distToAxis
    }

    val angle = Math.min(netTorqueNorm / moi, maxAngle)
    (axis, angle)
  }

  /** Returns a list of pairs (Atom, DenseVector) containing the net force for
    * each atom in molB.     */
  private def getForces(molA: Molecule, molB: Molecule, avgBondEnergy: Double) = {
    /*molB.Atoms
      .filter(atomB => !(ignoreHydrogen && atomB.isElement("H")))
      .map(atomB => (atomB, molToAtomForce(molA, atomB, coverage)))
*/
    val all =
      for (a <- molA.Atoms.filter(atomA => !(ignoreHydrogen && atomA.isElement("H"))) ;
        b <- molB.Atoms.filter(atomB => !(ignoreHydrogen && atomB.isElement("H"))))
          yield (a, b, atomToAtomForce(a, b, avgBondEnergy))

    //val sorted = all.sortBy(t => t match { case (a, b, f) => a.distTo(b) })
    //val toTake = Math.round(sorted.size * coverage).toInt
    //val forces = sorted.take(toTake)
    val forces = all.filter(t => t match { case (a, b, f) =>
      val maxCover = b.distTo(molA.getGeometricCentre)          // can be cached!
      val minCover = optimalDistance(a, b) * 1.25;
      val cover = (maxCover-minCover)*coverage + minCover
      a.distTo(b) <= cover
    })

    var dict = Map[Atom, DenseVector[Double]]()
    for (t <- forces)
      t match {
        case (a, b, f) => dict += (b -> (
          dict.get(b) match {
            case None => f
            case Some(f2) => f + f2
          }))
      }
    dict.toList
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
  private def molToAtomForce(molA: Molecule, atomB: Atom) = {
    molA.Atoms
      .filter(atomA => !(ignoreHydrogen && atomA.isElement("H")))
      .map(atomA => atomToAtomForce(atomA, atomB, avgBondEnergy))   //TODO: use avgbond from argument!
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
  private def atomToAtomForce(atomA: Atom, atomB: Atom, avgBondEnergy: Double) = {
    val dif = atomA.coords - atomB.coords;                    // direction from b to a
    val dir = dif / norm(dif);                                // normalized to length 1

    // Atomic force:
    val atomicForceNorm = getAtomicForceNorm(atomA, atomB)
    val atomicForce = dir * atomicForceNorm

    // electric force:
    val electricForce = dir * getElectricForceNorm(atomA, atomB)

    // bond force:
    val bondForce = dir * getBondForceNorm(atomA, atomB, avgBondEnergy)

    // weighted result:
    val force = atomicForce * atomicForceWeight +
      electricForce * electricForceWeight +
      bondForce * bondForceWeight

    force
  }


  /** Calculates the atomic force norm f such that the force that a exerts on b
    * is f * a normalized vector pointing from b to a.
    * Positive = attraction, negative = repulsion */
  private def getAtomicForceNorm(atomA: Atom, atomB: Atom) = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1

    val actualDistance = atomA.distTo(atomB)
    val normalized = actualDistance/optimalDistance(atomA, atomB)

    val force = explog(normalized / (softness * 3 + 1))
    force
  }

  private def getElectricForceNorm(atomA: Atom, atomB: Atom): Double = {
    val chargeProduct = atomA.partialCharge*atomB.partialCharge
    val dist = atomA.distTo(atomB)

    //Math.exp(-dist) * (- Math.signum(chargeProduct))  // For different charges, return positive value, else negative
    // TODO: shouldn't we take into account the amount chargeProduct??
    if (chargeProduct < 0)
      explog(dist)
    else if (chargeProduct > 0)
      minusExpOverX(dist)
    else
      0.0

    /*if (chargeProduct < 0)                                  // different charges, attract
      Math.max(atomicForceNorm, 0.0)                        // the atomic force, if distance > optimal
    else if (chargeProduct > 0)                             // equal charges, reject
      Math.min(atomicForceNorm, 0.0)                        // the atomic force, if distance < optimal
    else                                                    // no charge
      0.0*/
  }

  private def getBondForceNorm(atomA: Atom, atomB: Atom, avgBondEnergy: Double) = {
    val bondEnergy = BondEnergy(atomA.element, atomB.element)
    val dist = atomA.distTo(atomB)
    val diff = bondEnergy - avgBondEnergy
    if (diff > 0)
      explog(dist)
    else if (diff < 0)
      minusExpOverX(dist)
    else
      0.0
  }

  /* --- Scoring --- */
  //TODO: scoring should be independent of the number of atoms!
  private def getScore(molA: Molecule, molB: Molecule, onlyTargetRadius: Boolean) = {
    var totalScore = 0.0
    for (a <- molA.Atoms.filter(atomA => !(ignoreHydrogen && atomA.isElement("H")))){
      for(b <- molB.Atoms.filter(atomB => !(ignoreHydrogen && atomB.isElement("H")))) {
        val actualDistance = a.distTo(b)
        totalScore += (
          if (!onlyTargetRadius || actualDistance <= optimalDistance(a, b) * 1.25) {
            //TODO: put 1.1 in a constant or something...
            //expsquare(actualDistance * Math.sqrt(2.0) / optimalDistance(a, b))
            val x = (actualDistance / optimalDistance(a, b)) * 1.6290055996317214
            Math.exp(-Math.pow(x - 1, 2)) * Math.log(x) //  exp(-(x-1)^2)*log(x)  peaks at 1.6290055996317214
          } else
            0.0
          )
      }
    }
    totalScore
  }




  /* --- Distance Functions --- */

  private def optimalDistance(a: Atom, b: Atom) = a.radius + b.radius + 2 * surface

  private def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)

  private def explog(x: Double) = Math.exp(-Math.pow(x,2))*Math.log(x)    // Max of explog is reached at 1.327864011995167

  private def minusExpOverX(x: Double) = - Math.exp(-x) * 0.1 / x
}
