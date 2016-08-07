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


  val initialDeltaAngle = Math.toRadians(20) // 20 degrees in radians
  val initialDeltaSpace = 1.0
  val minDeltaSpace = 1.0       // minimum to be used only when the molecules are too far apart

  var avgBondEnergy = 0.0
  var iter =0

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
    }).minBy(p => p._2)
  }

  /**  Docks b into a from b's initial position and orientation, using force vectors.
    *   Try to minimize score
    */
  private def dockFromPos(molA: Molecule, b: Molecule,
                           pos: DenseVector[Double],
                           threshold: Double,
                           avgBondEnergy: Double, log: ![Any]) = {

    println(s"Docking from pos $iter")
    iter+=1

    // move the molecule to the starting position
    val molB = b.clone
    log!"reset"
    molB.translate(pos)
    log!new Translate(pos)

    var maxAngle = initialDeltaAngle
    var maxTranslate = initialDeltaSpace

    var currScore = Double.PositiveInfinity
    var lastScore = Double.PositiveInfinity

    var scoreChange = 0.0                             // how much has the score changed in the last iteration
    var lastScoreChange = Double.NegativeInfinity     // how much had the score changed in the previous interation

    var minTranslationApplied = false


    /*while (lastScore == Double.PositiveInfinity               // first and second iterations
      || scoreChange < -threshold                             // good rate of score improvement
      || scoreChange > 0                                      // decelerated
      || minTranslationApplied  ) {                           // min translation had to be applied because the molecules were too far apart
    */
    while (true) {
      val forcesAndScores = getForces(molA, molB, avgBondEnergy)       // it is a list of tuples (atomB, force, score)
      val forces = forcesAndScores.map(e => e._2)
      val totalForceAmount = forces.map(f => norm(f)).sum

      val (translation, m) = getTranslation(forces, totalForceAmount, maxTranslate)
      minTranslationApplied = m
      val (axis, angle) = getRotation(molB, maxAngle, forcesAndScores)

      molB.translate(translation)
      log!new Translate(translation)

      molB.rotate(molB.getGeometricCentre, axis, angle)
      log!new Rotate(molB.getGeometricCentre, axis, angle)

      lastScore = currScore
      currScore = forcesAndScores.map(e => e._3).sum
      lastScoreChange = scoreChange
      scoreChange = currScore - lastScore

      if (scoreChange > 0){     // decelerate
        maxAngle = maxAngle * 0.9
        maxTranslate = maxTranslate * 0.9
      }
      //println(s"score: $currScore, minTrans: $minTranslationApplied, maxAngle $maxAngle, maxTrans: $maxTranslate")
    }
    (molB, currScore)
  }

  /** Returns a pair (vector, boolean). The vector has the translation,
    * the boolean indicates if minimum translation had to be applied
    */
  private def getTranslation(forces: Seq[DenseVector[Double]],
                             totalForceAmount: Double,
                             maxTranslate: Double) = {
    val netForce = forces.reduce((a, b) => a + b)
    val minApplied = totalForceAmount < minDeltaSpace
    val translateDistance =
      if (minApplied)
        minDeltaSpace
      else
        Math.min(norm(netForce), maxTranslate)
    val translation = (netForce * translateDistance) / norm(netForce)
    (translation, minApplied)
  }

  private def getRotation(molB: Molecule, maxAngle: Double,
                          forces: Seq[(Atom, DenseVector[Double], Double)]) = {
    // torque in Euler axis/angle format is cross(r, force) - see http://web.mit.edu/8.01t/www/materials/modules/chapter21.pdf
    val torques = forces.map { case (atomB, force, score) =>
      val r = atomB.coords - molB.getGeometricCentre // radius vector from centre to atom
      linalg.cross(r, force)
    }
    val netTorque = torques.reduce((a, b) => a + b)     // add all torques together
    val netTorqueNorm = norm(netTorque)
    val axis = netTorque / netTorqueNorm // normalize
    val angle = Math.min(netTorqueNorm, maxAngle)
    (axis, angle)
  }

  /** Returns a list of pairs (Atom, DenseVector) containing the net force for
    * each atom in molB.     */
  private def getForces(molA: Molecule, molB: Molecule, avgBondEnergy: Double) = {
    molB.Atoms
      .filter(atomB => !(ignoreHydrogen && atomB.isElement("H")))
      .map(atomB => {
        val forceAndScore = molToAtomForce(molA, atomB)
        (atomB, forceAndScore._1, forceAndScore._2)
      })
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
  private def molToAtomForce(molA: Molecule, atomB: Atom) = {
    molA.Atoms
      .filter(atomA => !(ignoreHydrogen && atomA.isElement("H")))
      .map(atomA => atomToAtomForce(atomA, atomB, avgBondEnergy))
      .reduce((a, b) => (a._1 + b._1, a._2 + b._2))
  }

  /** calculates the force that atomA excerts on atomB
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
    val (atomicForceNorm, atomicScore) = getAtomicForceNorm(atomA, atomB)
    val atomicForce = dir * atomicForceNorm

    // electric force:
    val electricForce = dir * getElectricForceNorm(atomA, atomB)

    // bond force:
    val bondForce = dir * getBondForceNorm(atomA, atomB, avgBondEnergy)

    // weighted result:
    val force = atomicForce * atomicForceWeight +
      electricForce * electricForceWeight +
      bondForce * bondForceWeight

    val score = atomicScore
    (force, score)
  }


  /** Calculates the atomic force norm f such that the force that a exerts on b
    * is f * a normalized vector pointing from b to a.
    * Positive = attraction, negative = repulsion */
  private def getAtomicForceNorm(atomA: Atom, atomB: Atom) = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1

    val optimalDistance = atomA.radius + atomB.radius + 2 * surface // optimal distance
    val actualDistance = atomA.distTo(atomB)
    val normalized = actualDistance/optimalDistance
    val force = explog(normalized)//val force = explog(normalized*1.327864011995167)

    // Max of explog is reached at 1.327864011995167
    //val score = -explog(normalized * 1.327864011995167)

    // Max of expsquare is reached at sqrt(2)
    val score = -expsquare(normalized * Math.sqrt(2.0))

    (force, score)
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



  /* --- Distance Functions --- */

  private def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)

  private def explog(x: Double) = Math.exp(-Math.pow(x,2))*Math.log(x)

  private def minusExpOverX(x: Double) = - Math.exp(-x) * 0.1 / x
}
