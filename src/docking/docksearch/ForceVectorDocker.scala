package docking.docksearch

import breeze.linalg
import breeze.linalg._

import docking.dockscore.Scorer
import docking._
import opt._
import io.threadcso._
import model._

// Created by Ernesto on 08/06/2016.
// FIX THE DESIGN A BIT!!!
class ForceVectorDocker(val surface: Double, val maxDecays: Int = 10) extends Docker {
  val initialDeltaAngle = Math.toRadians(20) // 20 degrees in radians
  val initialDeltaSpace = 1.0
  val minDeltaSpace = 1.0       // minimum to be used only when the molecules are too far apart

  var iter =0

  override def dock(molA: Molecule, molB: Molecule, scorer: Scorer, log: ![Any]) = {
    // First, make sure both molecules are centered -> move B into the centre of A
    val centreVect = molA.getGeometricCentre - molB.getGeometricCentre
    log!new Translate(centreVect)
    molB.translate(centreVect)
    log!"save"

    val radius = molA.getRadius + molB.getRadius
    val initialConfigs = Geometry.sphereOrientations(radius, Math.toRadians(90))

    initialConfigs.map(pos => {
      dockFromPos(molA, molB, pos, scorer, 1.0e-5, log)
    }).maxBy(scorer.score)
  }

  /**  Docks b into a from b's initial position and orientation, using force vectors
    */
  private def dockFromPos(molA: Molecule, b: Molecule,
                           pos: DenseVector[Double], scorer: Scorer,
                           threshold: Double, log: ![Any]) = {

    println(s"Docking from pos $iter")
    iter+=1

    // move the molecule to the starting position
    val molB = b.clone
    log!"reset"
    molB.translate(pos)
    log!new Translate(pos)

    val state = new DockingState(molA, molB)

    var maxAngle = initialDeltaAngle
    var maxTranslate = initialDeltaSpace

    var currScore = scorer.score(state)
    var lastScore = Double.NegativeInfinity

    var scoreChange = 0.0                             // how much has the score changed in the last iteration
    var lastScoreChange = Double.NegativeInfinity     // how much had the score changed in the previous interation

    val radius = molA.getRadius + molB.getRadius

    var minTranslationApplied = false


    while (lastScore == Double.NegativeInfinity               // first and second iterations
      || scoreChange > threshold                              // good rate of score increase
      || scoreChange < 0                                      // decayed
      || minTranslationApplied  ) {                           // min translation had to be applied because the molecules were too far apart

      val forces = getForces(molA, molB)  // it is a list of pairs (atomB, force)
      val (translation, m) = getTranslation(forces, maxTranslate)
      minTranslationApplied = m
      val (axis, angle) = getRotation(molB, maxAngle, forces)

      molB.translate(translation)
      log!new Translate(translation)

      molB.rotate(molB.getGeometricCentre, axis, angle)
      log!new Rotate(molB.getGeometricCentre, axis, angle)

      lastScore = currScore
      currScore = scorer.score(state)
      lastScoreChange = scoreChange
      scoreChange = currScore - lastScore

      if (scoreChange < 0){     // decelerate
        maxAngle = maxAngle * 0.9
        maxTranslate = maxTranslate * 0.9
      }

    }
    state
  }

  /** Returns a pair (vector, boolean). The vector has the translation,
    * the boolean indicates if minimum translation had to be applied
    */
  private def getTranslation(atomForces: Seq[(Atom, DenseVector[Double])],
                             maxTranslate: Double) = {
    val forces = atomForces.map{ case (atomB, force) => force}
    val netForce = forces.reduce((a, b) => a + b)
    val totalForceNorm = forces.map(f => norm(f)).sum
    val (translateDistance, minApplied) =
      if (totalForceNorm < minDeltaSpace)
        (minDeltaSpace, true)
      else
        (Math.min(norm(netForce), maxTranslate), false)
    val translation = (netForce * translateDistance) / norm(netForce)
    (translation, minApplied)
  }

  private def getRotation(molB: Molecule, maxAngle: Double,
                          forces: Seq[(Atom, DenseVector[Double])]) = {
    // torque in Euler axis/angle format is cross(r, force) - see http://web.mit.edu/8.01t/www/materials/modules/chapter21.pdf
    val torques = forces.map { case (atomB, force) =>
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
  private def getForces(molA: Molecule, molB: Molecule) = {
    molB.Atoms.map(atomB => (atomB, molToAtomForce(molA, atomB) ))
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
  private def molToAtomForce(molA: Molecule, atomB: Atom): DenseVector[Double] = {
    molA.Atoms.map(atomA => atomToAtomForce(atomA, atomB)).reduce((a, b) => a+b)
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
  private def atomToAtomForce(atomA: Atom, atomB: Atom): DenseVector[Double] = {
    // Atomic force:
    val dif = atomA.coords - atomB.coords;                    // direction from b to a
    val dir = dif / norm(dif);                                // normalized to length 1
    val atomicForceNorm = forceFactor(atomA, atomB)
    val atomicForce = dir * atomicForceNorm

    // electric force:
    val chargeProduct = atomA.partialCharge*atomB.partialCharge
    val electricForceNorm =
      if (chargeProduct < 0)                                  // different charges, attract
        Math.min(atomicForceNorm, 0.0)                        // the atomic force, if distance > optimal
      else if (chargeProduct > 0)                             // equal charges, reject
        Math.max(atomicForceNorm, 0.0)                        // the atomic force, if distance < optimal
      else                                                    // no charge
        0.0
    val electricForce = dir * electricForceNorm

    val bondEnergy = BondEnergy(atomA.element, atomB.element)

    (atomicForce + electricForce) * (bondEnergy / 1000.0)
  }


  /** Calculates a factor f such that the force that a exerts on b
    * is f * a normalized vector pointing from b to a.
    * Positive = attraction, negative = repulsion */
  private def forceFactor(atomA: Atom, atomB: Atom): Double = {
    // This function is similar to SurfaceDistanceScorer, except that on ideal distance it returns 0.
    // The root is 1, so normalize such that optimal distance --> 1

    val optimal = atomA.radius + atomB.radius + 2 * surface // optimal distance
    val actual = atomA.distTo(atomB)
    val dNormalized = actual/optimal
    //expsquare(dNormalized)
    //Math.log(dNormalized)
    explog(dNormalized)
  }

  private def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)

  private def explog(x: Double) = Math.exp(-Math.pow(x,2))*Math.log(x)
}