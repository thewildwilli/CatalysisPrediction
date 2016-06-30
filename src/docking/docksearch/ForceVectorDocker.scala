package docking.docksearch

import breeze.linalg
import breeze.linalg._

import docking.dockscore.Scorer
import docking._
import opt._
import io.threadcso._
import model.{Atom, Molecule}

// Created by Ernesto on 08/06/2016.
// FIX THE DESIGN A BIT!!!
class ForceVectorDocker(val surface: Double, val maxDecays: Int = 10) extends Docker {
  val initialDeltaAngle = Math.toRadians(20) // 20 degrees in radians
  val initialDeltaSpace = 1.0

  override def dock(molA: Molecule, molB: Molecule, scorer: Scorer, log: ![Any]) = {

    val initialDistance = Math.max(molA.getRadius, molB.getRadius)

    val right = DenseVector(initialDistance, 0.0, 0.0)
    val left = DenseVector(-initialDistance, 0.0, 0.0)
    val over = DenseVector(0.0, initialDistance, 0.0)
    val under = DenseVector(0.0, -initialDistance, 0.0)

    List(over, under, left, right).map(pos => {
      log!Reset
      val b = molB.clone

      b.translate(pos)
      log!new Translate(pos)

      forceVectorDock(molA, b, scorer, log)
    }).maxBy(scorer.score)
  }

  /**  Docks b into a from b's initial position and orientation, using force vectors
    */
  private def forceVectorDock(molA: Molecule, molB: Molecule,
                              scorer: Scorer, log: ![Action]) = {

    // get the force of each atom
    // translate a rotate accordingly
    // for a test, just translate the entire molecule

    var deltaAngle = initialDeltaAngle
    var deltaSpace = initialDeltaSpace
    var lastScore = Double.NegativeInfinity
    val state = new DockingState(molA, molB)

    var decays = 0
    while (decays <= maxDecays) {
      val forces = getForces(molA, molB);  // it is a list of pairs (atomB, force)

      val netForce = forces.map{ case (atomB, force) => force}.reduce((a, b) => a + b)
      val translation = (netForce * deltaSpace) / norm(netForce)

      // torque in Euler axis/angle format is cross(r, force) - see http://web.mit.edu/8.01t/www/materials/modules/chapter21.pdf
      val torques = forces.map { case (atomB, force) =>
        val r = atomB.coords - molB.getGeometricCentre // radius vector from centre to atom
        linalg.cross(r, force)
      }
      val netTorque = torques.reduce((a, b) => a + b)     // add all torques together
      val axis = netTorque / norm(netTorque) // normalize

      molB.translate(translation)
      log!new Translate(translation)

      molB.rotate(molB.getGeometricCentre, axis, deltaAngle)
      log!new Rotate(molB.getGeometricCentre, axis, deltaAngle)

      val score = scorer.score(state)
      if (score < lastScore){
        deltaAngle = deltaAngle / 2.0
        deltaSpace = deltaSpace / 2.0
        decays += 1
      }
      lastScore = score
    }
    state
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

  /** calculates the force that atomA excerts on atomB */
  private def atomToAtomForce(atomA: Atom, atomB: Atom): DenseVector[Double] = {
    val dir = atomA.coords - atomB.coords;      // direction from b to a
    val normdir = dir / norm(dir);              // normalized to length 1
    normdir * forceFactor(atomA, atomB)
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
    expsquare(dNormalized)
  }

  private def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)
}
