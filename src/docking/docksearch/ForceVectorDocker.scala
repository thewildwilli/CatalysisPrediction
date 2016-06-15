package docking.docksearch

import breeze.linalg.{DenseVector, norm}
import docking.{Docker, DockingState}
import docking.dockscore.Scorer
import model.{Atom, Molecule}

// Created by Ernesto on 08/06/2016.
// FIX THE DESIGN A BIT!!!
class ForceVectorDocker(val surface: Double) {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer) = {
    // just to test
    forceVectorDock(molA, molB)
    molB
  }

  /**  Docks b into a from b's initial position and orientation, using force vectors
    */
  private def forceVectorDock(molA: Molecule, molB: Molecule) = {
    // get the force of each atom
    // translate a rotate accordingly
    // for a test, just translate the entire molecule
    for (i <- 0 until 50) {
      val molBForce = molB.Atoms.map(atomB => molToAtomForce(molA, atomB)).reduce((a, b) => a + b)
      molB.translate(molBForce)
    }
  }

  /** Calculates the total force that molA exerts on atomB: the sum
    * of the forces each atom in a exerts on atomB   */
  private def molToAtomForce(molA: Molecule, atomB: Atom): DenseVector[Double] = {
    molA.Atoms.map(atomA => atomToAtomForce(atomA, atomB)).reduce((a, b) => a+b);
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

    val radii = atomA.radius + atomB._radius + 2 * surface // optimal distance
    val d = atomA.distTo(atomB)
    val dNormalized = d/radii
    expsquare(dNormalized)
  }

  private def expsquare(x: Double) = Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)
}
