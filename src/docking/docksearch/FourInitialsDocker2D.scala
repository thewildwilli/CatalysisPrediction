package docking.docksearch

import breeze.linalg.DenseVector
import docking.{Decaying, Docker, DockingState}
import docking.dockscore.Scorer
import model.{Molecule}
import opt.{EnhHillClimbing}

/**   This docker places molecule B in for initial positions with respect
  *   to molecule A: to the right, left, over, under.
  *   Then uses hill climbing from each of these to find a solution.
  *   Hill climbing needs to both rotate and translate molecule B.
  */
object FourInitialsDocker2D extends Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    val aRadius = molA.getRadius
    val bRadius = molB.getRadius

    val right = molB.clone; right.translate(DenseVector(aRadius + bRadius, 0.0, 0.0))
    val left = molB.clone; left.translate(DenseVector(-aRadius - bRadius, 0.0, 0.0))
    val over = molB.clone; over.translate(DenseVector(0.0, aRadius + bRadius, 0.0))
    val under = molB.clone; under.translate(DenseVector(0.0, -aRadius - bRadius, 0.0))

    List(right, left, over, under).map(b => performDock(molA, b, scorer)).maxBy(scorer.score)
  }

  private def performDock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    val initialState = new AllNeighboursState2D(molA, molB)
    print("performDock...")
    EnhHillClimbing.optimize(initialState, 1000, scorer.score).asInstanceOf[AllNeighboursState2D]
  }

}


/** DefaultAngle: in radians. DefaultSpace: in Angstroms.
 */
class AllNeighboursState2D(molA: Molecule, molB: Molecule,
                           var deltaAngle: Double = Math.toRadians(20), var deltaSpace: Double = 0.1)
                           extends DockingState(molA, molB) with Decaying {

  /** Neighbours are all rotations of molecule bMol around atom bAtom
    * by DeltaAngle in both directions in the Z axis, plus all translations in X and Y.
    * That is, each state has 6 neighbours.
    * Rotations are done around b's geometric centre.
    * Molecule a is fixed.
    */
  override def getNeighbours = {
    val centre = molB.getGeometricCentre
    val fwZ = molB.clone; fwZ.rotateZ(centre, deltaAngle)   // forward rotation on Z axis
    val bwZ = molB.clone; bwZ.rotateZ(centre, -deltaAngle)  // backward rotation on Z axis

    val fwtX = molB.clone; fwtX.translate(DenseVector(deltaSpace, 0.0, 0.0))   // forward translation on X axis
    val bwtX = molB.clone; bwtX.translate(DenseVector(-deltaSpace, 0.0, 0.0))  // backward translation on X axis
    val fwtY = molB.clone; fwtY.translate(DenseVector(0.0, deltaSpace, 0.0))   // forward translation on Y axis
    val bwtY = molB.clone; bwtY.translate(DenseVector(0.0, -deltaSpace, 0.0))  // backward translation on Y axis
    List(fwtX, bwtX, fwtY, bwtY,fwZ, bwZ).map(b => new AllNeighboursState2D(molA, b, deltaAngle, deltaSpace))
  }

  /** Decay parameters by half */
  override def decayRate: Unit = {
    println("decaying...")
    deltaAngle = deltaAngle / 2.0
    deltaSpace = deltaSpace / 2.0
  }
}
