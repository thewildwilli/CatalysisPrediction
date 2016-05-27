package docking.docksearch

import breeze.linalg.DenseVector
import docking.{Docker, DockingState}
import docking.dockscore.Scorer
import model.{Atom, Molecule}
import opt.HillClimbing

/**   This docker places molecule B in six initial positions with respect
  *   to molecule A: to the right, left, on the front, on the back, over, under.
  *   Then uses hill climbing from each of these to find a solution.
  *   Hill climbing needs to both rotate and translate molecule B.
  */
object SixInitialsDocker extends Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    val aRadius = molA.getRadius
    val bRadius = molB.getRadius

    val right = molB.clone; right.translate(DenseVector(aRadius + bRadius, 0.0, 0.0))
    val left = molB.clone; left.translate(DenseVector(-aRadius - bRadius, 0.0, 0.0))
    val over = molB.clone; over.translate(DenseVector(0.0, aRadius + bRadius, 0.0))
    val under = molB.clone; under.translate(DenseVector(0.0, -aRadius - bRadius, 0.0))
    val behind = molB.clone; behind.translate(DenseVector(0.0, 0.0, aRadius + bRadius))
    val before = molB.clone; before.translate(DenseVector(0.0, 0.0, -aRadius - bRadius))

    List(right, left, over, under, behind, before).map(b => performDock(molA, b, scorer)).maxBy(scorer.score)
  }

  private def performDock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    val initialState = new AllNeighboursState(molA, molB)
    print("performDock...")
    HillClimbing.optimize(initialState, 1000, scorer.score).asInstanceOf[AllNeighboursState]
  }

}

class AllNeighboursState(molA: Molecule, molB: Molecule) extends DockingState(molA, molB) {

  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians
  final val DeltaSpace = 0.1 // Armstrongs

  /** Neighbours are all rotations of molecule bMol around atom bAtom
    * by DeltaAngle in both directions in all 3 axis, plus all translations.
    * That is, each state has 12 neighbours.
    * Rotations are done around b's geometric centre.
    * Molecule a is fixed.
    */
  override def getNeighbours = {
    val centre = molB.getGeometricCentre
    val fwX = molB.clone; fwX.rotateX(centre, DeltaAngle)   // forward rotation on X axis
    val bwX = molB.clone; bwX.rotateX(centre, -DeltaAngle)  // backward rotation on X axis
    val fwY = molB.clone; fwY.rotateY(centre, DeltaAngle)   // forward rotation on Y axis
    val bwY = molB.clone; bwY.rotateY(centre, -DeltaAngle)  // backward rotation on Y axis
    val fwZ = molB.clone; fwZ.rotateZ(centre, DeltaAngle)   // forward rotation on Z axis
    val bwZ = molB.clone; bwZ.rotateZ(centre, -DeltaAngle)  // backward rotation on Z axis

    val fwtX = molB.clone; fwtX.translate(DenseVector(DeltaSpace, 0.0, 0.0))   // forward translation on X axis
    val bwtX = molB.clone; bwtX.translate(DenseVector(-DeltaSpace, 0.0, 0.0))  // backward translation on X axis
    val fwtY = molB.clone; fwtY.translate(DenseVector(0.0, DeltaSpace, 0.0))   // forward translation on Y axis
    val bwtY = molB.clone; bwtY.translate(DenseVector(0.0, -DeltaSpace, 0.0))  // backward translation on Y axis
    val fwtZ = molB.clone; fwtZ.translate(DenseVector(0.0, 0.0, DeltaSpace))   // forward translation on Z axis
    val bwtZ = molB.clone; bwtZ.translate(DenseVector(0.0, 0.0, -DeltaSpace))  // backward translation on Z axis
    List(fwX, bwX, fwY, bwY, fwZ, bwZ, fwtX, bwtX, fwtY, bwtY, fwtZ, bwtZ).map(b => new AllNeighboursState(molA, b))
  }
}
