package docking.docksearch

import breeze.linalg.DenseVector
import docking._
import docking.dockscore.Scorer
import io.threadcso._
import model.{Molecule, Rotate, Translate}
import opt.{EnhHillClimbing}

/**   This docker places molecule B in six initial positions with respect
  *   to molecule A: to the right, left, on the front, on the back, over, under.
  *   Then uses hill climbing from each of these to find a solution.
  *   Hill climbing needs to both rotate and translate molecule B.
  */
class EhcDocker(val scorer: Scorer, val maxIters: Int) extends Docker {
  var InitialDeltaAngle = Math.toRadians(20) // 20 degrees in radians
  var InitialDeltaSpace = 1.0

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]) = {
    var deltaAngle = InitialDeltaAngle
    var deltaSpace = InitialDeltaSpace

    /* Neighbours are translations of molB and rotations of molB around
      * its own centre   */
    val finalState = EnhHillClimbing.optimize[DockingState](
      new DockingState(molA, molB), (s) => {
        val centre = s.b.getGeometricCentre
        List(
          new Rotate(centre, DenseVector(1.0, 0, 0),  deltaAngle),   // forward rotation on X axis
          new Rotate(centre, DenseVector(1.0, 0, 0), -deltaAngle),   // backward rotation on X axis
          new Rotate(centre, DenseVector(0.0, 1, 0),  deltaAngle),   // forward rotation on Y axis
          new Rotate(centre, DenseVector(0.0, 1, 0), -deltaAngle),   // backward rotation on Y axis
          new Rotate(centre, DenseVector(0.0, 0, 1),  deltaAngle),   // forward rotation on Z axis
          new Rotate(centre, DenseVector(0.0, 0, 1), -deltaAngle),    // backward rotation on Z axis

          new Translate(DenseVector(deltaSpace, 0.0, 0.0)),                  // forward translation on X axis
          new Translate(DenseVector(-deltaSpace, 0.0, 0.0)),                 // backward translation on X axis
          new Translate(DenseVector(0.0, deltaSpace, 0.0)),                  // forward translation on Y axis
          new Translate(DenseVector(0.0, -deltaSpace, 0.0)),                 // backward translation on Y axis
          new Translate(DenseVector(0.0, 0.0, deltaSpace)),                  // forward translation on Z axis
          new Translate(DenseVector(0.0, 0.0, -deltaSpace))                  // backward translation on Z axis
        )
    }, DockingState.transition, () => {
      deltaAngle = deltaAngle / 2.0
      deltaSpace = deltaSpace / 2.0
    }, scorer.score, maxIters, log)
    (finalState.b, scorer.score(finalState))
  }

}
