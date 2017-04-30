package docking.docksearch

import breeze.linalg.DenseVector
import docking._
import docking.dockscore.Scorer
import model.{Molecule, Rotate, Translate}
import opt.{EnhHillClimbing}

/**   This docker places molecule B in for initial positions with respect
  *   to molecule A: to the right, left, over, under.
  *   Then uses hill climbing from each of these to find a solution.
  *   Hill climbing needs to both rotate and translate molecule B.
  */
class FourInitialsDocker2D(val scorer: Scorer) extends Docker {
  var InitialDeltaAngle = Math.toRadians(20) // 20 degrees in radians
  var InitialDeltaSpace = 1.0

  def dock(molA: Molecule, molB: Molecule, log: DockLog) = {
    val aRadius = molA.getRadius
    val bRadius = molB.getRadius

    val right = DenseVector(aRadius + bRadius, 0.0, 0.0)
    val left = DenseVector(-aRadius - bRadius, 0.0, 0.0)
    val over = DenseVector(0.0, aRadius + bRadius, 0.0)
    val under = DenseVector(0.0, -aRadius - bRadius, 0.0)

    val bestState = List(over, under, left, right).map(pos =>
      performDock(molA, molB, pos, scorer, log)).maxBy(scorer.score)

    (bestState.b, scorer.score(bestState))
  }

  private def performDock(molA: Molecule, molB: Molecule, startPos: DenseVector[Double],
                          scorer: Scorer, log: DockLog): DockingState = {
    log.reset

    // move molB to its starting position
    val t = new Translate(startPos)
    val initState = DockingState.transition(new DockingState(molA, molB), t)
    log.action(t)

    var deltaAngle = InitialDeltaAngle
    var deltaSpace = InitialDeltaSpace

    /* Neighbours are translations of molB and rotations of molB around
      * its own centre   */
    EnhHillClimbing.optimize[DockingState](initState, (s) => {
      val centre = s.b.getGeometricCentre
      List(
        new Rotate(centre, DenseVector(0.0, 0, 1),  deltaAngle),   // forward rotation on Z axis
        new Rotate(centre, DenseVector(0.0, 0, 1), -deltaAngle),    // backward rotation on Z axis

        new Translate(DenseVector(deltaSpace, 0.0, 0.0)),                  // forward translation on X axis
        new Translate(DenseVector(-deltaSpace, 0.0, 0.0)),                 // backward translation on X axis
        new Translate(DenseVector(0.0, deltaSpace, 0.0)),                  // forward translation on Y axis
        new Translate(DenseVector(0.0, -deltaSpace, 0.0))                 // backward translation on Y axis
      )
    }, DockingState.transition, () => {
      deltaAngle = deltaAngle / 2.0
      deltaSpace = deltaSpace / 2.0
    }, scorer.score, 50, log)
  }

}