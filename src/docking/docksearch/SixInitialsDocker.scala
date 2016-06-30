package docking.docksearch

import breeze.linalg.DenseVector
import docking._
import docking.dockscore.Scorer
import io.threadcso._
import model.Molecule
import opt.{Action, EnhHillClimbing, HillClimbing}

import scala.collection.mutable.ListBuffer

/**   This docker places molecule B in six initial positions with respect
  *   to molecule A: to the right, left, on the front, on the back, over, under.
  *   Then uses hill climbing from each of these to find a solution.
  *   Hill climbing needs to both rotate and translate molecule B.
  */
object SixInitialsDocker extends Docker {
  var InitialDeltaAngle = Math.toRadians(20) // 20 degrees in radians
  var InitialDeltaSpace = 0.1

  def dock(molA: Molecule, molB: Molecule, scorer: Scorer, log: ![Any]): DockingState = {
    val aRadius = molA.getRadius
    val bRadius = molB.getRadius

    // initial positions:
    val right = molB.clone; right.translate(DenseVector(aRadius + bRadius, 0.0, 0.0))
    val left = molB.clone; left.translate(DenseVector(-aRadius - bRadius, 0.0, 0.0))
    val over = molB.clone; over.translate(DenseVector(0.0, aRadius + bRadius, 0.0))
    val under = molB.clone; under.translate(DenseVector(0.0, -aRadius - bRadius, 0.0))
    val behind = molB.clone; behind.translate(DenseVector(0.0, 0.0, aRadius + bRadius))
    val before = molB.clone; before.translate(DenseVector(0.0, 0.0, -aRadius - bRadius))

    // for each initial position, get 10 initial orientations: 4 rotations of 90 degrees in each of the 3 axis
    val initials = new ListBuffer[Molecule]()
    for(i <- List(right, left, over, under, behind, before)) {
      val c = i.getGeometricCentre; val ang = Math.toRadians(90)
      initials.append(i); initials.append(i.clone.rotateX(c, ang)); initials.append(i.clone.rotateX(c, 2*ang)); initials.append(i.clone.rotateX(c, 3*ang))
      initials.append(i.clone.rotateY(c, ang)); initials.append(i.clone.rotateY(c, 2*ang)); initials.append(i.clone.rotateY(c, 3*ang))
      initials.append(i.clone.rotateZ(c, ang)); initials.append(i.clone.rotateZ(c, 2*ang)); initials.append(i.clone.rotateZ(c, 3*ang))
    }

    initials.map(b => performDock(molA, b, scorer)).maxBy(scorer.score)
  }

  private def performDock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    print("performDock...")

    var deltaAngle = InitialDeltaAngle
    var deltaSpace = InitialDeltaSpace

    /* Neighbours are translations of molB and rotations of molB around
      * its own centre   */
    EnhHillClimbing.optimize[DockingState](new DockingState(molA, molB), (s) => {
      val centre = s.b.getGeometricCentre
      List(
        new Rotate(centre, DenseVector(1.0, 0, 0, 1),  deltaAngle),   // forward rotation on X axis
        new Rotate(centre, DenseVector(1.0, 0, 0, 1), -deltaAngle),   // backward rotation on X axis
        new Rotate(centre, DenseVector(0.0, 1, 0, 1),  deltaAngle),   // forward rotation on Y axis
        new Rotate(centre, DenseVector(0.0, 1, 0, 1), -deltaAngle),   // backward rotation on Y axis
        new Rotate(centre, DenseVector(0.0, 0, 1, 1),  deltaAngle),   // forward rotation on Z axis
        new Rotate(centre, DenseVector(0.0, 0, 1, 1), -deltaAngle),    // backward rotation on Z axis

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
    }, scorer.score, 50, null)
  }

}

class AllNeighboursState(molA: Molecule, molB: Molecule,
                          var deltaAngle: Double = Math.toRadians(20), var deltaSpace: Double = 0.1)
                          extends DockingState(molA, molB) with Decaying {

  /** Decay parameters by half */
  override def decayRate: Unit = {
    println("decaying...")
    deltaAngle = deltaAngle / 2.0
    deltaSpace = deltaSpace / 2.0
  }
}
