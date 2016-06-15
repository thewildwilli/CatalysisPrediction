package docking.docksearch

import breeze.linalg.DenseVector
import docking.{Decaying, Docker, DockingState}
import docking.dockscore.Scorer
import model.Molecule
import opt.EnhHillClimbing

import scala.collection.mutable.ListBuffer

/**   This docker places molecule B in six initial positions with respect
  *   to molecule A: to the right, left, on the front, on the back, over, under.
  *   Then uses hill climbing from each of these to find a solution.
  *   Hill climbing needs to both rotate and translate molecule B.
  */
object SixInitialsDocker extends Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
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
      val c = i.getGeometricCentre; val ang = Math.toRadians(90);
      initials.append(i); initials.append(i.clone.rotateX(c, ang)); initials.append(i.clone.rotateX(c, 2*ang)); initials.append(i.clone.rotateX(c, 3*ang));
      initials.append(i.clone.rotateY(c, ang)); initials.append(i.clone.rotateY(c, 2*ang)); initials.append(i.clone.rotateY(c, 3*ang));
      initials.append(i.clone.rotateZ(c, ang)); initials.append(i.clone.rotateZ(c, 2*ang)); initials.append(i.clone.rotateZ(c, 3*ang));
    }

    initials.map(b => performDock(molA, b, scorer)).maxBy(scorer.score)
  }

  private def performDock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    val initialState = new AllNeighboursState(molA, molB)
    print("performDock...")
    EnhHillClimbing.optimize(initialState, 1000, scorer.score).asInstanceOf[AllNeighboursState]
  }

}

class AllNeighboursState(molA: Molecule, molB: Molecule,
                          var deltaAngle: Double = Math.toRadians(20), var deltaSpace: Double = 0.1)
                          extends DockingState(molA, molB) with Decaying {

  /** Neighbours are all rotations of molecule bMol around atom bAtom
    )* by DeltaAngle in both directions in all 3 axis, plus all translations.
    * That is, each state has 12 neighbours.
    * Rotations are done around b's geometric centre.
    * Molecule a is fixed.
    */
  override def getNeighbours = {
    val centre = molB.getGeometricCentre
    val fwX = molB.clone; fwX.rotateX(centre, deltaAngle)   // forward rotation on X axis
    val bwX = molB.clone; bwX.rotateX(centre, -deltaAngle)  // backward rotation on X axis
    val fwY = molB.clone; fwY.rotateY(centre, deltaAngle)   // forward rotation on Y axis
    val bwY = molB.clone; bwY.rotateY(centre, -deltaAngle)  // backward rotation on Y axis
    val fwZ = molB.clone; fwZ.rotateZ(centre, deltaAngle)   // forward rotation on Z axis
    val bwZ = molB.clone; bwZ.rotateZ(centre, -deltaAngle)  // backward rotation on Z axis

    val fwtX = molB.clone; fwtX.translate(DenseVector(deltaSpace, 0.0, 0.0))   // forward translation on X axis
    val bwtX = molB.clone; bwtX.translate(DenseVector(-deltaSpace, 0.0, 0.0))  // backward translation on X axis
    val fwtY = molB.clone; fwtY.translate(DenseVector(0.0, deltaSpace, 0.0))   // forward translation on Y axis
    val bwtY = molB.clone; bwtY.translate(DenseVector(0.0, -deltaSpace, 0.0))  // backward translation on Y axis
    val fwtZ = molB.clone; fwtZ.translate(DenseVector(0.0, 0.0, deltaSpace))   // forward translation on Z axis
    val bwtZ = molB.clone; bwtZ.translate(DenseVector(0.0, 0.0, -deltaSpace))  // backward translation on Z axis
    List(fwX, bwX, fwY, bwY, fwZ, bwZ, fwtX, bwtX, fwtY, bwtY, fwtZ, bwtZ).map(b =>
      new AllNeighboursState(molA, b, deltaAngle, deltaSpace))
  }

  /** Decay parameters by half */
  override def decayRate: Unit = {
    println("decaying...")
    deltaAngle = deltaAngle / 2.0
    deltaSpace = deltaSpace / 2.0
  }
}
