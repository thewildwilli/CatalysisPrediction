package docking.docksearch

import breeze.linalg.DenseVector
import docking.{DockLog, Docker, DockingState}
import docking.dockscore.Scorer
import model.{Molecule, Rotate, Translate}
import opt.HillClimbing

// Created by Ernesto on 23/05/2016.
class AtomPairDocker(val scorer: Scorer) extends Docker {
  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians

  def dock(molA: Molecule, molB: Molecule, log: DockLog) = {
    var maxScore = Double.NegativeInfinity
    var bestMatch = null.asInstanceOf[Molecule]
    var i=0

    for (x <- molA.atoms; y <- molB.atoms) {
      val optimized = dockPair2D(molA, x.id, molB, y.id, scorer, log)
      val score = scorer.score(optimized)
      if (score > maxScore) {
        maxScore = score
        bestMatch = optimized.b
      }

      printf("docked %d of %d pairs, best score: %.4f %n", {i+=1;i}, molA.atoms.size*molB.atoms.size, maxScore)
    }
    (bestMatch, maxScore)
  }

  private def dockPair(molA: Molecule, x: Int, molB: Molecule, y: Int,
                       scorer: Scorer, log: DockLog) = {
    log.reset

    //    translate molB so that the pair (atomA, atomB) overlaps:
    val t = new Translate(molA(x).coords - molB(y).coords)
    val initState = DockingState.transition(new DockingState(molA, molB), t)
    log.action(t)

    //    call hill climbing - neigbouring states are rotations
    HillClimbing.optimize[DockingState](initState, (s) => {
      val atomBCoords = s.b(y).coords
      List(
        new Rotate(atomBCoords, DenseVector(1.0, 0, 0),  DeltaAngle),   // forward rotation on X axis
        new Rotate(atomBCoords, DenseVector(1.0, 0, 0), -DeltaAngle),   // backward rotation on X axis
        new Rotate(atomBCoords, DenseVector(0.0, 1, 0),  DeltaAngle),   // forward rotation on Y axis
        new Rotate(atomBCoords, DenseVector(0.0, 1, 0), -DeltaAngle),   // backward rotation on Y axis
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1),  DeltaAngle),   // forward rotation on Z axis
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1), -DeltaAngle)    // backward rotation on Z axis
      )
    }, DockingState.transition, scorer.score, 50, log)
  }

  private def dockPair2D(molA: Molecule, x: Int, molB: Molecule, y: Int,
                         scorer: Scorer, log: DockLog) = {
    log.reset

    //    translate molB so that the pair (atomA, atomB) overlaps:
    val t = new Translate(molA(x).coords - molB(y).coords)
    val initState = DockingState.transition(new DockingState(molA, molB), t)
    log.action(t)

    //    call hill climbing - neigbouring states are rotations
    HillClimbing.optimize[DockingState](initState, (s) => {
      val atomBCoords = s.b(y).coords
      List(
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1),  DeltaAngle),   // forward rotation on Z axis
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1), -DeltaAngle)    // backward rotation on Z axis
      )
    }, DockingState.transition, scorer.score, 50, log)
  }
}