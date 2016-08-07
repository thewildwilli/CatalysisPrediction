package docking.docksearch

import io.threadcso._

import breeze.linalg.DenseVector
import docking.{Docker, DockingState, Rotate, Translate}
import docking.dockscore.Scorer
import model.Molecule
import opt.HillClimbing


// Created by Ernesto on 23/05/2016.
class SurfaceAtomPairsWithTranslationDocker(val scorer: Scorer) extends Docker {

  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians
  final val DeltaSpace = 0.1 // Angstrongs

  def dock(molA: Molecule, molB: Molecule, log: ![Any] = null) = {
    var maxScore = Double.NegativeInfinity
    var bestMatch = null.asInstanceOf[Molecule]
    var i=0

    for (x <- molA.Atoms.indices; y <- molB.Atoms.indices) {
      if (molA(x).isSurface && molB(y).isSurface) {
        val optimized = dockPair(molA, x, molB, y, scorer)
        val score = scorer.score(optimized)
        if (score > maxScore) {
          maxScore = score
          bestMatch = optimized.b
        }

        printf("docked %d of %d pairs, best score: %.4f %n", {
          i += 1
          i
        }, molA.Atoms.count(a => a.isSurface) * molB.Atoms.count(a => a.isSurface), maxScore)
      }
    }
    (bestMatch, maxScore)
  }

  private def dockPair2D(molA: Molecule, x: Int, molB: Molecule, y: Int, scorer: Scorer) = {
    //    translate molB so that the pair (atomA, atomB) overlaps
    molB.translate(molA(x).coords - molB(y).coords)

    /** Neighbours are all rotations of molecule bMol around atom bAtom
      * by DeltaAngle in both directions in all 3 axis. That is, each
      * state has 6 neighbours. Molecule a is fixed.
      */
    HillClimbing.optimize[DockingState](new DockingState(molA, molB), (s) => {
      val atomBCoords = s.b(y).coords
      List(
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1, 1),  DeltaAngle),   // forward rotation on Z axis
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1, 1), -DeltaAngle),   // backward rotation on Z axis

        new Translate(DenseVector(DeltaSpace, 0.0, 0.0)),                  // forward translation on X axis
        new Translate(DenseVector(-DeltaSpace, 0.0, 0.0)),                 // backward translation on X axis
        new Translate(DenseVector(0.0, DeltaSpace, 0.0)),                  // forward translation on Y axis
        new Translate(DenseVector(0.0, -DeltaSpace, 0.0))                  // backward translation on Y axis
      )
    }, DockingState.transition, scorer.score, 50)
  }

  private def dockPair(molA: Molecule, x: Int, molB: Molecule, y: Int, scorer: Scorer) = {
    //    translate molB so that the pair (atomA, atomB) overlaps
    molB.translate(molA(x).coords - molB(y).coords)

    /** Neighbours are all rotations of molecule bMol around atom bAtom
      * by DeltaAngle in both directions in all 3 axis. That is, each
      * state has 6 neighbours. Molecule a is fixed.
      */
    HillClimbing.optimize[DockingState](new DockingState(molA, molB), (s) => {
      val atomBCoords = s.b(y).coords
      List(
        new Rotate(atomBCoords, DenseVector(1.0, 0, 0, 1),  DeltaAngle),   // forward rotation on X axis
        new Rotate(atomBCoords, DenseVector(1.0, 0, 0, 1), -DeltaAngle),   // backward rotation on X axis
        new Rotate(atomBCoords, DenseVector(0.0, 1, 0, 1),  DeltaAngle),   // forward rotation on Y axis
        new Rotate(atomBCoords, DenseVector(0.0, 1, 0, 1), -DeltaAngle),   // backward rotation on Y axis
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1, 1),  DeltaAngle),   // forward rotation on Z axis
        new Rotate(atomBCoords, DenseVector(0.0, 0, 1, 1), -DeltaAngle),   // backward rotation on Z axis

        new Translate(DenseVector(DeltaSpace, 0.0, 0.0)),                  // forward translation on X axis
        new Translate(DenseVector(-DeltaSpace, 0.0, 0.0)),                 // backward translation on X axis
        new Translate(DenseVector(0.0, DeltaSpace, 0.0)),                  // forward translation on Y axis
        new Translate(DenseVector(0.0, -DeltaSpace, 0.0)),                 // backward translation on Y axis
        new Translate(DenseVector(0.0, 0.0, DeltaSpace)),                  // forward translation on Y axis
        new Translate(DenseVector(0.0, 0.0, -DeltaSpace))                  // backward translation on Y axis
      )
    }, DockingState.transition, scorer.score, 50)
  }
}
