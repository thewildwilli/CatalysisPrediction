package docking.docksearch

import breeze.linalg.DenseVector
import docking.{Docker, DockingState}
import docking.dockscore.Scorer
import model.{Atom, Molecule}
import opt.{HillClimbing, State}

// Created by Ernesto on 23/05/2016.
object SurfaceAtomPairsWithTranslationDocker extends Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    var maxScore = Double.NegativeInfinity
    var bestMatch = null.asInstanceOf[DockingState]
    var i=0;

    for (x <- 0 until molA.Atoms.size; y <- 0 until molB.Atoms.size) {
      if (molA(x).isSurface && molB(y).isSurface) {
        val optimized = dockPair(molA, x, molB, y, scorer)
        val score = scorer.score(optimized)
        if (score > maxScore) {
          maxScore = score
          bestMatch = optimized
        }

        printf("docked %d of %d pairs, best score: %.4f %n", {
          i += 1;
          i
        }, molA.Atoms.filter(a => a.isSurface).size * molB.Atoms.filter(a => a.isSurface).size, maxScore);
      }
    }
    bestMatch
  }

  private def dockPair(molA: Molecule, x: Int, molB: Molecule, y: Int, scorer: Scorer) = {
    //    translate molB so that the pair (atomA, atomB) overlaps
    molB.translate(molA(x).coords - molB(y).coords)

    //    call hill climbing - neighbouring states are rotations
    val initialState = new SurfaceAtomPairsWithTranslationDocker.AtomPairState2D(molA, x, molB, y)
    HillClimbing.optimize(initialState, 50, scorer.score).asInstanceOf[SurfaceAtomPairsWithTranslationDocker.AtomPairState2D]
  }


  class AtomPairState(molA: Molecule, x: Int, molB: Molecule, y: Int) extends DockingState(molA, molB) {

    final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians
    final val DeltaSpace = 0.1 // Angstroms

    /** Neighbours are all rotations of molecule bMol around atom bAtom
      * by DeltaAngle in both directions in all 3 axis. That is, each
      * state has 6 neighbours. Molecule a is fixed.
      */
    override def getNeighbours = {
      val centre = molB(y).coords
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
      List(fwX, bwX, fwY, bwY, fwZ, bwZ, fwtX, bwtX, fwtY, bwtY, fwtZ, bwtZ).map(b =>
        new SurfaceAtomPairsWithTranslationDocker.AtomPairState(molA, x, b, y))
    }
  }

  class AtomPairState2D(molA: Molecule, x: Int, molB: Molecule, y: Int) extends DockingState(molA, molB) {

    final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians
    final val DeltaSpace = 0.1 // Angstrongs

    /** Neighbours are all rotations of molecule bMol around atom bAtom
      * by DeltaAngle in both directions in all 3 axis. That is, each
      * state has 6 neighbours. Molecule a is fixed.
      */
    override def getNeighbours = {
      val centre = molB(y).coords
      val fwZ = molB.clone; fwZ.rotateZ(centre, DeltaAngle)   // forward rotation on Z axis
      val bwZ = molB.clone; bwZ.rotateZ(centre, -DeltaAngle)  // backward rotation on Z axis

      val fwtX = molB.clone; fwtX.translate(DenseVector(DeltaSpace, 0.0, 0.0))   // forward translation on X axis
      val bwtX = molB.clone; bwtX.translate(DenseVector(-DeltaSpace, 0.0, 0.0))  // backward translation on X axis
      val fwtY = molB.clone; fwtY.translate(DenseVector(0.0, DeltaSpace, 0.0))   // forward translation on Y axis
      val bwtY = molB.clone; bwtY.translate(DenseVector(0.0, -DeltaSpace, 0.0))  // backward translation on Y axis

      List(fwZ, bwZ, fwtX, bwtX, fwtY, bwtY).map(b =>
        new SurfaceAtomPairsWithTranslationDocker.AtomPairState2D(molA, x, b, y))
    }
  }

}
