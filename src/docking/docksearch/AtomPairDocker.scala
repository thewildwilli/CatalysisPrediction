package docking.docksearch

import docking.DockingState
import docking.dockscore.Scorer
import model.{Atom, Molecule}
import opt.{HillClimbing}

// Created by Ernesto on 23/05/2016.
object AtomPairDocker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    var maxScore = Double.NegativeInfinity
    var bestMatch = null.asInstanceOf[AtomPairState]
    var i=0;

    for (atomA <- molA.Atoms; atomB <- molB.Atoms) {
      val optimized = dockPair(molA, atomA, molB, atomB, scorer)
      val score = scorer.score(optimized)
      if (score > maxScore) {
        maxScore = score
        bestMatch = optimized
      }

      printf("docked %d of %d pairs, best score: %.4f %n", {i+=1;i}, molA.Atoms.size*molB.Atoms.size, maxScore);
    }
    bestMatch
  }

  private def dockPair(molA: Molecule, atomA: Atom, molB: Molecule, atomB: Atom, scorer: Scorer) = {
    //    translate molB so that the pair (atomA, atomB) overlaps
    molB.translate(atomA.coords - atomB.coords)

    //    call hill climbing - neigbouring states are rotations
    val initialState = new AtomPairState(molA, atomA, molB, atomB)
    HillClimbing.optimize(initialState, 50, scorer.score).asInstanceOf[AtomPairState]
  }
}

class AtomPairState(molA: Molecule, atomA: Atom, molB: Molecule, atomB: Atom) extends DockingState(molA, molB) {

  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians

  /** Neighbours are all rotations of molecule bMol around atom bAtom
    * by DeltaAngle in both directions in all 3 axis. That is, each
    * state has 6 neighbours. Molecule a is fixed.
    */
  override def getNeighbours = {
    val fwX = molB.clone; fwX.rotateX(atomB.coords, DeltaAngle)   // forward rotation on X axis
    val bwX = molB.clone; bwX.rotateX(atomB.coords, -DeltaAngle)  // backward rotation on X axis
    val fwY = molB.clone; fwY.rotateY(atomB.coords, DeltaAngle)   // forward rotation on Y axis
    val bwY = molB.clone; bwY.rotateY(atomB.coords, -DeltaAngle)  // backward rotation on Y axis
    val fwZ = molB.clone; fwZ.rotateZ(atomB.coords, DeltaAngle)   // forward rotation on Z axis
    val bwZ = molB.clone; bwZ.rotateZ(atomB.coords, -DeltaAngle)  // backward rotation on Z axis
    List(fwX, bwX, fwY, bwY, fwZ, bwZ).map(b => new AtomPairState(molA, atomA, b, atomB))
    // WARNING: atomB is NOT from the cloned molecule
  }
}


class AtomPairState2D(molA: Molecule, atomA: Atom, molB: Molecule, atomB: Atom) extends DockingState(molA, molB) {

  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians
  override def getNeighbours = {
    val fwZ = molB.clone; fwZ.rotateZ(atomB.coords, DeltaAngle)   // forward rotation on Y axis
    val bwZ = molB.clone; bwZ.rotateZ(atomB.coords, -DeltaAngle)  // backward rotation on Y axis
    List(fwZ, bwZ).map(b => new AtomPairState(molA, atomA, b, atomB))
    // WARNING: atomB is NOT from the cloned molecule
  }
}
