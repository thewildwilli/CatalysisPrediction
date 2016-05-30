package docking.docksearch

import docking.{Docker, DockingState}
import docking.dockscore.Scorer
import model.{Atom, Molecule}
import opt.HillClimbing

// Created by Ernesto on 23/05/2016.
object AtomPairDocker extends Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState = {
    var maxScore = Double.NegativeInfinity
    var bestMatch = null.asInstanceOf[AtomPairState]
    var i=0;

    for (x <- 0 until molA.Atoms.size; y <- 0 until molB.Atoms.size) {
      val optimized = dockPair(molA, x, molB, y, scorer)
      val score = scorer.score(optimized)
      if (score > maxScore) {
        maxScore = score
        bestMatch = optimized
      }

      printf("docked %d of %d pairs, best score: %.4f %n", {i+=1;i}, molA.Atoms.size*molB.Atoms.size, maxScore);
    }
    bestMatch
  }

  private def dockPair(molA: Molecule, x: Int, molB: Molecule, y: Int, scorer: Scorer) = {
    //    translate molB so that the pair (atomA, atomB) overlaps
    molB.translate(molA(x).coords - molB(y).coords)

    //    call hill climbing - neigbouring states are rotations
    val initialState = new AtomPairState(molA, x, molB, y)
    HillClimbing.optimize(initialState, 50, scorer.score).asInstanceOf[AtomPairState]
  }
}

class AtomPairState(molA: Molecule, x: Int, molB: Molecule, y: Int) extends DockingState(molA, molB) {

  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians

  /** Neighbours are all rotations of molecule bMol around atom bAtom
    * by DeltaAngle in both directions in all 3 axis. That is, each
    * state has 6 neighbours. Molecule a is fixed.
    */
  override def getNeighbours = {
    val fwX = molB.clone; fwX.rotateX(fwX(y).coords, DeltaAngle)   // forward rotation on X axis
    val bwX = molB.clone; bwX.rotateX(bwX(y).coords, -DeltaAngle)  // backward rotation on X axis
    val fwY = molB.clone; fwY.rotateY(fwY(y).coords, DeltaAngle)   // forward rotation on Y axis
    val bwY = molB.clone; bwY.rotateY(bwY(y).coords, -DeltaAngle)  // backward rotation on Y axis
    val fwZ = molB.clone; fwZ.rotateZ(fwZ(y).coords, DeltaAngle)   // forward rotation on Z axis
    val bwZ = molB.clone; bwZ.rotateZ(bwZ(y).coords, -DeltaAngle)  // backward rotation on Z axis
    List(fwX, bwX, fwY, bwY, fwZ, bwZ).map(b => new AtomPairState(molA, x, b, y))
  }
}


class AtomPairState2D(molA: Molecule, x: Int, molB: Molecule, y: Int) extends DockingState(molA, molB) {

  final val DeltaAngle = Math.toRadians(20) // 20 degrees in radians
  override def getNeighbours = {
    val fwZ = molB.clone; fwZ.rotateZ(fwZ(y).coords, DeltaAngle)   // forward rotation on Y axis
    val bwZ = molB.clone; bwZ.rotateZ(bwZ(y).coords, -DeltaAngle)  // backward rotation on Y axis
    List(fwZ, bwZ).map(b => new AtomPairState(molA, x, b, y))
  }
}
