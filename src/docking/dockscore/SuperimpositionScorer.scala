package docking.dockscore

import docking.DockingState
import opt.State

// Created by Ernesto on 23/05/2016.
object SuperimpositionScorer extends Scorer {

  /** Returns 1 when all atoms in a are superimposed to some atom in b
    * Uses exp(-distance) to score each atom. O(n2) */
  def score(state: State) = {
    val s = state.asInstanceOf[DockingState]
    var score = 0.0
    for (atomA <- s.a.Atoms) {
      var bestA = 0.0
      val closestBAtom = s.b.Atoms.minBy(atomB => atomA.distTo(atomB))
      score += Math.exp(-atomA.distTo(closestBAtom))
    }
    score / (s.a.Atoms.size * s.b.Atoms.size)     // score as a fraction. 1 = perfect superposition.
  }
}
