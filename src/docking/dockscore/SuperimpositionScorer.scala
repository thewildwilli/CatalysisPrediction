package docking.dockscore

import docking.DockingState

// Created by Ernesto on 23/05/2016.
object SuperimpositionScorer {

  /** Returns 1 when all atoms in a are superimposed to some atom in b
    * Uses exp(-distance) to score each atom. O(n2) */
  def score(s: DockingState) = {
    val a = s.a; val b = s.b
    for (atomA <- s.a.Atoms) {
      var bestA = 0.0
      for (atomB <- s.b.Atoms){
        val score = Math.exp(-atomA.distTo(atomB))
        bestA = Math.max(bestA, score)
      }
    }
  }
}
