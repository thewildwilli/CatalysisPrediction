package docking.docksearch

import docking.DockingState
import model.Molecule

// Created by Ernesto on 23/05/2016.
object AtomPairDocker {
  def dock(a: Molecule, b: Molecule): DockingState = {
    for (atomA <- a.Atoms; atomB <- b.Atoms) {

      //    translate b so that the pair overlaps
      //    call hill climbing - neigbouring states are rotations
    }
    null
  }


}

/*class AtomPairState(a, b) extends DockingState(a, b) {
  override def getNeighbours = null
}*/
