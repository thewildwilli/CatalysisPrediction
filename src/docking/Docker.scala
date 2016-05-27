package docking

import docking.dockscore.Scorer
import model.Molecule

// Created by Ernesto on 27/05/2016.
trait Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer): DockingState
}
