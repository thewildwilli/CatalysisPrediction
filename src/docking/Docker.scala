package docking

import docking.dockscore.Scorer
import io.threadcso._
import model.Molecule
import opt.Action

// Created by Ernesto on 27/05/2016.
trait Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer,
           log: ![Any]): DockingState
}
