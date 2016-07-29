package docking

import docking.dockscore.Scorer
import io.threadcso._
import model.Molecule
import opt.Action

trait Docker {
  def dock(molA: Molecule, molB: Molecule, scorer: Scorer,
           log: ![Any]): DockingState
}
