package docking

import docking.dockscore.Scorer
import io.threadcso._
import model.Molecule
import opt.Action

trait Docker {
  /**
    * Returns the docked molecule and the score
    * @param molA
    * @param molB
    * @param log
    * @return
    */
  def dock(molA: Molecule, molB: Molecule, log: ![Any]): (Molecule, Double)
}
