package docking

import model.Molecule

trait Docker {
  /**
    * Returns the docked molecule and the score
    * @param molA
    * @param molB
    * @param log
    * @return
    */
  def dock(molA: Molecule, molB: Molecule, log: DockLog): (Molecule, Double)
}
