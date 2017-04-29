package docking.initials

import model.{Molecule, Transform}

trait InitialsGenerator {
  def apply(molA: Molecule, molB: Molecule): Seq[Transform]
}
