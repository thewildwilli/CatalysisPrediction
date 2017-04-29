package docking.initials

import model.{Molecule, Transform}

trait InitialsGenerator {
  def apply(m: Molecule, radius: Double): Seq[Transform]
}
