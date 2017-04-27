package docking.docksearch.initials

import io.threadcso.!
import model.Molecule

trait InitialsGenerator {
  def apply[A](m: Molecule, radius: Double, log: ![Any], cmd: Molecule => A): Seq[A]
}
