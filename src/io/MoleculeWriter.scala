package io

import model.Molecule

/**
  * Created by Ernesto on 11-Feb-17.
  */
trait MoleculeWriter {
  def write(m: Molecule)
}
