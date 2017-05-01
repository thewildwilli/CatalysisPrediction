package jmolint

import breeze.linalg.DenseVector
import model._


// Created by Ernesto on 13/07/2016.
object JmolMoleculeReader {
  def read(panel: JmolPanel, modelIndex: Int, readSurface: Boolean = true): Molecule = {
    val ms =  panel.viewer.ms
    var atoms = Map[Int, Atom]()  // Atoms by id

    if (readSurface)
      ms.calculateSurface(null, -1)  // copied from AtomCollection.java line 409.

    // Create atoms
    for (jmolAtom <- ms.at)
      if (jmolAtom.mi == modelIndex) {
        atoms += (jmolAtom.i -> new Atom(
          jmolAtom.i,
          jmolAtom.getElementSymbol,
          jmolAtom.x, jmolAtom.y, jmolAtom.z,
          jmolAtom.getPartialCharge,

          jmolAtom.getAtomName,
          jmolAtom.group.getSeqcodeString,
          jmolAtom.group.getGroup3,
          List(),
          if (readSurface) jmolAtom.getSurfaceDistance100() == 0 else false
        ))
      }

    // Run again through the list and create bonds
    for (jmolAtom <- ms.at)
      if (jmolAtom.mi == modelIndex && jmolAtom.bonds != null)
        for (bond <- jmolAtom.bonds) {
          val other = bond.getOtherAtom(jmolAtom)
          atoms(jmolAtom.i).addBond(other.i)
        }

    new Molecule(atoms.values)
  }
}
