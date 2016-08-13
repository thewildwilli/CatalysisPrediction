package jmolint

import model._
import scala.collection.mutable.ArrayBuffer


// Created by Ernesto on 13/07/2016.
object JmolMoleculeReader {
  def read(panel: JmolPanel, modelIndex: Int): Molecule = {
    val ms =  panel.viewer.ms
    var atoms = Map[Int, Atom]()  // Atoms by id

    // Create atoms
    for (jmolAtom <- ms.at)
      if (jmolAtom.mi == modelIndex)
        atoms += (jmolAtom.i -> new Atom(
          jmolAtom.i,
          jmolAtom.getElementSymbol,
          jmolAtom.x, jmolAtom.y, jmolAtom.z,
          jmolAtom.getPartialCharge,

          jmolAtom.getAtomName,
          jmolAtom.group.getSeqcodeString,
          jmolAtom.group.getGroup3
        ))

    // Run again through the list and create bonds
    for (jmolAtom <- ms.at)
      if (jmolAtom.mi == modelIndex)
        for (bond <- jmolAtom.bonds) {
          val other = bond.getOtherAtom(jmolAtom)
          atoms(jmolAtom.i).addBond(other.i)
        }

    new Molecule(atoms.values)
  }

  /*def read(m: Model) = {
    m.
    val atomList = ArrayBuffer[Atom]()
    for (i <- 0 until atomCol.getAtomCount()){
      val point = atomCol.getAtomPoint3f(i)
      atomList.append(new Atom(atomCol.getElementSymbol(i), point.x, point.y, point.z))
    }
    new Molecule(atomList)
  }*/
}
