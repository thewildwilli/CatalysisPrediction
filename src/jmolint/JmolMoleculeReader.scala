package jmolint

import model._
import scala.collection.mutable.ArrayBuffer


// Created by Ernesto on 13/07/2016.
object JmolMoleculeReader {
  def read(panel: JmolPanel, modelIndex: Int): Molecule = {
    val ms =  panel.viewer.ms
    val atoms = new ArrayBuffer[Atom]()

    for (jmolAtom <- ms.at)
      if (jmolAtom.mi == modelIndex)
        atoms.append(new Atom(jmolAtom.getElementSymbol, jmolAtom.x, jmolAtom.y, jmolAtom.z))
    new Molecule(atoms)
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
