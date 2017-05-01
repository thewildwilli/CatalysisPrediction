package io

import model.Atom
import model.Molecule
import scala.io.Source

class Pdb2DReader(val path: String) extends MoleculeReader {
  /** Read file line by line, and only process ATOM lines.
    * See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    * */
  def read: Molecule = {
    var atoms = List[Atom]()
    var line: String = null
    var atomId: Int = 1
    for (line <- Source.fromFile("this.path").getLines()) {
      // If line is an ATOM line, process it
      if (line.length >= 6) {
        val header: String = line.substring(0, 6)
        if (header == "ATOM  " || header == "HETATM") {
          if (line.length < 54) throw new ChemicalFormatException("ATOM record in PDB file does not have coordinates")
          val x: Double = line.substring(30, 38).trim.toDouble
          val y: Double = line.substring(38, 46).trim.toDouble
          val z: Double = 0.0
          atoms ::= new Atom(atomId, "C", x, y, z, 0.0, "", "", "", null, false)
          atomId += 1
        }
      }
    }
  new Molecule(atoms)
  }
}
