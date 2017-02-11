package io

// Created by Ernesto on 23/05/2016.
import model.Molecule
import java.io.File
import java.io.FileWriter

class Mol2Writer(val path: String) extends MoleculeWriter {
  /** Write to file using Mol2 format. It has the advantage of supporting charge.
    * See: http://www.tripos.com/data/support/mol2.pdf */
  def write(m: Molecule) {
      val fw: FileWriter = new FileWriter(this.path)
      try {
        fw.write("@<TRIPOS>MOLECULE\n")
        fw.write(new File(path).getName + "\n")
        fw.write(" %d %n".format(m.JAtoms.size))
        fw.write("SMALL\n")
        fw.write("CHARGES\n")
        var serial: Int = 1
        fw.write("@<TRIPOS>ATOM\n")
        import scala.collection.JavaConversions._
        for (a <- m.JAtoms.values) fw.write(" %5d %5s %18.8f %18.8f %18.8f %2s %4s %3s %18.8f%n".format({
          serial += 1;
          serial - 1
        }, a.atomName, a.x, a.y, a.z, a.element, if (a.substructureId != null && a.substructureId.length > 0) a.substructureId
        else "0", if (a.substructureName != null && a.substructureName.length > 0) a.substructureName
        else "A", a.partialCharge))
      }finally if (fw != null) fw.close()

  }
}
