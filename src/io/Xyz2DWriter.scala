package io

// Created by Ernesto on 23/05/2016.
import model.Molecule
import java.io.FileWriter

class Xyz2DWriter(val path: String) extends MoleculeWriter {
  /** Write to file using xyz format. Element name is not supported so C is always used.
    * See: http://openbabel.org/wiki/XYZ_(format) */
  def write(m: Molecule) {
    val fw: FileWriter = new FileWriter(this.path)
    try {
      fw.write(m.atoms.size + "\n") // First line: Atom count
      fw.write(this.path + "\n") // Second line: title
      import scala.collection.JavaConversions._
      for (a <- m.JAtoms.values) fw.write("%s %.3f %.3f %.3f%n".format(a.element, +a.x, a.y, 0.0))
    }
    finally if (fw != null) fw.close()

  }
}
