package prog

// Created by Ernesto on 26/05/2016.
import io._
import model.Molecule

object ConvertTo2D {
  def main(args: Array[String]) {
    var path: String = args(0)
    val reader: MoleculeReader = new Pdb2DReader(path)
    val m: Molecule = reader.read
    if (path.endsWith(".pdb") || path.endsWith(".PDB")) path = path.replace(".pdb", ".xyz")
    else path += ".xyz"
    new XyzWriter(path).write(m)
  }
}
