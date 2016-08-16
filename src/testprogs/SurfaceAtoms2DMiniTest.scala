package testprogs

import io.{MoleculeReader, Pdb2DReader, XyzWriter}
import model.Molecule

object SurfaceAtoms2DMiniTest {
  def main(args: Array[String]): Unit ={
    val reader: MoleculeReader = new Pdb2DReader(args(0))
    val m: Molecule = reader.read
    m.computeSurfaceAtoms2D()
    for (a <- m.atoms)
      if (a.isSurface)
        a.setElement("N")
    new XyzWriter(args(1)).write(m)
  }
}
