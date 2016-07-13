package test

import io.PdbReader
import jmolint.{JmolFrame, JmolMoleculeReader}
import org.scalatest._
import model.{Atom, Molecule}

// Created by Ernesto on 13/07/2016.
class JmolMoleculeReaderTest extends FlatSpec{
    "A JmolMoleculeReaderTest " should "load atom coordinates the same way a PdbReader does" in {
      val molAFromPdb = new PdbReader("Hexamer.pdb").read
      val molBFromPdb = new PdbReader("HexamerComplement.pdb").read

      val frame = new JmolFrame(500, 500, false)
      val jmolPanel = frame.getPanel
      jmolPanel.openFiles(List("Hexamer.pdb", "HexamerComplement.pdb"))


      val molAFromJmol = JmolMoleculeReader.read(jmolPanel, 0)
      val molBFromJmol = JmolMoleculeReader.read(jmolPanel, 1)


      assertResult(molAFromPdb.Atoms.size)(molAFromJmol.Atoms.size)
      for (i <- molAFromPdb.Atoms.indices){
        assert(floatAreSimilar(molAFromPdb.Atoms(i).x, molAFromJmol.Atoms(i).x))
        assert(floatAreSimilar(molAFromPdb.Atoms(i).y, molAFromJmol.Atoms(i).y))
        assert(floatAreSimilar(molAFromPdb.Atoms(i).z, molAFromJmol.Atoms(i).z))
      }

      assertResult(molBFromPdb.Atoms.size) (molBFromJmol.Atoms.size)
      for (i <- molBFromPdb.Atoms.indices){
        assert(floatAreSimilar(molBFromPdb.Atoms(i).x, molBFromJmol.Atoms(i).x))
        assert(floatAreSimilar(molBFromPdb.Atoms(i).y, molBFromJmol.Atoms(i).y))
        assert(floatAreSimilar(molBFromPdb.Atoms(i).z, molBFromJmol.Atoms(i).z))
      }
    }

  def floatAreSimilar(a: Double, b: Double) = Math.abs(a-b) < 1.0e-5
}
