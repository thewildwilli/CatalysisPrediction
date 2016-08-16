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


      assertResult(molAFromPdb.atoms.size)(molAFromJmol.atoms.size)
      for (i <- molAFromPdb.atomMap.keys){
        assert(floatAreSimilar(molAFromPdb(i).x, molAFromJmol(i).x))
        assert(floatAreSimilar(molAFromPdb(i).y, molAFromJmol(i).y))
        assert(floatAreSimilar(molAFromPdb(i).z, molAFromJmol(i).z))
      }

      assertResult(molBFromPdb.atoms.size) (molBFromJmol.atoms.size)
      for (i <- molBFromPdb.atomMap.keys){
        assert(floatAreSimilar(molBFromPdb(i).x, molBFromJmol(i).x))
        assert(floatAreSimilar(molBFromPdb(i).y, molBFromJmol(i).y))
        assert(floatAreSimilar(molBFromPdb(i).z, molBFromJmol(i).z))
      }
    }

  def floatAreSimilar(a: Double, b: Double) = Math.abs(a-b) < 1.0e-5
}
