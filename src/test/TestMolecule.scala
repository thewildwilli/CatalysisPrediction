package test

import breeze.linalg.DenseVector
import model.{Atom, Molecule}
import org.scalatest._

class TestMolecule extends FlatSpec{
  "A Molecule " should "calculate its geometric centre (1)" in {
    val a = new Molecule()
    a.Atoms.add(new Atom(-1, -1, 0))
    a.Atoms.add(new Atom(1, 1, 0))
    assert(a.getGeometricCentre == DenseVector(0.0, 0.0, 0.0))
  }

  it should "calculate its geometric centre (2)" in {
    val a = new Molecule()
    a.Atoms.add(new Atom(-3, 2, 5))
    a.Atoms.add(new Atom(3, 1, 2))
    a.Atoms.add(new Atom(4, -5, 6))
    assert(a.getGeometricCentre == DenseVector(4.0/3, -2.0/3, 13.0/3))
  }

  it should "calculate its radius (1)" in {
    val a = new Molecule()
    a.Atoms.add(new Atom(-1, 0, 0))
    a.Atoms.add(new Atom(5, 0, 0))
    assert(a.getRadius == 3.0)
  }

}
