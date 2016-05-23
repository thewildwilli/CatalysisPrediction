// Created by Ernesto on 23/05/2016.

package model

import breeze.linalg.DenseVector

import collection.JavaConverters._

class Molecule(val Atoms: scala.collection.mutable.Set[Atom]) {

  def this(){ this (scala.collection.mutable.Set[Atom]()) }

  /** Deep copy */
  override def clone = {
    val result = new Molecule()
    for (a <- this.Atoms)
      result.Atoms.add(a.clone())
    result
  }

  def translate(v: DenseVector[Double]): Unit = {
    for (a <- this.Atoms)
      a.translate(v)
  }

  def JAtoms = Atoms.asJava


}
