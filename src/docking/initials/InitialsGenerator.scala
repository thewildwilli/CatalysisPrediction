package docking.initials

import breeze.linalg.norm
import model.{Molecule, Transform, Translate, VanDerWaalsRadii}

trait InitialsGenerator {
  def apply(molA: Molecule, molB: Molecule): Seq[Transform]

  /** Returns a translation that will bring molB close to molA
    * without overlapping, assuming they initially were not overlapping)
    */
  def approach(molA: Molecule, molB: Molecule) = {
    // get closest pair of atoms and translate B to get as close as possible to A
    val (a, b)  = (for (a <- molA.atoms; b <- molB.atoms) yield (a,b)).minBy(pair => pair._1.distTo(pair._2))
    val distance = Math.max(VanDerWaalsRadii(a.element), VanDerWaalsRadii(b.element))
    val bToA = a.coords - b.coords
    new Translate(bToA - bToA * distance / norm(bToA))
  }
}
