package docking.docksearch.forcevector

import model.{Atom, Molecule}

object HBond {
  val hBondDistance = 1.8
  def getVectToHBondSpot(h: Atom, other: Atom, hmol: Molecule) = {
    if (canHBond(h, other)){
      val hNeighbour = hmol(h.bonds(0));
      if (hNeighbour.partialCharge < 0 && hNeighbour.isOneOf("F", "O", "N")) {
        val (_, _, dirToNeighbour) = h.distDifDir(hNeighbour)
        val toSpot = -(dirToNeighbour) * hBondDistance // 1.8 Angtrom in the opposite direction
        val spot = h.coords + toSpot
        spot - other.coords
      } else
        null
    }
    else
      null
  }

  def canHBond(h: Atom, other: Atom) =
    h.isH && h.partialCharge > 0 && other.partialCharge < 0 &&
      h.bonds.size == 1 && other.isOneOf("F", "O", "N")

}
