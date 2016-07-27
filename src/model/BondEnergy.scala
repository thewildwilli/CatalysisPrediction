package model

import java.util.NoSuchElementException

/**
  *  This object holds a double entry table of bond dissociation energies. Given 2 elements,
  *  the table gives the energy needed to break a single bond of atoms of those two elements.
  *  Only a few elements supported so far.
  */
object BondEnergy {
  val energies = Map(
    "C"  -> Map("C" -> 607.000),
    "Cl" -> Map("C" -> 338.000, "Cl" -> 242.580),
    "H"  -> Map("C" -> 337.200, "Cl" -> 431.800, "F" -> 568.600),
    "F"  -> Map("C" -> 536.000, "Cl" -> 250.540, "F" -> 156.900, "H" -> 568.600),
    "N"  -> Map("C" -> 770.000, "Cl" -> 389.000, "F" -> 301.000, "H" -> 314.000, "N" -> 945.330),
    "O"  -> Map("C" -> 1076.50, "Cl" -> 272.000, "F" -> 222.000, "H" -> 428.000, "N" -> 630.570, "O" -> 498.340),
    "P"  -> Map("C" -> 513.000, "Cl" -> 289.000, "F" -> 439.000, "H" -> 343.000, "N" -> 617.000, "O" -> 596.000, "P" -> 490.000)
  )

  def apply(element1: String, element2: String) = {
    try {
      energies(element1)(element2)
    } catch {
      case e:NoSuchElementException =>
        try {
          energies(element2)(element1)
        }
        catch {
          case e: NoSuchElementException => 0.0
        }
    }
  }
}
