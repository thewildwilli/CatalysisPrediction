package model

object VanDerWaalsRadii {
  val radii = Map(
    "H" -> 1.2,
    "C" -> 1.7,
    "N" -> 1.55,
    "O" -> 1.52,
    "F" -> 1.47,
    "P" -> 1.8,
    "Cl" -> 1.75
  )

  def apply(element: String) = radii.getOrElse(element, 1.4)

}
