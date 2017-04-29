package docking.initials
import model._
import breeze.linalg

class GlobeInitialsGenerator(initialConfigLevel: Integer,
                             angRad: Double) extends InitialsGenerator{

  override def apply(m: Molecule, radius: Double) = {
    val xone = linalg.DenseVector(1.0, 0.0, 0.0)
    var result = List[Transform]()
    val orientations = Geometry.sphereOrientations(1, angRad)
    for (pos <- orientations) {
      var firstOrientation = true
      for (o <- orientations if initialConfigLevel >= 1 || firstOrientation) {
        firstOrientation = false
        var first3dRotation = true
        for (secondAngle <- 0.0 until Math.toRadians(360) by angRad
             if initialConfigLevel >= 2 || first3dRotation ) {
          first3dRotation = false

          val translateVector = pos * radius
          val newCentre = m.getGeometricCentre + translateVector    // centre after translation
          val axis = linalg.cross(o, xone)
          val firstAngle = Math.asin(linalg.norm(axis))
          result ::= new MultiTransform(
            new Translate(translateVector),                             // translate out
            new Rotate(newCentre, axis, firstAngle),  // rotate to the orientation:
            new Rotate(newCentre, o, secondAngle)     // rotate on the axis of the orientation
          )
        }
      }
    }
    result
  }
}
