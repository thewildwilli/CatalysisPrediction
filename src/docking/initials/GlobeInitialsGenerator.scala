package docking.initials
import model._
import breeze.linalg
import breeze.linalg.norm

class GlobeInitialsGenerator(initialConfigLevel: Integer,
                             angRad: Double) extends InitialsGenerator{

  override def apply(molA: Molecule, molB: Molecule) = {
    val radius = molA.getRadius + molB.getRadius;
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
          val newCentre = molB.getGeometricCentre + translateVector    // centre after translation
          val axis = linalg.cross(o, xone)
          val firstAngle = Math.asin(linalg.norm(axis))

          val transform = new MultiTransform(
            new Translate(translateVector),                             // translate out
            new Rotate(newCentre, axis, firstAngle),  // rotate to the orientation:
            new Rotate(newCentre, o, secondAngle)     // rotate on the axis of the orientation
          )

          // need to apply these transforms for the next step
          val bCopy = molB.clone
          transform.applyTo(bCopy)

          // get closest pair of atoms and translate B to get as close as possible to A
          val (a, b)  = (for (a <- molA.atoms; b <- bCopy.atoms) yield (a,b)).minBy(pair => pair._1.distTo(pair._2))
          val distance = Math.max(VanDerWaalsRadii(a.element), VanDerWaalsRadii(b.element))
          val bToA = a.coords - b.coords
          val approach = new Translate(bToA - bToA * distance / norm(bToA))

          result ::= new MultiTransform(transform.transforms :+ approach :_*)
        }
      }
    }
    result
  }
}
