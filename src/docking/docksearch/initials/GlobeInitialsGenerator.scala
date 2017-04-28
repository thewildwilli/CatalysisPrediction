package docking.docksearch.initials
import io.threadcso.!
import model.{Geometry, Molecule}

import model.{Rotate, Translate}
import breeze.linalg

class GlobeInitialsGenerator(initialConfigLevel: Integer,
                             angRad: Double) extends InitialsGenerator{

  override def apply[A](m: Molecule, radius: Double, log: ![Any], cmd: Molecule => A) = {
    val xone = linalg.DenseVector(1.0, 0.0, 0.0)
    var l = List[A]()
    val orientations = Geometry.sphereOrientations(1, angRad)
    for (pos <- orientations) {
      var firstOrientation = true
      for (o <- orientations if initialConfigLevel >= 1 || firstOrientation) {
        firstOrientation = false
        var first3dRotation = true
        for (secondAngle <- 0.0 until Math.toRadians(360) by angRad
             if initialConfigLevel >= 2 || first3dRotation ) {
          first3dRotation = false

          log ! "reset"
          val bCopy = m.clone
          bCopy.translate(pos * radius)
          log ! new Translate(pos * radius)

          // rotate to the orientation:
          val axis = linalg.cross(o, xone)
          val firstAngle = Math.asin(linalg.norm(axis))
          bCopy.rotate(bCopy.getGeometricCentre, axis, firstAngle)
          log ! new Rotate(bCopy.getGeometricCentre, axis, firstAngle)

          // rotate on the axis of the orientation
          bCopy.rotate(bCopy.getGeometricCentre, o, secondAngle)
          log ! new Rotate(bCopy.getGeometricCentre, o, secondAngle)

          l ::= cmd(bCopy)
        }
      }
    }
    l
  }
}
