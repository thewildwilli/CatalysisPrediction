package docking.docksearch

import breeze.linalg
import breeze.linalg._
import docking.{Docker, Rotate, Translate}
import io.threadcso.!
import model.{Geometry, Molecule}
import profiling.Profiler

/**
  * @param initialConfigLevel: 0 => only translation, 1 => orientations in 2D, 2 => orientations in 3D.
  *                         For angle=90 degrees, 0 gives 6 initial configurations, 1 gives 36 and 2 gives 144.
 */
class MultipleInitialsDocker(val docker: Docker, angRad: Double,
                             initialConfigLevel: Integer) extends Docker {

  val xone = DenseVector(1.0, 0.0, 0.0)

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]): (Molecule, Double) = {
    val centreVect = molA.getGeometricCentre - molB.getGeometricCentre
    log!new Translate(centreVect)
    molB.translate(centreVect)
    log!"save"

    val radius = molA.getRadius + molB.getRadius
    val best = forOrientations(molB, radius, log, bCopy =>
      Profiler.time("dock") { docker.dock(molA, bCopy, log)}).maxBy(p => p._2)
    Profiler.report
    best
  }

  private def forOrientations[A](m: Molecule, radius: Double, log: ![Any], cmd: Molecule => A): Seq[A] = {
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
          val firstAngle = Math.asin(norm(axis))
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

  private def dockFromPos(molA: Molecule, molB: Molecule, log: ![Any]) =
    Profiler.time("dock") {  docker.dock(molA, molB, log) }
    /*var l = List[A]()
    val orientations = Geometry.sphereOrientations(1, angRad)
    for (pos <- orientations) {
      for (o <- orientations) {
        for (secondAngle <- 0.0 until Math.toRadians(360) by angRad) {
          l ::= cmd
        }
      }
    }
    l*/
}
