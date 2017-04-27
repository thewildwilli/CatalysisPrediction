package docking.docksearch

import breeze.linalg
import breeze.linalg._
import docking.docksearch.initials.GlobeInitialsGenerator
import docking.{Docker, Rotate, Translate}
import io.threadcso.!
import model.{Geometry, Molecule}
import profiling.Profiler

import scala.util.Random

/**
  * @param initialConfigLevel: 0 => only translation, 1 => orientations in 2D, 2 => orientations in 3D.
  *                         For angle=90 degrees, 0 gives 6 initial configurations, 1 gives 36 and 2 gives 144.
 */
class MultipleInitialsDocker(val docker: Docker, angRad: Double,
                             initialConfigLevel: Integer,
                             randomInitial: Boolean) extends Docker {

  val xone = DenseVector(1.0, 0.0, 0.0)
  val initials = new GlobeInitialsGenerator(initialConfigLevel, angRad)

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]): (Molecule, Double) = {
    val centreVect = molA.getGeometricCentre - molB.getGeometricCentre
    log!new Translate(centreVect)
    molB.translate(centreVect)
    if (randomInitial)
    doRandomRotation(molB, log)

    log!"save"

    val radius = molA.getRadius + molB.getRadius
    val best = initials(molB, radius, log, bCopy => {
      Profiler.time("dock") {  docker.dock(molA, bCopy, log) }
    }).maxBy(p => p._2)
    best
  }

  private def doRandomRotation(m: Molecule, log: ![Any]): Unit = {
    val axis = DenseVector(Random.nextDouble(), Random.nextDouble(), Random.nextDouble())
    val angle = Random.nextDouble() * 2 * Math.PI
    m.rotate(m.getGeometricCentre, axis, angle)
    log! new Rotate(m.getGeometricCentre, axis, angle)
  }
}
