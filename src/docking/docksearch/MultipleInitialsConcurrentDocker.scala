package docking.docksearch

import breeze.linalg
import breeze.linalg._
import docking.{Docker, Rotate, Translate}
import io.threadcso._
import model.{Geometry, Molecule}
import profiling.Profiler

import scala.util.Random


/**
  * @param initialConfigLevel: 0 => only translation, 1 => orientations in 2D, 2 => orientations in 3D.
  *                         For angle=90 degrees, 0 gives 6 initial configurations, 1 gives 36 and 2 gives 144.
 */
class MultipleInitialsConcurrentDocker(val createDocker : () => Docker,
                                       angRad: Double,
                                       initialConfigLevel: Integer,
                                       randomInitial: Boolean,
                                       val workers: Int = 8) extends Docker {

  val xone = DenseVector(1.0, 0.0, 0.0)
  type DockTask = (Molecule, Molecule, ![Any])
  type DockResult = (Molecule, Double)

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]): (Molecule, Double) = {
    val centreVect = molA.getGeometricCentre - molB.getGeometricCentre
    log!new Translate(centreVect)
    molB.translate(centreVect)
    if (randomInitial)
      doRandomRotation(molB, log)

    log!"save"

    val tasksChan = N2N[DockTask](1, workers, "taskschan")
    val resultsChan = N2N[DockResult](workers, 1, "resultschan")
    var result = (null.asInstanceOf[Molecule], 0.0)
    val workerProcs = || (for (i <- 0 until workers) yield doTasks(tasksChan, resultsChan))
    val system = (pushTasks(tasksChan, molA, molB, log) || workerProcs ||
      proc { result = collectResults(resultsChan) })
    Profiler.time("dock") { system() }
    result
  }

  private def pushTasks(tasksChan: ![DockTask], molA: Molecule, molB: Molecule, log: ![Any]) = proc {
    val radius = molA.getRadius + molB.getRadius;
    forOrientations(molB, radius, log, bCopy => tasksChan!(molA, bCopy, log))
    tasksChan.closeOut
  }

  private def doTasks(tasksChan: ?[DockTask], resultsChan: ![DockResult]) = proc {
    repeat {
      val (molA, molB, log) = tasksChan?;
      val docker = createDocker()
      val result = docker.dock(molA, molB, log)
      resultsChan!result
    }
    tasksChan.closeIn
    resultsChan.closeOut
  }

  private def collectResults(resultsChan: ?[DockResult]) = {
    var bestDock = null.asInstanceOf[Molecule]
    var bestScore = Double.NegativeInfinity
    repeat {
      val (molB, score) = resultsChan?;
      if (score > bestScore) {
        bestScore = score;
        bestDock = molB
      }
    }
    resultsChan.closeIn
    (bestDock, bestScore)
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

  private def doRandomRotation(m: Molecule, log: ![Any]): Unit = {
    val axis = DenseVector(Random.nextDouble(), Random.nextDouble(), Random.nextDouble())
    val angle = Random.nextDouble() * 2 * Math.PI
    m.rotate(m.getGeometricCentre, axis, angle)
    log! new Rotate(m.getGeometricCentre, axis, angle)
  }
}
