package docking.docksearch

import breeze.linalg._
import docking.docksearch.initials.GlobeInitialsGenerator
import docking.{Docker}
import io.threadcso._
import model.{Molecule, Rotate, Translate, Transform}
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

  val initials = new GlobeInitialsGenerator(initialConfigLevel, angRad)

  type DockTask = (Molecule, Molecule, Transform, ![Any])
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
    for (transform <- initials(molB, radius)) {
      tasksChan!(molA, molB, transform, log)
    }
    tasksChan.closeOut
  }

  private def doTasks(tasksChan: ?[DockTask], resultsChan: ![DockResult]) = proc {
    repeat {
      val (molA, molB, transform, log) = tasksChan?;

      log ! "reset"
      log ! transform

      val bCopy = molB.clone
      transform.applyTo(bCopy)

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

  private def doRandomRotation(m: Molecule, log: ![Any]): Unit = {
    val axis = DenseVector(Random.nextDouble(), Random.nextDouble(), Random.nextDouble())
    val angle = Random.nextDouble() * 2 * Math.PI
    m.rotate(m.getGeometricCentre, axis, angle)
    log! new Rotate(m.getGeometricCentre, axis, angle)
  }
}
