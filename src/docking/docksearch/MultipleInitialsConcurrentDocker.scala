package docking.docksearch

import docking.initials.InitialsGenerator
import docking.{DockLog, Docker}
import io.threadcso._
import model._
import profiling.Profiler


class MultipleInitialsConcurrentDocker(val createDocker : () => Docker,
                                       val initials: InitialsGenerator,
                                       val workers: Int = 8) extends Docker {

  type DockTask = (Molecule, Molecule, Transform, DockLog)
  type DockResult = (Molecule, Double)

  override def dock(molA: Molecule, molB: Molecule, log: DockLog): (Molecule, Double) = {
    val translateCentre = new Translate(molA.getGeometricCentre - molB.getGeometricCentre)
    translateCentre.applyTo(molB)
    log.action(translateCentre)

    log.save

    val tasksChan = N2N[DockTask](1, workers, "taskschan")
    val resultsChan = N2N[DockResult](workers, 1, "resultschan")
    var result = (null.asInstanceOf[Molecule], 0.0)
    val workerProcs = || (for (i <- 0 until workers) yield doTasks(tasksChan, resultsChan))
    val system = pushTasks(tasksChan, molA, molB, log) || workerProcs ||
      proc { result = collectResults(resultsChan) }
    Profiler.time("dock") { system() }
    result
  }

  private def pushTasks(tasksChan: ![DockTask], molA: Molecule, molB: Molecule, log: DockLog) = proc {
    for (transform <- initials(molA, molB)) {
      tasksChan!(molA, molB, transform, log)
    }
    tasksChan.closeOut
  }

  private def doTasks(tasksChan: ?[DockTask], resultsChan: ![DockResult]) = proc {
    repeat {
      val (molA, molB, transform, log) = tasksChan?;

      log.reset
      log.action(transform)

      val bCopy = molB.clone
      transform.applyTo(bCopy)

      val docker = createDocker()
      val result = docker.dock(molA, bCopy, log)
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
}
