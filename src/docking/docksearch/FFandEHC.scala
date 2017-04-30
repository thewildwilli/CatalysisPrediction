package docking.docksearch

import docking.{DockLog, Docker}
import docking.dockscore.Scorer
import docking.docksearch.forcevector.{DockingParams, ForceVectorDocker}
import model.Molecule

// Created by Ernesto on 27/08/2016.
class FFandEHC(val ffparams: DockingParams,
               val ehcScorer: Scorer,
               val maxEhcIters: Integer) extends Docker {

  val ff = new ForceVectorDocker(ffparams)
  val ehc = new EhcDocker(ehcScorer, maxEhcIters)

  override def dock(molA: Molecule, molB: Molecule, log: DockLog) = {
    val (docked, _) = ff.dock(molA, molB, log)
    ehc.dock(molA, docked, log)
  }

}
