package docking.docksearch

import docking.Docker
import docking.dockscore.Scorer
import docking.docksearch.forcevector.ForceVectorDocker
import io.threadcso._
import model.Molecule

// Created by Ernesto on 27/08/2016.
class FFandEHC(val surface: Double = 1.4,
               val permeability: Double = 0.5,
               val maxDecelerations: Int = 10,
               val ignoreAHydrogens: Boolean = false,
               val threshold: Double = 1.0e-5,
               val geometricForceWeight: Double = .25,
               val electricForceWeight: Double = .25,
               val hydrogenBondsForceWeight: Double = .25,
               val bondForceWeight: Double = .25,
               val ehcScorer: Scorer,
               val maxEhcIters: Integer) extends Docker {

  val ff = new ForceVectorDocker(surface, permeability, maxDecelerations, ignoreAHydrogens,
    threshold, geometricForceWeight, electricForceWeight, hydrogenBondsForceWeight, bondForceWeight)
  val ehc = new EhcDocker(ehcScorer, maxEhcIters)

  override def dock(molA: Molecule, molB: Molecule, log: ![Any]) = {
    val (docked, _) = ff.dock(molA, molB, log)
    ehc.dock(molA, docked, log)
  }

}
