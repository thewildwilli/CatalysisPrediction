package docking.docksearch.forcevector

import model.Molecule

/**
  * Created by Ernesto on 10-Dec-16.
  */
object DockingParamsHeuristic {
  def estimate(molA: Molecule, molB: Molecule) = {
    val result = new DockingParams()
    if (molA.atoms.size < 20) {
      // Tiny molecule
      result.threshold = 1.0e-5
      result.surface = 0
      result.permeability = 0
      result.geometricForceWeight = 1/10.0
      result.electricForceWeight = 3/10.0
      result.hydrogenBondsForceWeight = 6/10.0
      result.bondForceWeight = 0/7.0
      result.ignoreAHydrogens = false
    } else if (molA.atoms.size < 300) {
      // Small molecule
      result.threshold = 1.0e-5
      result.surface = 0
      result.permeability = 0
      result.geometricForceWeight = 2.5/10.0
      result.electricForceWeight = 0.5/10.0
      result.hydrogenBondsForceWeight = 7/10.0
      result.bondForceWeight = 0.0
      result.ignoreAHydrogens = false
    } else if (molA.atoms.size < 2000) {
      // Large molecule
      result.threshold = 1.0e-5
      result.surface = 1.4
      result.permeability = 0.9
      result.geometricForceWeight = 5/10.0
      result.electricForceWeight = 0.0
      result.hydrogenBondsForceWeight = 5/10.0
      result.bondForceWeight = 0.0
      result.ignoreAHydrogens = true
    } else {
      // huge molecule
      result.threshold = 1.0e-2
      result.surface = 1.4
      result.permeability = 0.9
      result.geometricForceWeight = 100/100.0
      result.electricForceWeight = 0.0
      result.hydrogenBondsForceWeight = 0/100.0
      result.bondForceWeight = 0.0
      result.ignoreAHydrogens = true
    }
    result
  }
}
