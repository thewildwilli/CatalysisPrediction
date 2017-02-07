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
      result.geometricForceWeight = 1/7.0
      result.electricForceWeight = 1/7.0
      result.hydrogenBondsForceWeight = 5/7.0
      result.bondForceWeight = 0.0
      result.ignoreAHydrogens = false
    } else if (molA.atoms.size < 300) {
      // Small molecule
      result.threshold = 1.0e-2
      result.surface = 0
      result.permeability = 0
      result.geometricForceWeight = 1/101.0
      result.electricForceWeight = 0.0
      result.hydrogenBondsForceWeight = 100/101.0
      result.bondForceWeight = 0.0
      result.ignoreAHydrogens = false
    } else {
      // Large molecule
      result.threshold = 1.0e-5
      result.surface = 1.4
      result.permeability = 0.9
      result.geometricForceWeight = 1/2.0
      result.electricForceWeight = 0.0
      result.hydrogenBondsForceWeight = 1/2.0
      result.bondForceWeight = 0.0
      result.ignoreAHydrogens = true
    }
    result
  }
}
