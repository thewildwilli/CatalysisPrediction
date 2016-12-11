package docking.docksearch.forcevector

/**
  * Created by Ernesto on 10-Dec-16.
  */
class DockingParams {
  var surface: Double = 1.4
  var permeability: Double = 0.5
  var maxDecelerations: Int = 10
  var ignoreAHydrogens: Boolean = false
  var threshold: Double = 1.0e-5
  var geometricForceWeight: Double = .25
  var electricForceWeight: Double = .25
  var hydrogenBondsForceWeight: Double = .25
  var bondForceWeight: Double = .25
}
