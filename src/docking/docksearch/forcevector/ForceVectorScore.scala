package docking.docksearch.forcevector

import breeze.linalg.{DenseVector, norm}
import docking.dockscore.Scorer
import docking.docksearch.forcevector.Functions.explog2
import model.{BondEnergy, Molecule}
import opt.State
import HBond._
import Functions._
import docking.DockingState

class ForceVectorScore(val params: DockingParams,
                       val minCoverage: Double = 1.25,
                       val scoreOnlyTargetRadius: Boolean = false) extends Scorer {

  var avgBondEnergy = Double.NaN
  /**
    * Molecules assumed to be nonempty!
    *
    * @param molA
    * @param molB
    * @param onlyTargetRadius
    * @return
    */
  def getScore(molA: Molecule, molB: Molecule, onlyTargetRadius: Boolean) = {
    if (avgBondEnergy.isNaN)
      avgBondEnergy = getAvgBondEnergy(molA, molB)

    def hScore(vectToSpot: DenseVector[Double]) = {
      if (vectToSpot != null)
        gaussian(norm(vectToSpot)) //* a.partialCharge * b.partialCharge * -1.0
      else
        0.0
    }

    var totalGeometricScore = 0.0
    var totalElectricScore = 0.0
    var totalHBondScore = 0.0
    var totalBondStrengthScore = 0.0

    var aCount = 0
    var bCount = 0
    for (a <- molA.surfaceAtoms(params.ignoreAHydrogens) ){
      aCount += 1
      for(b <- molB.surfaceAtoms) {
        bCount += 1
        val actualDistance = a.distTo(b)
        val opt = optimalDistance(a, b, params.surface)
        if (!onlyTargetRadius || actualDistance <= opt * minCoverage) {
          val n = actualDistance / opt

          if (params.geometricForceWeight > 0)
            if (!a.isElement("H") && !b.isElement("H"))
              totalGeometricScore += explog2(n * explog2.maxX) / explog2.maxY

          if (params.electricForceWeight > 0) {
            val chargeProduct = a.partialCharge * b.partialCharge
            totalElectricScore += (
              if (chargeProduct < 0)
                (explog2(n * explog2.maxX) / explog2.maxY) * a.partialCharge * b.partialCharge * -1.0
              else if (chargeProduct > 0 && n <= 2)
                minusExpOverX(n) * chargeProduct
              else
                0.0
              )
          }

          if (params.hydrogenBondsForceWeight > 0)
            totalHBondScore += hScore(getVectToHBondSpot(a, b, molA))// + hScore(getVectToHBondSpot(b, a, molB))

          if (params.bondForceWeight > 0) {
            val bondEnergy = BondEnergy(a.element, b.element)
            val frac = bondEnergy / avgBondEnergy
            totalBondStrengthScore += (
              if (frac > 1)
                (explog2(n * explog2.maxX) / explog2.maxY) * frac
              else if (frac < 1 && n <= 2)
                minusExpOverX(n) * (1 - frac)
              else
                0.0
              )
          }
        }
      }
    }
    val totalScore =
      totalGeometricScore * params.geometricForceWeight +
        totalElectricScore * params.electricForceWeight +
        totalHBondScore * params.hydrogenBondsForceWeight +
        totalBondStrengthScore * params.bondForceWeight
    //println(s"SCORES: geo: ${totalGeometricScore * geometricForceWeight}, electric: ${totalElectricScore * electricForceWeight}, hbond: ${totalHBondScore * hydrogenBondsForceWeight}, bondstrength: ${totalBondStrengthScore * bondForceWeight}")
    totalScore / (Math.min(aCount, bCount))
  }

  def getAvgBondEnergy(molA: Molecule, molB: Molecule) =
    (for (a <- molA.atoms; b <- molB.atoms) yield BondEnergy(a.element, b.element)).sum /
      (molA.atoms.size * molB.atoms.size)

  override def score(state: State): Double = {
    val s = state.asInstanceOf[DockingState]
    getScore(s.a, s.b, scoreOnlyTargetRadius)
  }
}
