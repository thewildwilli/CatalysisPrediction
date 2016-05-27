package docking.dockscore

import docking.DockingState
import model.Atom
import opt.State

/** radius measured in Armstorngs
  * This scoring function will give best score when an atom is at a distance
  * of 2*radius from another.
  * The problem is that atoms that are in the same place will be scored
  * negative infinite.
  * */
class ConstantSurfaceScorer(val radius: Double) extends Scorer {

  val maxPairScore = scoreDist(2*radius)

  def score(state: State) = {
    val s = state.asInstanceOf[DockingState]
    var score = 0.0
    for (atomA <- s.a.Atoms; atomB <- s.b.Atoms) {
      score += scoreDist(atomA.distTo(atomB)) / maxPairScore      // normalize to 0..1
    }
    score / (s.a.Atoms.size * s.b.Atoms.size)     // score as a fraction. 1 = perfect superposition.
  }

  def scoreDist(d: Double) = {
    val f1 = Math.exp(-Math.pow(d-2*radius, 2))   // exp(-((d-2*r)^2))
    val f2 = Math.log(d+0.5)/d
    f1+f2
  }

  override def toString: String = "ConstantSurface scorer with radius " + radius
}
