package docking.dockscore

import docking.DockingState
import model.Atom
import opt.State

/** radius measured in Angstroms.
  * This scoring function will give best score when an atom is at a distance
  * of 2*radius from another.
  * The problem is that atoms that are in the same place will be scored
  * negative infinite.
  * */
class SurfaceDistanceScorer extends Scorer {

  def score(state: State) = {
    val s = state.asInstanceOf[DockingState]
    var score = 0.0
    for (atomA <- s.a.Atoms; atomB <- s.b.Atoms) {
      score += scoreDist(atomA, atomB)
    }
    score / (s.a.Atoms.size * s.b.Atoms.size)     // score as a fraction. 1 = perfect superposition.
  }

  def scoreDist(atomA: Atom, atomB: Atom) = {
    val radii = atomA.radius + atomB._radius
    val d = atomA.distTo(atomB)
    val dNormalized = (d/radii)*1.41421356237310
    expsquare(dNormalized)
  }


  /* Older function:
  def gauss(x: Double) = Math.exp(-Math.pow(x, 2))
  def softstep(x: Double) = {
    val skew = 5.0
    val trans = -1.0
    val dragdown = 20
    val base = 1/(1+Math.exp(- ((x-trans)*skew) ))-0.5
    base*(dragdown+1) - (dragdown/2-1)
  }*/

  // Max is reached at sqrt(2) so normalize input so that radiusA+radiusB -> sqrt(2)
  def expsquare(x: Double) = 10*Math.exp(-Math.pow(x,2))*(Math.pow(x,2)-1)

  override def toString: String = "Surface distance scorer "
}
