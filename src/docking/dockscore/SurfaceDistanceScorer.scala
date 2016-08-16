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
class SurfaceDistanceScorer(val surface: Double = 1.4) extends Scorer {

  def score(state: State) = {
    val s = state.asInstanceOf[DockingState]
    var score = 0.0
    for (atomA <- s.a.atoms; atomB <- s.b.atoms) {
      score += scoreDist(atomA, atomB)
    }
    score / (s.a.atoms.size * s.b.atoms.size).toDouble
  }

  def scoreDist(atomA: Atom, atomB: Atom) = {
    val optimal = atomA.radius + atomB.radius + 2 * surface   // optimal distance
    val actual = atomA.distTo(atomB)

    // Max is reached at sqrt(2) so normalize input so that optimal distance -> sqrt(2)
    val dNormalized = actual * Math.sqrt(2.0) / optimal
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

  private def expsquare(x: Double) = 10* Math.exp(-Math.pow(x,2)) * (Math.pow(x,2)-1)

  override def toString: String = "Surface distance scorer "
}
