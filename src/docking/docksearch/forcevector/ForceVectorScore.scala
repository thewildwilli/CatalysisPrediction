package docking.docksearch.forcevector

import docking.dockscore.Scorer
import opt.State
import docking.DockingState

class ForceVectorScore(val params: DockingParams) extends Scorer {
  val docker = new ForceVectorDocker(params)

  override def score(state: State): Double = {
    val s = state.asInstanceOf[DockingState]
    docker.score(s.a, s.b)
  }
}
