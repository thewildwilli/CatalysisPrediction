package docking.dockscore

import opt.State

trait Scorer {
  def score(state: State): Double
}
