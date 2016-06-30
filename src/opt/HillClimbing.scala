package opt

import io.threadcso._

// Created by Ernesto on 23/05/2016.
object HillClimbing{

  def optimize[S](init: S, actions: (S) => Seq[Action],
                           transition: (S, Action) => S,
                           scoring: S => Double,
                           maxIters: Int,
                           actionsOut: ![Action] = null): S = {
    var bestState = init
    var bestScore = Double.NegativeInfinity

    for (i <- 0 to maxIters) {
      var bestAction = null.asInstanceOf[Action]
      var bestNeighbour = null.asInstanceOf[S]
      var bestNeighbourScore = Double.NegativeInfinity
      for (a <- actions(bestState)) {
        val n = transition(bestState, a)
        val score = scoring(n)
        if (score > bestNeighbourScore) {
          bestAction = a
          bestNeighbourScore = score
          bestNeighbour = n
        }
      }

      if (bestNeighbourScore < bestScore) // reached local maximum and now going down
          return bestState

      bestState = bestNeighbour             // move to best neighbour
      bestScore = bestNeighbourScore
      if (actionsOut != null) actionsOut!bestAction // log actions chosen
    }
    bestState
  }
}
