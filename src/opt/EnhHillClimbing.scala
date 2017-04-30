package opt

import docking.{DockLog}

// Created by Ernesto on 23/05/2016.
object EnhHillClimbing{
  val maxDown = 10; // How many times are we allowed to go down before stopping?

  def optimize[S](init: S, actions: (S) => Seq[Action],
                  transition: (S, Action) => S,
                  decelerate: () => Unit,
                  scoring: S => Double,
                  maxIters: Int,
                  log: DockLog): S = {

    var currState = init
    var currScore = Double.NegativeInfinity

    var bestState = init
    var bestScore = Double.NegativeInfinity

    var downCount = 0

    for (i <- 0 to maxIters) {
      var bestAction = null.asInstanceOf[Action]
      var bestNeighbour = null.asInstanceOf[S]
      var bestNeighbourScore = Double.NegativeInfinity
      for (a <- actions(currState)) {
        val n = transition(currState, a)
        val score = scoring(n)
        if (score > bestNeighbourScore) {
          bestAction = a
          bestNeighbourScore = score
          bestNeighbour = n
        }
      }

      log.action(bestAction)
      log.other(bestNeighbour)

      if (bestNeighbourScore < currScore) { // reached local maximum and now going down
        if (downCount > maxDown)
          return bestState
        downCount += 1

        decelerate() // went one step downhill
      }

      currState = bestNeighbour             // move to best neighbour
      currScore = bestNeighbourScore

      if (currScore > bestScore) {
        bestScore = currScore
        bestState = currState
      }
    }
    bestState
  }
}
