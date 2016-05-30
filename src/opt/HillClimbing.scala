package opt

// Created by Ernesto on 23/05/2016.
object HillClimbing{
  val maxDown = 10; // How many times are we allowed to go down before stopping?

  def optimize(initState: State, maxIters: Int, scoring: State => Double): State = {
    var currState = initState
    var currScore = Double.NegativeInfinity

    var bestState = initState
    var bestScore = Double.NegativeInfinity

    var downCount = 0

    for (i <- 0 to maxIters) {
      var bestNeighbour = null.asInstanceOf[State]
      var bestNeighbourScore = Double.NegativeInfinity
      for (n <- currState.getNeighbours) {
        val score = scoring(n);
        if (score > bestNeighbourScore) {
          bestNeighbourScore = score
          bestNeighbour = n
        }
      }

      if (bestNeighbourScore < currScore) { // reached local maximum and now going down
        if (downCount > maxDown)
          return bestState
        downCount += 1
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
