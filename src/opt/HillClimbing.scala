package opt

// Created by Ernesto on 23/05/2016.
object HillClimbing {

  def optimize(initState: State, maxIters: Int, scoring: State => Double): State = {
    var currState = initState
    var maxScore = Double.NegativeInfinity

    for (i <- 0 to maxIters){
      var bestNeighbour = null.asInstanceOf[State]
      var bestNeighbourScore = Double.NegativeInfinity
      for (n <- currState.getNeighbours){
        val score = scoring(n);
        if (score > bestNeighbourScore){
          bestNeighbourScore = score
          bestNeighbour = n
        }
      }

      if (bestNeighbourScore < maxScore) // reached local maxima and now going down
        return currState

      currState = bestNeighbour
      maxScore = bestNeighbourScore
    }
    currState
  }
}
