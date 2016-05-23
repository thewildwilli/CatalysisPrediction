package opt

// Created by Ernesto on 23/05/2016.
object HillClimbing {
  def optimize(initState: State, maxIters: Int, scoring: State => Double): State = {
    var currState = initState
    for (i <- 0 to maxIters){
      val nextState = currState.getNeighbours.maxBy(scoring)
      if (currState == nextState)
        return currState
      currState = nextState
    }
    currState
  }
}
