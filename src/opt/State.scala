package opt

// Created by Ernesto on 23/05/2016.
trait State {
  def getNeighbours: List[State]
}
