package docking
import opt.{Action, State}
import model.{Molecule, Transform}

// Created by Ernesto on 23/05/2016.
class DockingState(val a: Molecule, val b: Molecule) extends State {

}

object DockingState {
  def transition(state: DockingState, action: Action) = {
    val bcopy = state.b.clone
    action match {
      case t: Transform => t.applyTo(bcopy)
    }
    new DockingState(state.a, bcopy)
  }
}

class TransformAction(val transform: Transform) extends Action
