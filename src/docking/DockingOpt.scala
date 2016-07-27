package docking
import breeze.linalg.DenseVector
import opt.{State, Action}
import model.Molecule

// Created by Ernesto on 23/05/2016.
class DockingState(val a: Molecule, val b: Molecule) extends State {

}

object DockingState {
  def transition(state: DockingState, action: Action) = {
    val bcopy = state.b.clone
    action match {
      case t: Translate => bcopy.translate(t.v)
      case r: Rotate => bcopy.rotate(r.c, r.axis, r.angRad)
    }
    new DockingState(state.a, bcopy)
  }
}

/**
  * This action defines a translation with vector v
  * @param v: 3-dimension translation vector
  */
class Translate(val v: DenseVector[Double]) extends Action

/**
  * This action defines a rotation around centre c with respect
  * to axis, by angRad
  * @param c: centre, 3 dimensions
  * @param axis: axis vector, 3 dimensions
  * @param angRad: angle
  */
class Rotate(val c: DenseVector[Double], val axis: DenseVector[Double],
             val angRad: Double) extends Action {
  if (c.length != 3)
    throw new IllegalArgumentException("Centre must be a 3-dimensional vector")
  if (axis.length != 3)
    throw new IllegalArgumentException("Axis must be a 3-dimensional vector")
}