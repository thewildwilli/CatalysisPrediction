package model
import breeze.linalg.DenseVector
import opt.Action

trait Transform extends Action {
  def applyTo(m: Molecule);
}

/**
  * This action defines a translation with vector v
  * @param v: 3-dimension translation vector
  */
class Translate(val v: DenseVector[Double]) extends Transform {
  override def applyTo(m: Molecule): Unit = m.translate(v)
}

/**
  * This action defines a rotation around centre c with respect
  * to axis, by angRad
  * @param c: centre, 3 dimensions
  * @param axis: axis vector, 3 dimensions
  * @param angRad: angle
  */
class Rotate(val c: DenseVector[Double], val axis: DenseVector[Double],
             val angRad: Double) extends Transform {
  if (c.length != 3)
    throw new IllegalArgumentException("Centre must be a 3-dimensional vector")
  if (axis.length != 3)
    throw new IllegalArgumentException("Axis must be a 3-dimensional vector")

  override def applyTo(m: Molecule): Unit = m.rotate(c, axis, angRad)
}

class MultiTransform(val transforms: Transform*) extends Transform {
  override def applyTo(m: Molecule): Unit = {
    for (t <- transforms) t.applyTo(m)
  }
}
