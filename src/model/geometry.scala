package model

import breeze.linalg._

/**
  *   Provides general 3D geometry routines.
  *   Transformations are expressed in 4x4 matrices.
  *   3D vectors must be in the format [x, y, z, 1]
  */
object Geometry {



  /** Returns transformation matrix */
  def translate(v: DenseVector[Double]) = DenseMatrix (
    (1.0, 0.0, 0.0, v(0)),
    (0.0, 1.0, 0.0, v(1)),
    (0.0, 0.0, 1.0, v(2)),
    (0.0, 0.0, 0.0, 1.0)
    )

  /**
    * https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis-angle
    * @param axis: a 4-element [x, y, z, 1] vector.
    * @param angRad: rotation angle, in radians
    * @return a 4x4 transformation matrix
    */
  def rotate(axis: DenseVector[Double], angRad: Double) = {
    val crossprod = DenseMatrix(
      ( 0.0     , -axis(2)  ,  axis(1) ),
      ( axis(2) ,  0.0      , -axis(0) ),
      (-axis(1) ,  axis(0)  ,  0.0     )
    )

    //I3.cos(ang) + (1 - cos(ang)) . e . transpose(t) + [ex] . sin(ang) where [ex] is the cross product matrix
    val rotMatrix = DenseMatrix.eye[Double](3) * Math.cos(angRad) +
      (1 - Math.cos(angRad)) * axis * axis.t +
      crossprod * Math.sin(angRad)

    // Convert rotMatrix to a 4x4 matrix adding 0s and 1 in the diagonal:
    DenseMatrix.horzcat(
      DenseMatrix.vertcat(rotMatrix, DenseMatrix((0.0, 0.0, 0.0))),
      DenseMatrix(0.0, 0.0, 0.0, 1.0)
    )
  }

  def compose(transforms: DenseMatrix[Double]*) = {
    // multiply transforms in reverse order
    transforms.reduce((a, b) => b*a)
  }

  def transform(v: DenseVector[Double], t: DenseMatrix[Double]) =
    t * v




}
