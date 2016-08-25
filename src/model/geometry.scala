package model

import breeze.linalg
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
    *
    * @param axis: a 3-element [x, y, z] vector (additional elements are ignored)
    * @param angRad: rotation angle, in radians
    * @return a 4x4 transformation matrix
    */
  def rotate(axis: DenseVector[Double], angRad: Double) = {
    // Normalize the axis and add 1.0 at the end:
    val axisNorm = norm(axis)
    val normAxis = DenseVector[Double](axis(0)/axisNorm, axis(1)/axisNorm, axis(2)/axisNorm, 1.0)

    val crossprod = DenseMatrix(
      ( 0.0         , -normAxis(2)  ,  normAxis(1) ),
      ( normAxis(2) ,  0.0          , -normAxis(0) ),
      (-normAxis(1) ,  normAxis(0)  ,  0.0     )
    )

    //I3.cos(ang) + (1 - cos(ang)) . e . transpose(t) + [ex] . sin(ang) where [ex] is the cross product matrix
    val rotMatrix = DenseMatrix.eye[Double](3) * Math.cos(angRad) +
      (1 - Math.cos(angRad)) * normAxis * normAxis.t +
      crossprod * Math.sin(angRad)

    // Convert rotMatrix to a 4x4 matrix adding 0s and 1 in the diagonal:
    val res = DenseMatrix.horzcat(
      DenseMatrix.vertcat(rotMatrix, DenseMatrix((0.0, 0.0, 0.0))),
      DenseMatrix(0.0, 0.0, 0.0, 1.0)
    )

    if (axisNorm > 0)
      res
    else
      DenseMatrix.eye[Double](4)
  }

  def compose(transforms: DenseMatrix[Double]*) = {
    // multiply transforms in reverse order
    transforms.reduce((a, b) => b*a)
  }

  def transform(v: DenseVector[Double], t: DenseMatrix[Double]) =
    t * v

  /**
    * Generates a sequence of 3-dimensional vectors representing points around a sphere
    * with centre in the origin.
    * These points are generated by starting with a (radius, 0, 0) and rotating it
    * first with respect to the Y axis by angle each time, and then towards the poles.
    */
  def sphereOrientations(radius: Double, angle: Double) = {
    var northPoleDone = false
    var southPoleDone = false   // this is to avoid repeating the poles
    val northPole = DenseVector[Double](0.0, radius, 0.0)
    val southPole = DenseVector[Double](0.0, -radius, 0.0)

    val orientations = for (xzAngle <- 0.0 until Math.toRadians(180) by angle;
         otherAngle <- 0.0 until Math.toRadians(360) by angle) yield {

      // start with (radius, 0, 0) and rotate it counterclockwise in the XZ plane (a 2D operation):
      val rotateXZ = DenseVector[Double](radius * Math.cos(xzAngle), 0.0, radius * Math.sin(xzAngle), 1.0)

      // now get a perpendicular vector in the XZ plane (still a 2D operation):
      val axis = DenseVector[Double](-rotateXZ(2), 0.0, rotateXZ(0), 1.0);

      // now rotate with respect to that axis, and discard the 4th element
      val rotated = rotate(axis, otherAngle) * rotateXZ
      val rotated3d = rotated(0 to 2)

      val isNorthPole = norm(rotated3d-northPole) < radius*0.001
      val isSouthPole = norm(rotated3d-southPole) < radius*0.001

      val result =
        if (isNorthPole && northPoleDone) null
        else if (isSouthPole && southPoleDone) null
        else rotated3d

      if (isNorthPole) northPoleDone = true
      if (isSouthPole) southPoleDone = true

      result
    }

    orientations.filter(o => o != null)
  }

  /**
    * Project p onto the line that passes through x1 and x2, and return
    * the fraction t such that: projection = x1 + (x2-x1) t.
    * For t = 0, p projects on x1, for t = 1, p projects on x2.
    * http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    */
  def project(x1: DenseVector[Double], x2: DenseVector[Double], p: DenseVector[Double]) = {
    - ((x1 - p) dot (x2 - x1)) / Math.pow(norm(x2 - x1), 2)
  }

  /**
    * Distance of point p to the line that passes through x1 and x2.
    * http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    */
  def distToLine (x1: DenseVector[Double], x2: DenseVector[Double], p: DenseVector[Double]) = {
    norm(linalg.cross(p - x1, p - x2)) / norm(x2-x1)
  }



}
