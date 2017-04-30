package docking.initials
import breeze.linalg.{DenseVector, norm}
import model._

import scala.util.Random

class RandomInitials(val n: Integer) extends InitialsGenerator {
  override def apply(molA: Molecule, molB: Molecule): Seq[Transform] =
    (1 to n).map(_ => getRandomInitial(molA, molB))

  def getRandomInitial(molA: Molecule, molB: Molecule) = {
    // we need a copy for the final steps
    val bCopy = molB.clone
    val radius = molA.getRadius + molB.getRadius

    // get a random orientation to translate to:
    val orientationDir = DenseVector(Random.nextDouble(), Random.nextDouble(), Random.nextDouble())
    val orientation = (orientationDir / norm(orientationDir)) * radius

    // move there:
    val translate = new Translate(orientation)
    translate.applyTo(bCopy)

    // do a random rotation:
    val (axis, angle) = Geometry.randomAxisAndAngle
    val rotate = new Rotate(bCopy.getGeometricCentre, axis, angle)
    rotate.applyTo(bCopy)

    // approach the target molecule:
    new MultiTransform(translate, rotate, approach(molA, bCopy))
  }
}
