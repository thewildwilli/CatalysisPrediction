package test

/**
  * Run this program from intellij Idea and you'll get very weird results
  */
object TestOrientations {
  def main(args: Array[String]): Unit = {
    var i = 0
    for (xzAngle <- 0.0 until Math.toRadians(180) by Math.toRadians(90);
         otherAngle <- 0.0 until Math.toRadians(360) by Math.toRadians(90)) {

      i += 1
      println(s"hi $i ${Math.toDegrees(xzAngle)} ${Math.toDegrees(otherAngle)}")
    }
  }
}
