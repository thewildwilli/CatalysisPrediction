package testprogs

import jmolint.{JmolFrame, JmolMoleculeReader}
import model.Molecule

// Created by Ernesto on 26/08/2016.
object MeasureRMSD {

  val frame = new JmolFrame(500, 500, false)
  val jmolPanel = frame.getPanel

  /**
    * The first argument is the total number of models
    * Loads the files in the args, then computs RMSD of first model against all others, and finally reports min, max and average.
    */
  def main(args: Array[String]): Unit = {
    val modelCount = args(0).toInt
    jmolPanel.openFiles(args.slice(1, args.length), false)
    val refMol: Molecule = JmolMoleculeReader.read(jmolPanel, 0)

    var rmsds = List[Double]()

    for (i <- 1 until modelCount) {
      val other = JmolMoleculeReader.read(jmolPanel, i)
      try {
        val r = refMol.rmsd(other)
        rmsds ::= r
        println(s"RMSD model 0 to model $i: $r")
      } catch {
        case e:Exception => println(s"Could not measure RMSD to $i: ${e.getMessage}")
      }
    }

    println(s"Average RMSD = ${rmsds.sum / rmsds.size}")
    println(s"Max RMSD = ${rmsds.max}")
    println(s"Min RMSD = ${rmsds.min}")

    sys.exit(1)
  }

}
