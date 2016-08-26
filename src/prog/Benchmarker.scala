package prog

import profiling.Profiler
import prog.DockMain.DockArgs

// Created by Ernesto on 25/08/2016.
object Benchmarker {

  var repeats = 10
  var path = "benchmarkcmds.txt"

  def main(args: Array[String]): Unit = {
    parseBenchmarkerArgs(args)


    val cmds = scala.io.Source.fromFile(path).getLines.toSeq
    var i = 1
    for (line <- cmds ) {
      val cmd = line.trim
      if (!cmd.startsWith("#") && !cmd.isEmpty) {
        DockMain.parseArgs(cmd.split(" "))
        Profiler.clear
        var rmsdAvg = 0.0
        for (i <- 0 until repeats) {
          val (_, rmsd, _) = DockMain.doMainDock(DockArgs.fullPathA, DockArgs.fullPathB, DockArgs.fullPathOut)
          rmsdAvg += rmsd
        }
        rmsdAvg = rmsdAvg / repeats
        val avgTime = Profiler.getTimes("dock").toDouble / repeats
        println(s"Line $i: Average RMSD: $rmsdAvg, Average time: ${avgTime}ms. Line: $cmd")
      }
      i += 1
    }
    sys.exit(1)
  }

  def parseBenchmarkerArgs(args: Array[String]) = {
    var i = 0

    while (i < args.length) {
      args(i) match {
        case "-r" => repeats = args(i + 1).toInt;i += 2
        case "-path" => path = args(i + 1);i += 2
        case _ => sys.error("wrong args")
      }
    }
  }

}
