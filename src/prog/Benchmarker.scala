package prog

import profiling.Profiler
import collection.mutable.Map

import prog.DockMain.DockArgs

// Created by Ernesto on 25/08/2016.
object Benchmarker {

  var repeats = 10
  var path = "benchmarkcmds.txt"
  var printClosestRefCount = false

  def main(args: Array[String]): Unit = {
    parseBenchmarkerArgs(args)


    val cmds = scala.io.Source.fromFile(path).getLines.toSeq
    var i = 1
    for (line <- cmds ) {
      val cmd = line.trim
      if (cmd.isEmpty)
        println()
      else if (!cmd.startsWith("#")) {
        val dockArgs = DockMain.parseArgs(cmd.split(" "))
        print(s"Line $i: ")
        Profiler.clear
        var rmsdAvg = 0.0

        // Count how many times each reference was the closest:
        val closestRefMap = Map[String, Int]().withDefaultValue(0)

        for (i <- 0 until repeats) {
          val (_, (closestRef, rmsd), _) = DockMain.doMainDock(dockArgs)
          print(".")
          rmsdAvg += rmsd
          closestRefMap(closestRef) += 1
        }
        rmsdAvg = rmsdAvg / repeats
        val avgTime = Profiler.getTimes("dock").toDouble / repeats
        println(s"Average time: ${avgTime}ms, Average RMSD: $rmsdAvg. Line: $cmd")
        if (printClosestRefCount){
          println("Times result was closest to each reference:")
          for ( (refPath, count) <- closestRefMap)
            println(s"$refPath: $count")
        }
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
        case "-closestRefCount" => printClosestRefCount = true; i += 1
        case _ => sys.error("wrong args")
      }
    }
  }

}
