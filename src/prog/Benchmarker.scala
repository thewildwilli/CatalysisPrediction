package prog

import collection.mutable.Map
import java.io.{File, FileWriter}
import java.text.SimpleDateFormat
import java.util.Date

import profiling.Profiler

// Created by Ernesto on 25/08/2016.
object Benchmarker {

  var repeats = 10
  var commandsPath = "benchmarkcmds.txt"
  var printClosestRefCount = false

  var initials: Iterable[Int] = List[Int]()

  val csvPath = new SimpleDateFormat("YYYY.MM.dd_HH.mm.ss").format(new Date()) + ".csv"

  def main(args: Array[String]): Unit = {
    parseBenchmarkerArgs(args)
    val lines = scala.io.Source.fromFile(commandsPath).getLines.toSeq

    var i = 1
    for (line <- lines ) if (!line.startsWith("#")) {
      val (experiment, cmd) = getExperimentAndCmd(line)
      if (cmd.isEmpty)
        println()
      else {
        if (initials.size > 0) {
          for (initNumber <- initials) {
            val extendedCmd = cmd + s" -initials random $initNumber"
            runExperiment(experiment, i, extendedCmd)
          }
        }
        else runExperiment(experiment, i, cmd)
      }
      i += 1
    }
    sys.exit(0)
  }

  private def runExperiment(experiment: String, lineNumber: Int,  cmd: String): Unit = {
    val dockArgs = DockMain.parseArgs(cmd.split(" "))
    print(s"Line $lineNumber: ")
    Profiler.clear
    var rmsdAvg = 0.0

    // Count how many times each reference was the closest:
    val closestRefMap = Map[String, Int]().withDefaultValue(0)
    var rmsds = List[Double]()
    var scores = List[Double]()

    for (_ <- 0 until repeats) {
      val (_, (closestRef, rmsd), score) = DockMain.doMainDock(dockArgs)
      print(".")
      rmsds ::= rmsd
      scores ::= score
      closestRefMap(closestRef) += 1
    }
    val avgTime = Profiler.getTimes("dock").toDouble / repeats

    val initNumber = if (dockArgs.initials == "globe") dockArgs.initConfigLevel else dockArgs.initNumber
    reportToConosole(cmd, avgTime, rmsds, closestRefMap)
    reportToCsvFile(csvPath, experiment, dockArgs.initials, initNumber, avgTime, rmsds, closestRefMap, scores)
  }

  private def getExperimentAndCmd(line: String) = {
    val trimmed = line.trim()
    val pos = trimmed.indexOf(':')
    if (pos < 0)
      ("", trimmed)
    else
      (trimmed.substring(0, pos).trim(), trimmed.substring(pos+1).trim())
  }

  private def reportToConosole(cmd: String, avgTime: Double, rmsds: Seq[Double],
                               closestRefMap: Map[String, Int]) = {
    val rmsdAvg = rmsds.sum / repeats
    println(s"Average time: ${avgTime}ms, Average RMSD: $rmsdAvg. Line: $cmd")
    if (printClosestRefCount){
      println("Times result was closest to each reference:")
      for ( (refPath, count) <- closestRefMap)
        println(s"$refPath: $count")
    }
  }

  private def reportToCsvFile(csvPath: String, experiment: String,
                              initType: String, initNumber: Int,
                              avgTime: Double, rmsds: Seq[Double],
                              closestRefMap: Map[String, Int],
                              scores: Seq[Double]): Unit ={

    def writeCSVHeader(csvPath: String): Unit = {
      println(s"creating $csvPath")
      val fw: FileWriter = new FileWriter(csvPath, true)
      try {
        fw.write(String.join(",", "Experiment", "InitType", "InitNumber", "TimeAvg(ms)", "RMSDAvg",
          "RMSDStdDev", "RMSDMin", "RMSDMax", "ScoreAvg", "ScoreStdDev") + "\n")
      } finally {fw.close()}
    }

    if (! new File(csvPath).exists()) writeCSVHeader(csvPath)

    val rmsdAvg = rmsds.sum / rmsds.length
    val rmsdStdDev = stdDev(rmsds)
    val rmsdMin = rmsds.min
    val rmsdMax = rmsds.max
    val scoreAvg = scores.sum / scores.length
    val scoreStdDev = stdDev(scores)

    // write row
    val fw: FileWriter = new FileWriter(csvPath, true)
    try {
      fw.write(String.join(",", experiment, initType, initNumber.toString,
        avgTime.toString, rmsdAvg.toString, rmsdStdDev.toString,
        rmsdMin.toString, rmsdMax.toString, scoreAvg.toString, scoreStdDev.toString) + "\n")
    } finally {fw.close()}
  }

  def parseBenchmarkerArgs(args: Array[String]) = {
    var i = 0

    while (i < args.length) {
      args(i) match {
        case "-r" => repeats = args(i + 1).toInt;i += 2
        case "-path" => commandsPath = args(i + 1);i += 2
        case "-closestRefCount" => printClosestRefCount = true; i += 1
        case "-initials" =>
          initials = args(i+1).split(",").map(s => s.toInt);
          i += 2
        case _ => sys.error("wrong args")
      }
    }
  }

  def variance(xs: Iterable[Double]) = {
    val mean = xs.sum / xs.size
    xs.map(x => Math.pow(x-mean, 2)).sum / xs.size
  }

  def stdDev(xs: Iterable[Double]) = Math.sqrt(variance(xs))

}
