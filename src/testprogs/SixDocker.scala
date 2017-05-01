package testprogs

// Created by Ernesto on 27/05/2016.
import docking.Docker
import docking.dockscore.SurfaceDistanceScorer
import docking.dockscore.Scorer
import docking.docksearch.SixInitialsDocker
import model.Molecule

import io.{PdbReader, XyzWriter}

object SixDocker {
  def main(args: Array[String]) {
    // Input 2 models, dock them, output them together. We might also need to output the entire rotate / translate operation.
    System.out.println("loading")
    val molA: Molecule = new PdbReader(args(0)).read
    val molB: Molecule = new PdbReader(args(1)).read
    val scorer: Scorer = new SurfaceDistanceScorer(1.4)
    val docker: Docker = new SixInitialsDocker(scorer)
    System.out.println(String.format("docking with docker %s, scorer %s ", docker.getClass.getName, scorer.toString))
    val finalState: Tuple2[Molecule, Double] = docker.dock(molA, molB, null)
    System.out.println("Final score: " + finalState._1)
    val molBDocked: Molecule = finalState._1
    // Now create a molecule with both A's and B's atoms
    System.out.println("outputting to " + args(2))
    val result: Molecule = molA.clone
    import scala.collection.JavaConversions._
    for (bAtom <- molBDocked.JAtoms.values) result.JAtoms.put(bAtom.id, bAtom)
    new XyzWriter(args(2)).write(result)
  }
}
