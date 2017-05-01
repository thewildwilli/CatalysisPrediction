package testprogs

// Created by Ernesto on 27/05/2016.
import docking.{DockLog}
import docking.dockscore.SurfaceDistanceScorer
import docking.dockscore.Scorer
import docking.docksearch.FourInitialsDocker2D
import model.Molecule
import io._

object FourDocker {
  def main(args: Array[String]) {
    // Input 2 models, dock them, output them together. We might also need to output the entire rotate / translate operation.
    System.out.println("loading")
    val molA: Molecule = new Pdb2DReader(args(0)).read
    val molB: Molecule = new Pdb2DReader(args(1)).read
    val scorer: Scorer = new SurfaceDistanceScorer(0)
    val docker = new FourInitialsDocker2D(scorer)
    System.out.println(String.format("docking with docker %s, scorer %s ", docker.getClass.getName, scorer.toString))
    val finalState: Tuple2[Molecule, Double] = docker.dock(molA, molB, new DockLog(null, false))
    System.out.println("Final score: " + finalState._2)
    val molBDocked: Molecule = finalState._1
    // Now create a molecule with both A's and B's atoms
    System.out.println("outputting to " + args(2))
    val result: Molecule = molA.clone
    import scala.collection.JavaConversions._
    for (bAtom <- molBDocked.JAtoms.values) result.JAtoms.put(bAtom.id, bAtom)
    new XyzWriter(args(2)).write(result)
  }
}
