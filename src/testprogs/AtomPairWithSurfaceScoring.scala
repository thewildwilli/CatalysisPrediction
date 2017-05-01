package testprogs

// Created by Ernesto on 26/05/2016.
import docking.Docker
import docking.dockscore.SurfaceDistanceScorer
import docking.dockscore.Scorer
import docking.docksearch.SurfaceAtomPairsWithTranslationDocker
import io.{Pdb2DReader, XyzWriter}
import model.Molecule

object AtomPairWithSurfaceScoring {
  def main(args: Array[String]) {
    // Input 2 models, dock them, output them together. We might also need to output the entire rotate / translate operation.
    System.out.println("loading")
    val molA: Molecule = new Pdb2DReader(args(0)).read
    val molB: Molecule = new Pdb2DReader(args(1)).read
    molA.computeSurfaceAtoms2D()
    molB.computeSurfaceAtoms2D()
    val scorer: Scorer = new SurfaceDistanceScorer(0)
    val docker: Docker = new SurfaceAtomPairsWithTranslationDocker(scorer)
    System.out.println("docking with scorer: " + scorer.toString)
    val molBDocked: Molecule = docker.dock(molA, molB, null)._1
    // Now create a molecule with both A's and B's atoms
    System.out.println("outputting to " + args(2))
    val result: Molecule = molA.clone
    import scala.collection.JavaConversions._
    for (bAtom <- molBDocked.JAtoms.values) result.JAtoms.put(bAtom.id, bAtom)
    new XyzWriter(args(2)).write(result)
  }
}
