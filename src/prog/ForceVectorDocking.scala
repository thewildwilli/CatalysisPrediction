package prog

import docking.{Docker, DockingState}
import docking.dockscore.{Scorer, SurfaceDistanceScorer}
import io.{Pdb2DReader, XyzWriter}
import model.{Atom, Molecule}
import opt.State

// Created by Ernesto on 08/06/2016.
class ForceVectorDocking {
  def main(args: Array[String]) {
    // Input 2 models, dock them, output them together. We might also need to output the entire rotate / translate operation.
    /*System.out.println("ForceVectorDocking...")
    val molA: Molecule = new Pdb2DReader(args(0)).read
    val molB = new Molecule();
    molB.Atoms.append(new Atom("O", 100, 100, 100))

    System.out.println("outputting initial state to " + args(1))
    new XyzWriter(args(1)).write(molA.clone.importM(molB));


    val finalState: State = docker.dock(molA, molB, scorer)
    System.out.println("Final score: " + scorer.score(finalState))
    val molBDocked: Molecule = (finalState.asInstanceOf[DockingState]).b

    // Now create a molecule with both A's and B's atoms
    System.out.println("outputting to " + args(2))
    val result: Molecule = molA.clone
    result.setElement("O")
    import scala.collection.JavaConversions._
    for (bAtom <- molBDocked.JAtoms) result.JAtoms.add(bAtom)

    new XyzWriter(args(2)).write(result)*/
  }
}
