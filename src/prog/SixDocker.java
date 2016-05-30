package prog;// Created by Ernesto on 27/05/2016.

import docking.Docker;
import docking.dockscore.SurfaceDistanceScorer;
import docking.dockscore.Scorer;
import docking.docksearch.SixInitialsDocker$;
import io.ChemicalFormatException;
import io.Pdb2DReader;
import io.XyzWriter;
import model.Atom;
import model.Molecule;

import java.io.IOException;

public class SixDocker {
    public static void main(String[] args) throws IOException, ChemicalFormatException {
        // Input 2 models, dock them, output them together. We might also need to output the entire rotate / translate operation.
        System.out.println("loading");
        Molecule molA = new Pdb2DReader(args[0]).read();
        Molecule molB = new Pdb2DReader(args[1]).read();

        Scorer scorer = new SurfaceDistanceScorer();
        Docker docker = SixInitialsDocker$.MODULE$;
        System.out.println(String.format("docking with docker %s, scorer %s ", docker.getClass().getName(), scorer.toString()));
        Molecule molBDocked = docker.dock(molA, molB, scorer).b();

        // Now create a molecule with both A's and B's atoms
        System.out.println("outputting to " + args[2]);
        Molecule result = molA.clone();
        result.setElement("O");
        for (Atom bAtom : molBDocked.JAtoms())
            result.JAtoms().add(bAtom);

        new XyzWriter(args[2]).write(result);
    }
}
