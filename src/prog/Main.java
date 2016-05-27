package prog;// Created by Ernesto on 26/05/2016.

import docking.dockscore.ConstantSurfaceScorer;
import docking.dockscore.Scorer;
import docking.docksearch.AtomPairDocker;
import io.ChemicalFormatException;
import io.MoleculeReader;
import io.Pdb2DReader;
import io.XyzWriter;
import model.Atom;
import model.Molecule;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException, ChemicalFormatException {
        // Input 2 models, dock them, output them together. We might also need to output the entire rotate / translate operation.
        System.out.println("loading");
        Molecule molA = new Pdb2DReader(args[0]).read();
        Molecule molB = new Pdb2DReader(args[1]).read();

        Scorer scorer = new ConstantSurfaceScorer(1.4);
        System.out.println("docking with scorer: " + scorer.toString());
        Molecule molBDocked = AtomPairDocker.dock(molA, molB, scorer).b();

        // Now create a molecule with both A's and B's atoms
        System.out.println("outputting to " + args[2]);
        Molecule result = molA.clone();
        result.setElement('O');
        for (Atom bAtom : molBDocked.JAtoms())
            result.JAtoms().add(bAtom);

        new XyzWriter(args[2]).write(result);
    }
}
