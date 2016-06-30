package testprogs;

import model.*;
import io.ChemicalFormatException;
import io.MoleculeReader;
import io.PdbReader;

import java.io.IOException;

// Created by Ernesto on 23/05/2016.
public class PdbReaderMiniTest {
    public static void main(String[] args) throws IOException, ChemicalFormatException {
        MoleculeReader reader = new PdbReader(args[0]);
        Molecule m = reader.read();
        for (Atom a: m.JAtoms() )
            System.out.println(a.toString());
    }
}
