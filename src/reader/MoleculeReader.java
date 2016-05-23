package reader;

import model.Molecule;

import java.io.IOException;

// Created by Ernesto on 23/05/2016.
public interface MoleculeReader {
    Molecule read() throws IOException, ChemicalFormatException;
}
