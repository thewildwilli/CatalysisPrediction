package io;

import model.Molecule;

import java.io.IOException;

// Created by Ernesto on 23/05/2016.
public interface MoleculeWriter {
    void write(Molecule m) throws IOException;
}