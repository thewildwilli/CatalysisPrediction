package io;// Created by Ernesto on 23/05/2016.

import model.Atom;
import model.Molecule;

import java.io.FileWriter;
import java.io.IOException;

public class Xyz2DWriter implements MoleculeWriter {
    private final String path;
    public Xyz2DWriter(String path){
        this.path = path;
    }

    /** Write to file using xyz format. Element name is not supported so C is always used.
     *  See: http://openbabel.org/wiki/XYZ_(format)*/
    @Override
    public void write(Molecule m) throws IOException {
        try (FileWriter fw = new FileWriter(this.path)) {
            fw.write(m.Atoms().size() + "\n");        // First line: Atom count
            fw.write(this.path  + "\n");            // Second line: title
            for (Atom a: m.JAtoms())
                fw.write(String.format("%s %.3f %.3f %.3f%n", a.element(), + a.x(), a.y(), 0.0));
        }
    }
}
