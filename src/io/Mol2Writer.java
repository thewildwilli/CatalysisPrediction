package io;// Created by Ernesto on 23/05/2016.

import model.Atom;
import model.Molecule;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class Mol2Writer implements MoleculeWriter {
    private final String path;
    public Mol2Writer(String path){
        this.path = path;
    }

    /** Write to file using Mol2 format. It has the advantage of supporting charge.
     *  See: http://www.tripos.com/data/support/mol2.pdf */
    @Override
    public void write(Molecule m) throws IOException {
        try (FileWriter fw = new FileWriter(this.path)) {
            fw.write("@<TRIPOS>MOLECULE\n");
            fw.write(new File(path).getName() + "\n");
            fw.write(String.format(" %d %n", m.JAtoms().size()));
            fw.write("SMALL\n");
            fw.write("CHARGES\n");

            int serial = 1;
            fw.write("@<TRIPOS>ATOM\n");
            for (Atom a: m.JAtoms())
                fw.write(String.format(
                        " %d %s %.8f %.8f %.8f %.4s 1 A %.8f%n",
                        serial++,
                        a.element(),
                        a.x(), a.y(), a.z(),
                        a.element(),
                        a.partialCharge()
                ));
        }
    }
}
