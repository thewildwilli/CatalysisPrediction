package io;

import model.Atom;
import model.Molecule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Pdb2DReader implements MoleculeReader {
    private final String path;
    public Pdb2DReader(String path){
        this.path = path;
    }

    /** Read file line by line, and only process ATOM lines.
     *  See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     *   */
    @Override
    public Molecule read() throws IOException, ChemicalFormatException {
        Molecule result = new Molecule();
        String line;
        int atomId = 1;
        try (BufferedReader br = new BufferedReader(new FileReader(this.path))) {
            while ((line = br.readLine()) != null) {
                // If line is an ATOM line, process it
                if (line.length() >= 6) {
                    String header = line.substring(0, 6);
                    if (header.equals("ATOM  ") || header.equals("HETATM")) {
                        if (line.length() < 54)
                            throw new ChemicalFormatException("ATOM record in PDB file does not have coordinates");
                        double x = Double.parseDouble(line.substring(30, 38).trim());
                        double y = Double.parseDouble(line.substring(38, 46).trim());
                        double z = 0.0;
                        result.JAtoms().put(atomId, new Atom(atomId, "C", x, y, z, 0.0, "", "", "", null));
                        atomId++;
                    }
                }
            }
        }
        return result;
    }
}
