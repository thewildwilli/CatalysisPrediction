package io;

import model.*;

import java.io.*;

/**
 * Created by Ernesto on 23/05/2016.
 */
public class PdbReader implements MoleculeReader {
    private final String path;
    public PdbReader(String path){
        this.path = path;
    }

    /** Read file line by line, and only process ATOM lines.
     *  See http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
     *   */
    @Override
    public Molecule read() throws IOException, ChemicalFormatException {
        Molecule result = new Molecule();
        String line;
        try (BufferedReader br = new BufferedReader(new FileReader(this.path))) {
            while ((line = br.readLine()) != null) {
                // If line is an ATOM line, process it
                if (line.length() >= 4 && line.substring(0, 4).equals("ATOM")) {
                    if (line.length() < 54)
                        throw new ChemicalFormatException("ATOM record in PDB file does not have coordinates");
                    double x = Double.parseDouble(line.substring(30, 38).trim());
                    double y = Double.parseDouble(line.substring(38, 46).trim());
                    double z = Double.parseDouble(line.substring(46, 54).trim());
                    result.Atoms().add(new Atom(x, y, z));
                }
            }
        }
        return result;
    }
}
