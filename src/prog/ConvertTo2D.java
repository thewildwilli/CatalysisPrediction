package prog;// Created by Ernesto on 26/05/2016.

import io.*;
import model.Molecule;

import java.io.IOException;
import java.nio.file.Paths;

public class ConvertTo2D {
    public static void main(String[] args) throws IOException, ChemicalFormatException {
        String path = args[0];
        MoleculeReader reader = new Pdb2DReader(path);
        Molecule m = reader.read();
        if (path.endsWith(".pdb") || path.endsWith(".PDB"))
            path = path.replace(".pdb", ".xyz");
        else
            path += ".xyz";
        new XyzWriter(path).write(m);
    }
}
