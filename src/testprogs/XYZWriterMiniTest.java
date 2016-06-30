package testprogs;// Created by Ernesto on 23/05/2016.

import io.ChemicalFormatException;
import io.MoleculeReader;
import io.PdbReader;
import io.XyzWriter;
import model.Molecule;

import java.io.IOException;

public class XYZWriterMiniTest {
    public static void main(String[] args) throws IOException, ChemicalFormatException {
        MoleculeReader reader = new PdbReader(args[0]);
        Molecule m = reader.read();
        new XyzWriter(args[0] + ".xyz").write(m);
    }
}
