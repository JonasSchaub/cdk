package org.openscience.cdk;

import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.BiosynfoniFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;


import java.io.*;
import java.nio.charset.StandardCharsets;

public class BiosynfoniFingerprintTest {
    private final SilentChemObjectBuilder chemObjectBuilder = new SilentChemObjectBuilder();

    private final SmilesParser smilesParser = new SmilesParser(chemObjectBuilder);

    private final BiosynfoniFingerprinter fingerprint = new BiosynfoniFingerprinter();

    @Test
    void testBiosynfoniFingerprint() throws CDKException, IOException {
        String filePath = "src/test/resources/data/cdd.csv";

        File file = new File(filePath);

        if (file.exists()) {
            file.delete();
        }
        Reader intput = new InputStreamReader(
                new FileInputStream("src/test/resources/data/coconut_sdf_2d_lite-05-2026.sdf"),
                StandardCharsets.UTF_8
        );
        IteratingSDFReader reader = new IteratingSDFReader(
                intput, SilentChemObjectBuilder.getInstance());
        SmilesGenerator smilesGenerator = new SmilesGenerator().unique();
        int Index = 0;
        int Limit = 100;
        while (reader.hasNext() && Index < Limit) {
            IAtomContainer aMolecule = reader.next();
            ICountFingerprint fp = fingerprint.getCountFingerprint(aMolecule);
            String smiles = smilesGenerator.create(aMolecule);
            createCSV(fp, smiles, Index, filePath);
            Index += 1;
        }
        IAtomContainer aMolecule = smilesParser.parseSmiles("CC(=O)O");

        IBitFingerprint fp = fingerprint.getBitFingerprint(aMolecule);

        for (int i = 0; i < fp.size(); i++) {
            if (fp.get(i)) {
                System.out.println(i);
            }

        }
    }

    /**
     * Writes the count fingerprint of a molecule to a CSV file.
     * The method appends one line per molecule containig:
     * -the molecule SMILES string
     * -all count fingperprint values
     * The CSV file is appended to incrementally, making it suitable
     * for large dataset processing
     *
     * @param countFingerprint the count fingerprint to export
     * @param smiles           the molecule identifier or SMILES string
     * @param index            of the processed molecule
     * @param filePath         path to the output CSV file defaults to "data/cdd" if null
     * @throws IOException  if writing to the file fails
     * @throws CDKException if fingerprint processing fails
     *
     */
    public void createCSV(ICountFingerprint countFingerprint, String smiles, int index, String filePath) throws CDKException, IOException {
        if (filePath == null) {
            filePath = "data/cdd";
        }
        StringBuilder header = new StringBuilder("name");
        int countSize = Math.toIntExact(countFingerprint.size()); // z.B. 39
        int[] countArray = new int[countSize];

        for (int i = 0; i < countFingerprint.numOfPopulatedbins(); i++) {
            countArray[i] = countFingerprint.getCount(i);
        }
        StringBuilder builder = new StringBuilder();

        smiles = smiles.replace("\"", "\"\""); // escaped quotes
        builder.append("\"").append(smiles).append("\"");

        for (int i = 0; i < countArray.length; i++) {
            builder.append(",").append(countArray[i]);
            header.append(",count_").append(i);
        }

        String newLine = builder.toString();
        File file = new File(filePath);
        boolean writeheader = !file.exists() || file.length() == 0;
        try (
                BufferedWriter writer = new BufferedWriter(new FileWriter(filePath, true))) {
            //header
            if (writeheader) {
                writer.write(String.valueOf(header));
                writer.newLine();
            }
            writer.write(newLine);
            writer.newLine();


        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }
}
