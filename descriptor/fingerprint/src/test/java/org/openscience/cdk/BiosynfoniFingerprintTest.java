package org.openscience.cdk;

import net.bytebuddy.implementation.bind.annotation.IgnoreForBinding;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.BiosynfoniFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;


import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

public class BiosynfoniFingerprintTest {

    public String filePathCSVout = "src/test/resources/data/cdd.csv";
    public int Limit = 50000;
    private final SilentChemObjectBuilder chemObjectBuilder = new SilentChemObjectBuilder();

    private final SmilesParser smilesParser = new SmilesParser(chemObjectBuilder);

    private final BiosynfoniFingerprinter fingerprint = new BiosynfoniFingerprinter();

    @Test
    void testBiosynfoniFingerprint() throws CDKException, IOException {


        File file = new File(filePathCSVout);

        if (file.exists()) {
            file.delete();
        }
        Reader intput = new InputStreamReader(
                new FileInputStream("C:\\Users\\micro\\IdeaProjects\\cdk\\descriptor\\fingerprint\\src\\test\\resources\\data\\coconut_sdf_2d_lite-05-2026.sdf"),
                StandardCharsets.UTF_8
        );
        IteratingSDFReader reader = new IteratingSDFReader(
                intput, SilentChemObjectBuilder.getInstance());
        SmilesGenerator smilesGenerator = new SmilesGenerator().unique();
        int Index = 0;

        while (reader.hasNext() && Index < Limit) {
            IAtomContainer aMolecule = reader.next();
            ICountFingerprint fp = fingerprint.getCountFingerprint(aMolecule);
            String smiles = smilesGenerator.create(aMolecule);
            createCSV(fp, smiles, Index, filePathCSVout);
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
            filePath = filePathCSVout;
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

    @Test
    /**
     * Tests the differnces between the orignial Biosynfoni fingerprint from python and the CDK based version
     * Method does not create the CSV Data
     */
    public void compareCSV() throws IOException {

        List<List<String>> thisFingerprint = loadCsv(filePathCSVout);
        List<List<String>> original0Fingerprint = loadCsv("C:\\Users\\micro\\IdeaProjects\\biosynfoni\\src\\dataPython.csv");
        int error = 0;

        for (int row = 0; row < Math.min(thisFingerprint.size(), Limit); row++) {

            List<String> thisCounts = thisFingerprint.get(row);
            List<String> originalCounts = original0Fingerprint.get(row);

            String thisSmiles = thisCounts.get(0);
            String originalSmiles = originalCounts.get(0);

            for (int i = 1; i < thisCounts.size(); i++) {

                if (!thisCounts.get(i).equals(originalCounts.get(i))) {
                    System.out.println(
                            "\nDifference at count: " + i +
                                    "\noriginalCount: " + originalCounts.get(i) +
                                    "\nthis Count: " + thisCounts.get(i) +
                                    "\nOccured in Molecule:" + thisSmiles + "\n" + originalSmiles
                    );
                    error++;
                }
            }
        }
        System.out.println("\nErrorCount; " + error);
    }

    /**
     * Loads a CSV into a List<List<String> format.
     * The outer List contains the rows
     * The inner List contains the content of a row
     * @param s filePath of the CSV to Load
     * @return List<List<String> containing each row of the CSV
     * @throws IOException
     */
    public List<List<String>> loadCsv(String s) throws IOException {
        List<List<String>> data = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(s));
        String row;

        while ((row = br.readLine()) != null) {

            String[] pieces = row.split(",");

            List<String> line = new ArrayList<>();

            for (String piece : pieces) {
                line.add(piece);
            }
            data.add(line);

        }
        br.close();

        return data;
    }
    @Test
    public void RingTest() throws CDKException {
        String testMol1 = "C1CNC2=CC=CC=C21";
        String testMol2 = "C(CC1)CCC1C12OC1CCCC2";


        BiosynfoniFingerprinter  fp = new BiosynfoniFingerprinter();
        SmilesParser smilesParser1 = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer mol1 = smilesParser1.parseSmiles(testMol1);
        IAtomContainer mol2 = smilesParser1.parseSmiles(testMol2);


        mol1 = fp.getAromaticity(mol1);
        mol2 = fp.getAromaticity(mol2);

        ICountFingerprint count2 = fp.getCountFingerprint(mol2);


        CycleFinder finder = Cycles.or(Cycles.relevant(), Cycles.vertexShort());

        for (int i = 0; i<count2.size(); i++) {
            System.out.println(count2.getCount(i)+","+i);
        }
        IRingSet rings2 =
                finder.find(mol2).toRingSet();
        System.out.println("\n New molecule      \n");
        for (IAtomContainer ring : rings2.atomContainers()) {

            System.out.println(
                    "Ring mit " +
                            ring.getAtomCount() +
                            " Atomen"
            );

        }

    }
}

