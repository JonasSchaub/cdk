package org.openscience.cdk;

import org.junit.jupiter.api.Disabled;
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
import java.util.*;

/**
 * playground for deep testing
 */
public class BiosynfoniFingerprintTest {

    public String filePathCSVout = "src/test/resources/data/cdd.csv";
    public String filePathCSVout2 = "src/test/resources/data/cdd2.csv";
    public int Limit = 20000;
    private final SilentChemObjectBuilder chemObjectBuilder = new SilentChemObjectBuilder();

    private final SmilesParser smilesParser = new SmilesParser(chemObjectBuilder);

    private final BiosynfoniFingerprinter fingerprint = new BiosynfoniFingerprinter(true, false);

    @Disabled
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
            String smiles = smilesGenerator.create(aMolecule);
            IAtomContainer anotherMolecule = smilesParser.parseSmiles(smiles);
            ICountFingerprint fp = fingerprint.getCountFingerprint(aMolecule);

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

    @Test
    void testCSVCreation() throws Exception {
        String[] smiles = {"O=C(O)C1OC(OC2C(O)CC3(C)C(CCC4(C)C3CC=C5C6CC(C)(C)CCC6(C(=O)OC7OC(CO)C(O)C(O)C7O)CCC54C)C2(C)C)C(O)C(OC8OC(CO)C(O)C(OC9OCC(O)C(O)C9O)C8O)C1O",
                "O=C1OC2=CC(=C3C(OC=C3C)=C2C(=C1CC(=O)N4CCN(CCC=5C=CN=CC5)CC4)C)C",
                "O=C1OC2CCC3(O)C1C(C2)C45C(OC)CCC6(C)CN(CC)C5C3C(O)C64",
                "O=C(OC1CC2C(=C)C(O)C3(C2)C1C45COC3(O)C(O)C5C(C)(C)CCC4O)C",
                "O=C(OC1C2C(=C)C(O)C31C4OC(OC)C5(CCCC(C)(CO)C5C4)C3CC2)C",
                "OC1=CC=C(C=C1O)C2OCC(C(O)C3=CC=C(O)C(OC)=C3)C2CO",
                "OCC1(C)CCCC23COC(O)(C(O)C12)C45CC(C(=C)C4O)CC(O)C35",
                "O=C(OCC(=CC(C=C)C(=C)C)C)C(C)CC",
                "OCC1OC(OC2=CC=C3C(OCC4C5=CC=6OCOC6C=C5OC34)=C2)C(O)C(O)C1O",
                "OC1C=C2C3CC(C)(C)CCC3(C)CCC2(C)C4(C)CCC5C(C)(C)C(O)CCC5(C)C14",
                "O=C1OC2C3=C(C)CCC(O)C3(C)CC(OC(=O)C4(OC4C)C)C2C1=C",
                "O=C1C=2C(O)=CC(O)=CC2OC(C3=CC=C(O)C(O)=C3)C1OC4OC(C)C(O)C(OC5OC(O)C(O)C(O)C5O)C4O",
                "O=C(O)C1(C)C(OC2OC(CO)C(O)C(O)C2NC(=O)C)CCC3(C4=C(CCC31)C5CCC(C(C)C(O)C(O)C(C)C(C)C)C5(C)CC4)C",
                "O=C(O)C(N1C(=O)C=2C=CC=CC2C1)C(C)CC",
                "O=C1C2=C(O)C=3C(OC4OC(C(O)C(O)C4O)C5(O)C=CC(N)CC5C(SSC(CO)C1)CC=6C=CC=CC6)=CC(OC)=C(C3C=C2C)CO",
                "O=C(OC1C(=O)C2=C(C(=O)CC3C(C(=O)CCC23C)(C)C)C4(C(=O)CC(C(C)CC(=O)CC(C(=O)O)C)C14C)C)C",
                "O=C1C2=C(C=NN1CC(=O)NCCC3=CC=C(OC)C=C3)C=CC(OC)=C2OC",
                "O=C1C=2C=CC=CC2N3C(=O)CCC3(C(=O)NC4=NC=CS4)N1CCCC",
                "O=C(NC=1C=CC=C(OC)C1)C=2C(=O)N3C4=C(C=CC=C4CC3)C2O"};

        for (String smile : smiles) {
            IAtomContainer aMolecule = smilesParser.parseSmiles(smile);

            ICountFingerprint count = fingerprint.getCountFingerprint(aMolecule);

            createCSV(count, smile, 0, filePathCSVout2);

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

        for (int i = 0; i < countSize; i++) {
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

    @Test @Disabled
    /**
     * Tests the differnces between the orignial Biosynfoni fingerprint from python and the CDK based version
     * Method does not create the CSV Data
     */
    public void compareCSV() throws IOException {

        List<List<String>> thisFingerprint = loadCsv(filePathCSVout);
        List<List<String>> original0Fingerprint = loadCsv("C:\\Users\\micro\\IdeaProjects\\biosynfoni\\src\\dataPython.csv");
        List<List<String>> thisFingerprint2 = loadCsv(filePathCSVout2);
        int error = 0;
        int lastError = -1;
        for (int row = 0; row < Math.min(thisFingerprint.size(), Limit); row++) {

            List<String> thisCounts = thisFingerprint.get(row);
            List<String> originalCounts = original0Fingerprint.get(row);

            String thisSmiles = thisCounts.get(0);
            String originalSmiles = originalCounts.get(0);

            for (int i = 1; i < thisCounts.size(); i++) {

                if (!thisCounts.get(i).equals(originalCounts.get(i))) {
                    System.out.println(
                            "\nMolecule" + row + "\n" +
                                    "\nDifference at count: " + i +
                                    "\noriginalCount: " + originalCounts.get(i) +
                                    "\nthis Count: " + thisCounts.get(i) +
                                    "\nOccured in Molecule:" + thisSmiles
                    );
                    if(lastError != row ) {
                        error++;
                        lastError = row;
                    }
                }
            }
        }
        System.out.println("\nErrorCount; " + error);
    }

    /**
     * Loads a CSV into a List<List<String> format.
     * The outer List contains the rows
     * The inner List contains the content of a row
     *
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
        String testMol1 = "C1CCC2C(C1)O2";
        String testMol2 = "C1CC2C1O2";
        String testMol3 = "C1CC2C(C1)O2";


        String[] smarts = {
                "[#6]",
                "[#6;!$([r4])]",
                "[#6;$([r4])]"
        };

        BiosynfoniFingerprinter fp = new BiosynfoniFingerprinter(false, false,smarts);
        SmilesParser smilesParser1 = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer mol1 = smilesParser1.parseSmiles(testMol1);
        IAtomContainer mol2 = smilesParser1.parseSmiles(testMol2);
        IAtomContainer mol3 = smilesParser1.parseSmiles(testMol3);


        ICountFingerprint count2 = fp.getCountFingerprint(mol2);


        System.out.println("\n");
        ICountFingerprint count1 = fp.getCountFingerprint(mol1);
        System.out.println("\n");
        ICountFingerprint count3 = fp.getCountFingerprint(mol3);


        for (int i = 0; i<count2.size(); i++) {
            System.out.println(count1.getCount(i)+", "+count2.getCount(i)+", " + count3.getCount(i)+" "+ i);

        }
//
//        Cycles cycles = Cycles.sssr(mol1);
//        IRingSet rings = cycles.toRingSet();
//
//        for (int atomIndex = 0; atomIndex < mol1.getAtomCount(); atomIndex++) {
//
//            IAtom atom = mol1.getAtom(atomIndex);
//
//            System.out.println(
//                    "\nAtom " + atomIndex +
//                            " (" + atom.getSymbol() + ")"
//            );
//
//            int ringNumber = 1;
//
//            for (IAtomContainer ring : rings.atomContainers()) {
//
//                if (ring.contains(atom)) {
//
//                    System.out.println(
//                            " -> in Ring " + ringNumber +
//                                    " | Ringgröße: " + ring.getAtomCount()
//                    );
//                }
//
//                ringNumber++;
//            }
//        }

    }
}

