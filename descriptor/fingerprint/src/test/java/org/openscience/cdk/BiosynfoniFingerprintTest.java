package org.openscience.cdk;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.BiosynfoniFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;


import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.logging.Logger;

import static org.junit.jupiter.api.Assertions.*;

public class BiosynfoniFingerprintTest {

    public String filePathCSVout = "src/test/resources/data/cdd.csv";
    public int Limit = 100;
    private final SilentChemObjectBuilder chemObjectBuilder = new SilentChemObjectBuilder();

    private final SmilesParser smilesParser = new SmilesParser(chemObjectBuilder);

    private final BiosynfoniFingerprinter fingerprint = new BiosynfoniFingerprinter(true,false);

    @Disabled @Test
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
     *#TODO Search for smallest smile that has same error
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
                            "Molecule"+ row +
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
        String testMol1 = "C1CCC2C(C1)O2";
        String testMol2 = "C1CCCCC1";
        String testMol3 = "C1CCC(CC1)O";


        String[] smarts = {
                "[#6]",
                "[#6;!$([r6])]",
                "[#6;$([r6])]"
        };

        BiosynfoniFingerprinter  fp = new BiosynfoniFingerprinter(false,false, smarts);
        SmilesParser smilesParser1 = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer mol1 = smilesParser1.parseSmiles(testMol1);
        IAtomContainer mol2 = smilesParser1.parseSmiles(testMol2);
        IAtomContainer mol3 = smilesParser1.parseSmiles(testMol3);


        ICountFingerprint count2 = fp.getCountFingerprint(mol2);
        ICountFingerprint count1 = fp.getCountFingerprint(mol1);
        ICountFingerprint count3 = fp.getCountFingerprint(mol3);


        for (int i = 0; i<count2.size(); i++) {
            System.out.println(count2.getCount(i)+", "+count1.getCount(i)+", " + count3.getCount(i)+" "+ i);

        }

        Cycles cycles = Cycles.sssr(mol1);
        IRingSet rings = cycles.toRingSet();

        for (int atomIndex = 0; atomIndex < mol1.getAtomCount(); atomIndex++) {

            IAtom atom = mol1.getAtom(atomIndex);

            System.out.println(
                    "\nAtom " + atomIndex +
                            " (" + atom.getSymbol() + ")"
            );

            int ringNumber = 1;

            for (IAtomContainer ring : rings.atomContainers()) {

                if (ring.contains(atom)) {

                    System.out.println(
                            " -> in Ring " + ringNumber +
                                    " | Ringgröße: " + ring.getAtomCount()
                    );
                }

                ringNumber++;
            }
        }

    }
}

class BiosynfoniTest {
    String[] testSmiles = {"CC(=O)CC=O",
            "C1CCCCC1",
            "COc1cc(O)c2c(c1)oc(cc2=O)-c1ccc(OC)c(c1)-c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)cc1",
            "[H][C@]1(CC[C@@H](O)[C@@H](C1)OC)C[C@@H](C)[C@]1([H])CC(=O)[C@H](C)\\C=C(C)\\[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)\\C=C\\C=C\\C=C(C)\\[C@H](C[C@]2([H])CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCC[C@@]2([H])C(=O)O1)OC",
            "[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)*)CO[C@@H]1O[C@H](CO)[C@H]([C@@H]([C@H]1O)O)O[C@@H]2O[C@H](CO[C@]3(O[C@]([C@@H]([C@H](C3)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])C(=O)O)[C@@H]([C@H](O)[C@H]2O)O"};

    SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
@Test
void allSmartsViable(){
    BiosynfoniFingerprinter Bfp = new BiosynfoniFingerprinter(false,false);

    for(BiosynfoniFingerprinter.DefaultBiosynfoniKey key  : BiosynfoniFingerprinter.DefaultBiosynfoniKey.values()){
        try{
            SmartsPattern.create(key.smarts,
                    DefaultChemObjectBuilder.getInstance());
        }catch(Exception e){
            fail(" Invalid smiles detected" +
                      "Key: " + key.name() +
                      "Label: " + key.label +
                      "SMARTS: " + key.smarts +
                      "Error: " + e.getMessage());
        }
    }
}
@Test
void substructureDetectionSimple(){
    IAtomContainer smallMol = null;

    try{ smallMol = smilesParser.parseSmiles("CC");
} catch (InvalidSmilesException e) {
        fail("Error: " + e.getMessage()+ "this should not happen");
    }
    BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter();

    try {
        bfp.getBitFingerprint(smallMol);
        bfp.getCountFingerprint(smallMol);
    }catch (CDKException exception){
        fail("failed creating fingerprint"+"\n Error: "+ exception.getMessage());
    }

}
@Test
void substurctureDetection(){
    BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter();

    try{
        IAtomContainerSet aMoleculeSet = createMolecules();
        for (IAtomContainer aMolecule : aMoleculeSet){
            bfp.getCountFingerprint(aMolecule);
            bfp.getBitFingerprint(aMolecule);
        }
    } catch (CDKException e) {
        fail("failed creating fingerprint" +"\n Error: " + e.getMessage());
    }

}
@Test
void testNoChiralityDifference() {
    try{for(IAtomContainer mol : createMolecules()) {
        IAtomContainer noChrial = mol.clone();

        noChrial.setStereoElements(new ArrayList<>());

        BitSet noChiralFp = new BiosynfoniFingerprinter().getBitFingerprint(noChrial).asBitSet();
        BitSet chiralFp = new BiosynfoniFingerprinter().getBitFingerprint(mol).asBitSet();


        assertEquals(noChiralFp, chiralFp);}

    }catch(CDKException cdkException){
        fail("Error:"+  cdkException.getMessage());
    }catch(CloneNotSupportedException cloneException){
        fail("Error:"+  cloneException.getMessage());
    }
}



private IAtomContainerSet createMolecules() throws InvalidSmilesException {
    IAtomContainerSet molecules = new AtomContainerSet();
    for(String smile : testSmiles) {
        IAtomContainer aMolecule = smilesParser.parseSmiles(smile);
        molecules.addAtomContainer(aMolecule);
    }
    return molecules;
}


}

