/*
 * Copyright (c) 2026 Marlon Raffelt (git@marlon.raffelt.email)
 *
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version. All we ask is that proper credit is given
 * for our work, which includes - but is not limited to - adding the above
 * copyright notice to the beginning of your source code files, and to any
 * copyright notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 U
 */
package org.openscience.cdk.fingerprint;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.BitSet;

/**
 * Test class for {@link BiosynfoniFingerprinter}.
 * The class uses modified copies from <a href="https://github.com/lucinamay/biosynfoni/blob/main/tests/main_test.py">
 * original python implementation tests</a> , {@link SubstructureFingerprinterTest} and own additions.
 * In alot of test substructure key 13,14,15 are tested along other keys, because during development issues occurred,
 * in ring size perception.
 * @author Marlon Raffelt (MaRa1778)
 */
class BiosynfoniFingerprinterTest {
    //Following tests are modified copies of the test included in python implementation
    /**
     * Smiles are from the test used by the <a href="https://github.com/lucinamay/biosynfoni/blob/main/tests/main_test.py">
     * original python implementation</a>
     */
    private static final String[] TEST_SMILES = {"CC(=O)CC=O",
            "C1CCCCC1",
            "COc1cc(O)c2c(c1)oc(cc2=O)-c1ccc(OC)c(c1)-c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)cc1",
            "[H][C@]1(CC[C@@H](O)[C@@H](C1)OC)C[C@@H](C)[C@]1([H])CC(=O)[C@H](C)\\C=C(C)\\[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)\\C=C\\C=C\\C=C(C)\\[C@H](C[C@]2([H])CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCC[C@@]2([H])C(=O)O1)OC",
            "[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)*)CO[C@@H]1O[C@H](CO)[C@H]([C@@H]([C@H]1O)O)O[C@@H]2O[C@H](CO[C@]3(O[C@]([C@@H]([C@H](C3)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])C(=O)O)[C@@H]([C@H](O)[C@H]2O)O"};

    /**
     * Smiles Parser used in nearly every method.
     */
    private final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

    /**
     * General {@link BiosynfoniFingerprinter} with no overlap filtering active.
     */
    private final BiosynfoniFingerprinter biosynfoniFingerprinter = new BiosynfoniFingerprinter();

    /**
     * Verifies that {@link BiosynfoniFingerprinter} uses viable SMARTS for substructure detection
     */
    @Test
    void allSmartsViable() {
        for (BiosynfoniFingerprinter.BiosynfoniKey key : BiosynfoniFingerprinter.BiosynfoniKey.values()) {
            Assertions.assertDoesNotThrow(() -> {
                SmartsPattern.create(key.getSmarts(),
                        DefaultChemObjectBuilder.getInstance());
            }, " Invalid SMARTS detected" +
                    "Key: " + key.name() +
                    "Label: " + key.getLabel() +
                    "SMARTS: " + key.getSmarts());
        }
    }

    /**
     * Verifies that {@link BiosynfoniFingerprinter} fingerprinter can generate both bit and
     * count fingerprints for a simple molecule.
     *
     * <p>The test uses ethane ({@code CC}) as a minimal input and verifies
     * that fingerprint generation completes successfully.</p>
     */
    @Test
    void substructureDetectionSimple() throws CDKException {
        IAtomContainer smallMol = this.smilesParser.parseSmiles("CC");
        Assertions.assertNotNull(this.biosynfoniFingerprinter.getBitFingerprint(smallMol));
        Assertions.assertNotNull(this.biosynfoniFingerprinter.getCountFingerprint(smallMol));
        //check for ethyl key not to be null
        Assertions.assertNotNull(this.biosynfoniFingerprinter.getBitFingerprint(smallMol).get(19));
    }

    /**
     * Verifies that the Biosynfoni fingerprinter can generate both bit and
     * count fingerprints for a representative set of molecules.
     *
     * <p>The test iterates over the molecules returned by
     * {@code createMolecules()} and verifies that fingerprint generation
     * completes without throwing an exception.</p>
     */
    @Test
    void substructureDetection() throws CDKException {
        IAtomContainerSet aMoleculeSet = this.createMolecules();
        for (IAtomContainer aMolecule : aMoleculeSet) {
            Assertions.assertDoesNotThrow(() -> this.biosynfoniFingerprinter.getCountFingerprint(aMolecule));
            Assertions.assertDoesNotThrow(() -> this.biosynfoniFingerprinter.getBitFingerprint(aMolecule));
        }
    }

    /**
     * Verifies that stereochemical information does not affect the {@link BiosynfoniFingerprinter} fingerprint.
     *
     * <p>For each molecule in the test set, all stereochemical annotations are
     * removed and the resulting binary fingerprint is compared with the
     * fingerprint of the original molecule. The fingerprints are expected to be
     * identical, demonstrating that the Biosynfoni fingerprinter is
     * stereochemistry-independent.</p>
     */
    @Test
    void testNoChiralityDifference() throws CDKException, CloneNotSupportedException {

        for (IAtomContainer mol : createMolecules()) {
            IAtomContainer noChiral = mol.clone();

            noChiral.setStereoElements(new ArrayList<>());

            BitSet noChiralFp = new BiosynfoniFingerprinter().getBitFingerprint(noChiral).asBitSet();
            BitSet chiralFp = new BiosynfoniFingerprinter().getBitFingerprint(mol).asBitSet();

            Assertions.assertEquals(noChiralFp, chiralFp);
        }
    }

    /**
     * Creates a set of test molecules from the predefined SMILES strings.
     *
     * @return an {@link IAtomContainerSet} containing all molecules parsed from
     *         {@code testSmiles}
     * @throws InvalidSmilesException if any SMILES string in {@code testSmiles}
     *         cannot be parsed
     */
    private IAtomContainerSet createMolecules() throws InvalidSmilesException {
        IAtomContainerSet molecules = new AtomContainerSet();
        for (String smile : TEST_SMILES) {
            IAtomContainer aMolecule = this.smilesParser.parseSmiles(smile);
            molecules.addAtomContainer(aMolecule);
        }
        return molecules;
    }

    /**
     * Verifies that the {@link BiosynfoniFingerprinter} fingerprint has the right size.
     */
    @Test
    void testSize() {
        Assertions.assertEquals(39, this.biosynfoniFingerprinter.getSize());
    }

    //following Tests are modified Copies of the test included in substructureFingerprinter

    /**
     * Verifies that {@link BiosynfoniFingerprinter} bit fingerprinter correctly detects selected
     * functional group features.
     *
     * <p>The test uses the same reference molecule as the
     * {@code testFunctionalGroupsBinary()} test in the CDK
     * {@code SubstructureFingerprinterTest}. It verifies the expected presence
     * and absence of selected fingerprint features in the resulting binary
     * fingerprint.</p>
     */
    @Test
    void testFunctionalGroupsBinary() throws CDKException {

        IAtomContainer mol1 = this.smilesParser.parseSmiles("c1ccccc1CCC");
        IBitFingerprint bitfp = this.biosynfoniFingerprinter.getBitFingerprint(mol1);

        //molecule should match a phenylmotif with 3,2,1 side chain carbons and substructure for ring size 6
        Assertions.assertNotNull(bitfp);
        Assertions.assertTrue(bitfp.get(13));
        Assertions.assertTrue(bitfp.get(14));
        Assertions.assertTrue(bitfp.get(15));
        Assertions.assertTrue(bitfp.get(34));
        //molcule should not match CoA
        Assertions.assertFalse((bitfp.get(1)));
    }

    /**
     * Verifies that {@link BiosynfoniFingerprinter} count fingerprinter correctly detects selected
     * functional group features.
     *
     * <p>The test uses the same reference molecule as the
     * {@code testFunctionalGroupsBinary()} test in the CDK
     * {@code SubstructureFingerprinterTest}. It verifies the expected presence
     * and absence of selected fingerprint features in the resulting binary
     * fingerprint.</p>
     */
    @Test
    void testFunctionalGroupsCount() throws CDKException {

        IAtomContainer mol1 = this.smilesParser.parseSmiles("c1ccccc1CCC");
        ICountFingerprint cfp = this.biosynfoniFingerprinter.getCountFingerprint(mol1);

        //molecule should match a phenylmotif with 3,2,1 side chain carbons and substructure for ring size 6
        Assertions.assertNotNull(cfp);
        Assertions.assertEquals(1, cfp.getCountForHash(13));
        Assertions.assertEquals(1, cfp.getCountForHash(14));
        Assertions.assertEquals(1, cfp.getCountForHash(15));
        Assertions.assertEquals(1, cfp.getCountForHash(34));
        //molcule should not match CoA
        Assertions.assertEquals(0, cfp.getCountForHash(1));
    }

    /**
     * Verifies that {@link BiosynfoniFingerprinter} fingerprint correctly detects ring-related
     * fingerprint features.
     *
     * <p>The test uses the same reference molecule as the
     * {@code testRingsBinary()} test in the CDK
     * {@code SubstructureFingerprinterTest}. It verifies the expected presence
     * and absence of selected ring-related fingerprint features in the
     * resulting binary fingerprint.</p>
     */
    @Test
    void testRingsBinary() throws CDKException {

        IAtomContainer mol1 = this.smilesParser.parseSmiles("C(C1C2CCC2)C1(C1)C2(CCCC2)CCC1C1CCCCCC1");
        IBitFingerprint bitfp = this.biosynfoniFingerprinter.getBitFingerprint(mol1);
        //molcule contains rings the size of 3,4,5,6,7 and should match the respective keys 31,32,33,34,35
        Assertions.assertNotNull(bitfp);
        Assertions.assertTrue(bitfp.get(31));
        Assertions.assertTrue(bitfp.get(32));
        Assertions.assertTrue(bitfp.get(33));
        Assertions.assertTrue(bitfp.get(34));
        Assertions.assertTrue(bitfp.get(35));
        // bigger ring sizes shall not match like 36 (8 ring),37(9 ring),38 (10 ring)
        Assertions.assertFalse(bitfp.get(36));
        Assertions.assertFalse(bitfp.get(37));
        Assertions.assertFalse(bitfp.get(38));

    }

    /**
     * Verifies that the {@link BiosynfoniFingerprinter} fingerprint correctly identifies
     * aromatic and non-aromatic fingerprint features.
     * <p>The SMILES are taken from the
     * {@code testCountableMACCSBinary2} test in the CDK Substructure Fingerprinter
     * test suite, and 2 removed aliphatic nitrogen. The test verifies that selected Biosynfoni fingerprint features
     * are counted correctly.</p>
     */
    @Test
    void testAromaticityBinary() throws CDKException {

        IAtomContainer mol1 = this.smilesParser.parseSmiles("CCc1c[nH]c2cc(-c3ccc(CC)cc3)ccc12");
        IBitFingerprint bitfp = this.biosynfoniFingerprinter.getBitFingerprint(mol1);
        //molecule shall match patterns with defined aromaticity like 9 or 10 and not
        Assertions.assertNotNull(bitfp);
        Assertions.assertTrue(bitfp.get(9));
        Assertions.assertTrue(bitfp.get(10));
        //molecule shall not match aliphatic substructures
        Assertions.assertFalse(bitfp.get(27));

        IAtomContainer mol2 = this.smilesParser.parseSmiles("C1=C(NC=N1)CC(C(=O)O)N");
        IBitFingerprint bfp2 = this.biosynfoniFingerprinter.getBitFingerprint(mol2);
        //checks if aromcity detection for aminoacid smarts works properly
        Assertions.assertNotNull(bfp2);
        Assertions.assertTrue(bfp2.get(3));
        Assertions.assertFalse(bfp2.get(4));
    }

    /**
     * Verifies that the {@link BiosynfoniFingerprinter} fingerprint correctly identifies
     * fingerprint features for a non-proteinogenic amino acid.
     * The SMILES is Tryptophan made with MolView
     */
    @Test
    void testNonStandardAminoAcidsBinary() throws Exception {

        IAtomContainer mol1 = this.smilesParser.parseSmiles("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O");
        IBitFingerprint bitfp = this.biosynfoniFingerprinter.getBitFingerprint(mol1);
        //check for proper detection of amino acid substructures
        Assertions.assertNotNull(bitfp);
        Assertions.assertTrue(bitfp.get(3));
        Assertions.assertFalse(bitfp.get(4));
    }

    /**
     * Verifies that the {@link BiosynfoniFingerprinter} bit fingerprint reports the expected
     * feature counts for a reference molecule.
     *
     * <p>The SMILES are taken from the
     * {@code testCountableMACCSBinary2} test in the CDK Substructure Fingerprinter
     * test suite.</p>
     */
    @Test
    void testRightBits() throws Exception {

        IAtomContainer aMolecule = this.smilesParser.parseSmiles("C([S](O)(=O)=O)C1=C(C=CC=C1)CCCC[N+](=O)[O-]");
        IBitFingerprint bitfp = this.biosynfoniFingerprinter.getBitFingerprint(aMolecule);

        //molecule shall match substructure motif phenylethyl, isoprene and ethyl,
        //because it contains a ring with 2 attached side chains
        Assertions.assertTrue(bitfp.get(13));
        Assertions.assertTrue(bitfp.get(16));
        Assertions.assertTrue(bitfp.get(19));
        Assertions.assertTrue(bitfp.get(34));
        //molecule contains a sulfate group and should match substructure key 22 (Sulfonate) and key 30 (hydroxyl)
        Assertions.assertTrue(bitfp.get(22));
        Assertions.assertTrue(bitfp.get(30));
        //molecule does not contain a nitrogen close to the ring
        Assertions.assertFalse(bitfp.get(12));
        //ring size is not 8 in the given mol
        Assertions.assertFalse(bitfp.get(35));
    }

    /**
     * Verifies that the {@link BiosynfoniFingerprinter} count fingerprint reports the expected
     * feature counts for a reference molecule.
     *
     * <p>The SMILES are taken from the
     * {@code testCountableMACCSBinary2} test in the CDK Substructure Fingerprinter
     * test suite. The test verifies that selected Biosynfoni fingerprint features
     * are counted correctly.</p>
     */
    @Test
    void testRightCounts() throws CDKException {

        IAtomContainer aMolecule = this.smilesParser.parseSmiles("C([S](O)(=O)=O)C1=C(C=CC=C1)CCCC[N+](=O)[O-]");
        ICountFingerprint cfp = this.biosynfoniFingerprinter.getCountFingerprint(aMolecule);
        //molecule contains multiple carbons connected to each other
        Assertions.assertEquals(11, cfp.getCountForHash(19));
        //molecule contains one phenyethyl motif
        Assertions.assertEquals(1, cfp.getCountForHash(13));
        //molecule contains 7 not 2 isoprene substructures
        Assertions.assertNotEquals(2, cfp.getCountForHash(16));
    }

//own test additions

    /**
     * Verifies that {@link BiosynfoniFingerprinter} fingerprints can be used
     * with the CDK {@link Tanimoto} similarity calculation.
     *
     * <p>The test computes both count and bit fingerprints for two structurally
     * related molecules represented as SMILES strings and verifies that the
     * resulting Tanimoto similarity scores are valid (i.e., not {@code NaN}).</p>
     *
     * <p>SMILES correspond to:</p>
     * <ul>
     *   <li><a href="https://pubchem.ncbi.nlm.nih.gov/compound/445858">Ferulic acid</a></li>
     *   <li><a href="https://pubchem.ncbi.nlm.nih.gov/compound/8468">Vanillic acid</a></li>
     * </ul>
     */
    @Test
    void testSimilarityScores() throws CDKException {

        IAtomContainer aMolecule = this.smilesParser.parseSmiles("COC1=C(C=CC(=C1)/C=C/C(=O)O)O"); //ferulic acid
        IAtomContainer aMolecule2 = this.smilesParser.parseSmiles("COC1=C(C=CC(=C1)C(=O)O)O"); //Vanillic acid

        ICountFingerprint cfp1 = this.biosynfoniFingerprinter.getCountFingerprint(aMolecule);
        IBitFingerprint bitfp1 = this.biosynfoniFingerprinter.getBitFingerprint(aMolecule);

        ICountFingerprint cfp2 = this.biosynfoniFingerprinter.getCountFingerprint(aMolecule2);
        IBitFingerprint bitfp2 = this.biosynfoniFingerprinter.getBitFingerprint(aMolecule2);

        double countScore = Tanimoto.calculate(cfp1, cfp2);
        double bitScore = Tanimoto.calculate(bitfp1, bitfp2);

        Assertions.assertFalse(Double.isNaN(countScore));
        Assertions.assertFalse(Double.isNaN( bitScore));
    }

    /**
     * Verifies that enabling intra- and inter-substructure filtering does not
     * increase fingerprint feature counts. Filtering is expected to reduce or
     * maintain the counts.
     *
     * <p>The test uses
     * <a href="https://coconut.naturalproducts.net/compounds/CNP0175721.3">Rapamycin</a>
     * as the test molecule.</p>
     */
    @Test
    void filterTest() throws CDKException {
        IAtomContainer mol = this.smilesParser.parseSmiles("CO[C@@H]1C[C@H]2CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCCC2C(=O)O[C@@H]([C@@H](C)C[C@@H]2CC[C@@H](O)[C@H](OC)C2)CC(=O)[C@@H](C)/C=C(/C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)/C=C/C=C/C=C\\1C");

        BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(false, false);
        BiosynfoniFingerprinter bfp2 = new BiosynfoniFingerprinter(true, false);
        BiosynfoniFingerprinter bfp3 = new BiosynfoniFingerprinter(false, true);
        BiosynfoniFingerprinter bfp4 = new BiosynfoniFingerprinter(true, true);

        ICountFingerprint cIntraFilter = bfp2.getCountFingerprint(mol);
        ICountFingerprint cUnfiltered = bfp.getCountFingerprint(mol);
        ICountFingerprint cInterFilter = bfp3.getCountFingerprint(mol);
        ICountFingerprint cBothFilters = bfp4.getCountFingerprint(mol);

        for (int i = 0; i < bfp.getSize(); i++) {
            Assertions.assertTrue(cUnfiltered.getCountForHash(i) >= cIntraFilter.getCountForHash(i));
            Assertions.assertTrue(cUnfiltered.getCountForHash(i) >= cInterFilter.getCountForHash(i));
            Assertions.assertTrue(cUnfiltered.getCountForHash(i) >= cBothFilters.getCountForHash(i));
        }
    }

    /**
     * Verifies that {@link BiosynfoniFingerprinter} fingerprint detects aromatic indoles correctly.
     * First checks are to see expected behavior and second checks are to see if SMILES encoding matters.
     * mol2 is the provided SMILES for <a href="https://pubchem.ncbi.nlm.nih.gov/compound/L-Tryptophan">L-Tryptophan</a>
     * mol1 is self changed to ensure aromaticity in molecule.
     */
    @Test
    void testIndoleSubstructure() throws CDKException {
        IAtomContainer mol1 = this.smilesParser.parseSmiles("c1ccc2c(c1)c(c[nH]2)C[C@@H](C(=O)O)N");
        IAtomContainer mol2 = this.smilesParser.parseSmiles("C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N");

        ICountFingerprint cfp = this.biosynfoniFingerprinter.getCountFingerprint(mol1);
        ICountFingerprint cfp2 = this.biosynfoniFingerprinter.getCountFingerprint(mol2);

        Assertions.assertEquals(1, cfp.getCountForHash(9));
        Assertions.assertEquals(1, cfp.getCountForHash(10));
        Assertions.assertNotEquals(1, cfp.getCountForHash(11));
        Assertions.assertEquals(1, cfp.getCountForHash(12));

        Assertions.assertEquals(cfp2.getCountForHash(9), cfp.getCountForHash(9));
        Assertions.assertEquals(cfp2.getCountForHash(10), cfp.getCountForHash(10));
        Assertions.assertEquals(cfp2.getCountForHash(11), cfp.getCountForHash(11));
        Assertions.assertEquals(cfp2.getCountForHash(12), cfp.getCountForHash(12));
    }

    /**
     * Verifies that the canonical atom index mapping is generated correctly.
     *
     * <p>The test uses the phenol molecule and the expected canonical atom
     * index mapping derived from the CDK {@code Canon.label()} implementation.
     * The expected mapping is the inverse of the canonical labels and therefore
     * represents the canonical position of each atom in the original molecule.</p>
     */
    @Test
    void testCanoicalisation() throws InvalidSmilesException {
        IAtomContainer mol1 = this.smilesParser.parseSmiles("OC1=CC=CC=C1");
        BiosynfoniFingerprinter bfp3 = new BiosynfoniFingerprinter(false, true);
        int[] bfpCanonicalIndex = bfp3.getCanonicalIndexMap(mol1);
        int[] canonicalIndex = {0, 6, 4, 3, 1, 5, 2};
        Assertions.assertNotNull(bfpCanonicalIndex);
        Assertions.assertArrayEquals(bfpCanonicalIndex, canonicalIndex);
    }
}
