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
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;

/**
 * Test class for {@link BiosynfoniFingerprinter}.
 * The class uses modified copies from <a href="https://github.com/lucinamay/biosynfoni/blob/main/tests/main_test.py">
 * original python implementation tests</a> , {@link SubstructureFingerprinterTest} and own additions.
 * In a lot of test substructure key 13,14,15 are tested along other keys, because during development issues occurred,
 * in ring size perception. These keys are tested specifically in {@link #testPhenylKeysRingSizePerception()}.
 * //ToDo fix local build and tests
 * @author Marlon Raffelt (MaRa1778)
 */
class BiosynfoniFingerprinterTest {
    //Following tests are modified copies of the test included in the python implementation
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
     * Verifies that {@link BiosynfoniFingerprinter} uses viable original and refined SMARTS for substructure detection
     */
    @Test
    void allSmartsViable() {
        for (BiosynfoniFingerprinter.BiosynfoniKey key : BiosynfoniFingerprinter.BiosynfoniKey.values()) {
            for (String smarts : new String[]{key.getOriginalSmarts(), key.getRefinedSmarts()}) {
                Assertions.assertDoesNotThrow(() -> {
                    SmartsPattern.create(smarts,
                            DefaultChemObjectBuilder.getInstance());
                }, " Invalid SMARTS detected" +
                        "Key: " + key.name() +
                        "Label: " + key.getLabel() +
                        "SMARTS: " + smarts);
            }
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

        //molecule should match a phenyl motif with 3,2,1 side chain carbons and substructure for ring size 6
        Assertions.assertNotNull(bitfp);
        Assertions.assertTrue(bitfp.get(13));
        Assertions.assertTrue(bitfp.get(14));
        Assertions.assertTrue(bitfp.get(15));
        Assertions.assertTrue(bitfp.get(34));
        //molcule should not match CoA
        Assertions.assertFalse((bitfp.get(0)));
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
        Assertions.assertEquals(0, cfp.getCountForHash(0));
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
        //the original pattern of key 10 (phenethylamine) matches the aromatic pyrrole ring of the indole, the refined
        // one does not, see BiosynfoniKey.D_PHENYL_C2N
        Assertions.assertNotNull(bitfp);
        Assertions.assertFalse(bitfp.get(10));
        Assertions.assertTrue(new BiosynfoniFingerprinter(false, false, true).getBitFingerprint(mol1).get(10));
        //molecule shall not match aliphatic substructures
        Assertions.assertFalse(bitfp.get(27));

        IAtomContainer mol2 = this.smilesParser.parseSmiles("C1=C(NC=N1)CC(C(=O)O)N");
        IBitFingerprint bfp2 = this.biosynfoniFingerprinter.getBitFingerprint(mol2);
        //checks if aromaticity detection for amino acid smarts works properly
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
        //molecule contains a sulfonic acid group and should match substructure key 22 (Sulfonate)
        Assertions.assertTrue(bitfp.get(22));
        //but not key 30 (hydroxyl): the OH group of the sulfonic acid is bonded to sulfur, not to carbon or nitrogen,
        // and the [O-] of the charge-separated nitro group is not a deprotonated hydroxyl group (in the original
        // implementation, it is matched as hydroxyl, which was the only hydroxyl match here)
        Assertions.assertFalse(bitfp.get(30));
        //molecule does not contain a nitrogen close to the ring
        Assertions.assertFalse(bitfp.get(12));
        //ring size is not 7 in the given mol
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

    //the following test are additions specific for this implementation

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
        //the pyrrole ring of the indole is not a phenethylamine unit in the refined fingerprint, but in the original one
        Assertions.assertEquals(0, cfp.getCountForHash(10));
        Assertions.assertEquals(1, new BiosynfoniFingerprinter(false, false, true).getCountFingerprint(mol1).getCountForHash(10));
        Assertions.assertNotEquals(1, cfp.getCountForHash(11));
        Assertions.assertEquals(1, cfp.getCountForHash(12));

        Assertions.assertEquals(cfp2.getCountForHash(9), cfp.getCountForHash(9));
        Assertions.assertEquals(cfp2.getCountForHash(10), cfp.getCountForHash(10));
        Assertions.assertEquals(cfp2.getCountForHash(11), cfp.getCountForHash(11));
        Assertions.assertEquals(cfp2.getCountForHash(12), cfp.getCountForHash(12));
    }

    /**
     * Verifies the ring size perception of the phenyl keys 13 (phenylC3), 14 (phenylC2) and 15 (phenylC1).
     *
     * <p>The side chain carbons of keys 13 and 14 must not be part of any six-membered ring. In the original pattern
     * of the reference implementation, this is expressed with {@code [#6;!$([r6])]}, which only excludes atoms whose
     * smallest ring has six members. The refined pattern uses {@code [#6;!$(*~1~*~*~*~*~*~1)]} instead, which excludes
     * atoms in any six-membered ring (see {@link BiosynfoniFingerprinter.BiosynfoniKey#D_PHENYL_C3}).</p>
     *
     * <p>Molecules used:</p>
     * <ul>
     *   <li>{@code c1ccccc1CCC}: propylbenzene, all three keys match.</li>
     *   <li>{@code c1ccccc1CC12CCC1CC2}: benzyl bicyclo[2.2.0]hexane. The second side chain carbon is part of the
     *       six-membered perimeter ring, but its smallest ring has four members. The original pattern matches, the
     *       refined one must not.</li>
     *   <li>{@code Cc1cccc2ccccc12}: 1-methylnaphthalene. The fused ring atoms are part of two rings and
     *       therefore do not match {@code [#6;R1]} or {@code [c;R1]}, so no phenyl key matches.</li>
     * </ul>
     */
    @Test
    void testPhenylKeysRingSizePerception() throws CDKException {
        BiosynfoniFingerprinter originalBfp = new BiosynfoniFingerprinter(false, false, true);
        ICountFingerprint propylbenzene = this.biosynfoniFingerprinter.getCountFingerprint(
                this.smilesParser.parseSmiles("c1ccccc1CCC"));
        Assertions.assertEquals(1, propylbenzene.getCountForHash(13));
        Assertions.assertEquals(1, propylbenzene.getCountForHash(14));
        Assertions.assertEquals(1, propylbenzene.getCountForHash(15));

        IAtomContainer bicyclic = this.smilesParser.parseSmiles("c1ccccc1CC12CCC1CC2");
        ICountFingerprint bicyclicFp = this.biosynfoniFingerprinter.getCountFingerprint(bicyclic);
        Assertions.assertEquals(0, bicyclicFp.getCountForHash(13));
        Assertions.assertEquals(0, bicyclicFp.getCountForHash(14));
        Assertions.assertEquals(1, bicyclicFp.getCountForHash(15));
        //the original pattern does not exclude the side chain carbon in the six-membered ring
        Assertions.assertTrue(SmartsPattern.create(BiosynfoniFingerprinter.BiosynfoniKey.D_PHENYL_C2.getOriginalSmarts())
                .matches(bicyclic));
        Assertions.assertEquals(1, originalBfp.getCountFingerprint(bicyclic).getCountForHash(14));

        for (BiosynfoniFingerprinter bfp : new BiosynfoniFingerprinter[]{this.biosynfoniFingerprinter, originalBfp}) {
            ICountFingerprint methylnaphthalene = bfp.getCountFingerprint(this.smilesParser.parseSmiles("Cc1cccc2ccccc12"));
            Assertions.assertEquals(0, methylnaphthalene.getCountForHash(13));
            Assertions.assertEquals(0, methylnaphthalene.getCountForHash(14));
            Assertions.assertEquals(0, methylnaphthalene.getCountForHash(15));
        }
    }

    /**
     * Verifies that the canonical atom index mapping is generated correctly.
     *
     * <p>The test uses phenol in two different atom orders. The mapping returns the canonical
     * position (0 to 6) of every atom, indexed by the original atom index. The positions follow a
     * depth-first traversal along the bonds that starts at the oxygen atom (smallest canonical
     * label): oxygen 0, ipso carbon 1, then around the ring ortho, meta, para, meta, ortho. So,
     * independent of the atom order, the para carbon atom gets position 4, the two
     * symmetry-equivalent ortho carbon atoms the positions 2 and 6, and the two symmetry-equivalent
     * meta carbon atoms the positions 3 and 5. Which of two equivalent atoms gets which of the two
     * positions is not checked, because both atoms are interchangeable.</p>
     */
    @Test
    void testCanonicalisation() throws InvalidSmilesException {
        //the mapping only depends on the molecule, not on the overlap filter settings
        BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(true, false);
        //atom indices: O 0, ipso 1, ortho 2 and 6, meta 3 and 5, para 4
        int[] map1 = bfp.getCanonicalIndexMap(this.smilesParser.parseSmiles("OC1=CC=CC=C1"));
        BiosynfoniFingerprinterTest.assertPhenolMapping(map1, 0, 1, new int[]{2, 6}, new int[]{3, 5}, 4);
        //atom indices: para 0, meta 1 and 5, ortho 2 and 4, ipso 3, O 6
        int[] map2 = bfp.getCanonicalIndexMap(this.smilesParser.parseSmiles("C1=CC=C(C=C1)O"));
        BiosynfoniFingerprinterTest.assertPhenolMapping(map2, 6, 3, new int[]{2, 4}, new int[]{1, 5}, 0);
    }

    /**
     * Checks the canonical atom index mapping of phenol, see {@link #testCanonicalisation()}.
     *
     * @param map    the canonical atom index mapping, indexed by the original atom index
     * @param oxygen original index of the oxygen atom
     * @param ipso   original index of the ipso carbon atom
     * @param ortho  original indices of the ortho carbon atoms
     * @param meta   original indices of the meta carbon atoms
     * @param para   original index of the para carbon atom
     */
    private static void assertPhenolMapping(int[] map, int oxygen, int ipso, int[] ortho, int[] meta, int para) {
        Assertions.assertNotNull(map);
        Assertions.assertEquals(7, map.length);
        Assertions.assertEquals(0, map[oxygen]);
        Assertions.assertEquals(1, map[ipso]);
        Assertions.assertEquals(4, map[para]);
        int[] orthoPositions = {map[ortho[0]], map[ortho[1]]};
        int[] metaPositions = {map[meta[0]], map[meta[1]]};
        Arrays.sort(orthoPositions);
        Arrays.sort(metaPositions);
        Assertions.assertArrayEquals(new int[]{2, 6}, orthoPositions);
        Assertions.assertArrayEquals(new int[]{3, 5}, metaPositions);
    }

    /**
     * Verifies that the intra-substructure overlap filter selects the matches along chains and rings, so that they
     * are packed as densely as in the original Python implementation. The expected counts were taken from the original
     * implementation, where they do not depend on the input atom order for these molecules. Every molecule is given in
     * two atom orders.
     *
     * <ul>
     *   <li>p-benzoquinone: three disjoint ethyl units (key 19) in the six-membered ring and two acetyl units
     *   (key 17)</li>
     *   <li>thiophene 1,1-dioxide: two disjoint ethyl units in the chain of four carbon atoms; selecting the middle
     *   carbon-carbon bond first, as with a numbering directly by the canonical labels, would block both other bonds
     *   and result in only one ethyl unit</li>
     * </ul>
     */
    @Test
    void testIntraSubOverlapDensePacking() throws CDKException {
        BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(true, false);
        for (String smiles : new String[]{"O=C1C=CC(=O)C=C1", "C1=CC(=O)C=CC1=O"}) {
            ICountFingerprint cfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(smiles));
            Assertions.assertEquals(3, cfp.getCountForHash(19), smiles);
            Assertions.assertEquals(2, cfp.getCountForHash(17), smiles);
        }
        for (String smiles : new String[]{"O=S1(=O)C=CC=C1", "C1=CS(=O)(=O)C=C1"}) {
            ICountFingerprint cfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(smiles));
            Assertions.assertEquals(2, cfp.getCountForHash(19), smiles);
        }
    }

    /**
     * Verifies that the count fingerprint does not depend on the input atom order for all four combinations of the
     * overlap filters, for the refined and the original fingerprint. Every molecule is given in two different atom orders: as written and as unique SMILES generated
     * by the CDK. Without the canonical atom numbering, the intra-substructure overlap filter would select different
     * matches depending on the atom order, e.g. for the chain of four carbon atoms in thiophene 1,1-dioxide.
     *
     * <p>Molecules: p-benzoquinone, thiophene 1,1-dioxide,
     * <a href="https://coconut.naturalproducts.net/compounds/CNP0175721.3">Rapamycin</a>, and
     * <a href="https://coconut.naturalproducts.net/compounds/CNP0151033">CNP0151033</a>.</p>
     */
    @Test
    void testFilterAtomOrderIndependence() throws CDKException {
        String[] smilesArray = {
                "O=C1C=CC(=O)C=C1",
                "O=S1(=O)C=CC=C1",
                "CO[C@@H]1C[C@H]2CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCCC2C(=O)O[C@@H]([C@@H](C)C[C@@H]2CC[C@@H](O)[C@H](OC)C2)CC(=O)[C@@H](C)/C=C(/C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)/C=C/C=C/C=C\\1C",
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4"};
        SmilesGenerator smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);
        boolean[] toggles = {false, true};
        for (boolean isOriginalFingerprint : toggles) {
            for (boolean intraSubOverlapToggle : toggles) {
                for (boolean interSubOverlapToggle : toggles) {
                    BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(intraSubOverlapToggle, interSubOverlapToggle,
                            isOriginalFingerprint);
                    for (String smiles : smilesArray) {
                        IAtomContainer mol = this.smilesParser.parseSmiles(smiles);
                        IAtomContainer reorderedMol = this.smilesParser.parseSmiles(smilesGenerator.create(mol));
                        ICountFingerprint cfp = bfp.getCountFingerprint(mol);
                        ICountFingerprint reorderedCfp = bfp.getCountFingerprint(reorderedMol);
                        for (int i = 0; i < bfp.getSize(); i++) {
                            Assertions.assertEquals(cfp.getCountForHash(i), reorderedCfp.getCountForHash(i),
                                    "key " + i + ", intra " + intraSubOverlapToggle + ", inter " + interSubOverlapToggle
                                            + ", original " + isOriginalFingerprint + ": " + smiles);
                        }
                    }
                }
            }
        }
    }

    /**
     * Verifies that the counts are identical to the ones of the
     * <a href="https://github.com/lucinamay/biosynfoni">original Python implementation</a> without overlap filters
     * and with the inter-substructure overlap filter.
     *
     * <p>The expected counts were generated with the original implementation (version "full_1103", commit 0ba4620,
     * RDKit 2025.03.3, reproduced with RDKit 2025.03.6) via {@code Biosynfoni(mol, intersub_overlap=...).fingerprint},
     * with {@code intersub_overlap=True} for the counts without filter and {@code intersub_overlap=False} for the
     * counts with the inter-substructure overlap filter. In the original implementation, the parameters state whether
     * overlap is allowed, i.e. {@code False} enables the filter; {@code intrasub_overlap} keeps its default
     * {@code True}. The intra-substructure overlap filter is not compared because the counts of the original
     * implementation depend on the input atom order with this filter. The original fingerprint is used (see
     * {@link BiosynfoniFingerprinter#isOriginalFingerprint()}); the chosen molecules contain neither nitro groups,
     * nitrate esters, or N-oxides nor explicit hydrogen atoms. The refined fingerprint of the same molecules is tested
     * in {@link #testRefinedCounts()}.</p>
     *
     * <p>Molecules: <a href="https://pubchem.ncbi.nlm.nih.gov/compound/L-Tryptophan">L-Tryptophan</a>,
     * <a href="https://coconut.naturalproducts.net/compounds/CNP0175721.3">Rapamycin</a>,
     * <a href="https://coconut.naturalproducts.net/compounds/CNP0151033">CNP0151033</a>, and
     * <a href="https://pubchem.ncbi.nlm.nih.gov/compound/Coenzyme-A">Coenzyme A</a>.</p>
     */
    @Test
    void testReferenceCounts() throws CDKException {
        String[] smilesArray = {
                "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
                "CO[C@@H]1C[C@H]2CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCCC2C(=O)O[C@@H]([C@@H](C)C[C@@H]2CC[C@@H](O)[C@H](OC)C2)CC(=O)[C@@H](C)/C=C(/C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)/C=C/C=C/C=C\\1C",
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4",
                "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O"};
        int[][] expectedUnfiltered = {
                {0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 7, 2, 0, 11, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 1, 1, 17, 27, 14, 47, 10, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 1, 18, 29, 3, 29, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 6, 0, 0, 1, 1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 13, 5, 14, 2, 12, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0}};
        int[][] expectedInterFiltered = {
                {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 12, 2, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 15, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        BiosynfoniFingerprinter unfilteredBfp = new BiosynfoniFingerprinter(false, false, true);
        BiosynfoniFingerprinter interFilteredBfp = new BiosynfoniFingerprinter(false, true, true);
        this.assertCounts(smilesArray, expectedUnfiltered, expectedInterFiltered, unfilteredBfp, interFilteredBfp);
    }

    /**
     * Verifies the counts of the refined fingerprint (default) for the molecules of {@link #testReferenceCounts()},
     * without overlap filters and with the inter-substructure overlap filter.
     *
     * <p>Without filters, the counts differ from the original ones (see {@link #testReferenceCounts()}) as follows:
     * L-Tryptophan: no phenethylamine unit (key 10), since the pyrrole ring of the indole is not matched; Rapamycin:
     * the pipecolic acid residue, whose ring is matched as lysine side chain by both amino acid patterns, counts as
     * standard amino acid (key 3) instead of non-standard one (key 4), since its nitrogen atom is acylated and its
     * carboxy group esterified, the cyclohexane ring is no phenyl unit (keys 13 to 15), and the hemiketal hydroxyl
     * group is counted (key 30);
     * CNP0151033: the glucopyranose ring is no open-chain hexose or pentose (keys 5 and 6); Coenzyme A: three
     * phosphorus atoms instead of twelve phosphorus-oxygen bonds (key 21). The counts were cross-checked against the
     * pattern variants of the Biosynfoni fragmenter of MORTAR, which give the same counts.</p>
     */
    @Test
    void testRefinedCounts() throws CDKException {
        String[] smilesArray = {
                "C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
                "CO[C@@H]1C[C@H]2CC[C@@H](C)[C@@](O)(O2)C(=O)C(=O)N2CCCCC2C(=O)O[C@@H]([C@@H](C)C[C@@H]2CC[C@@H](O)[C@H](OC)C2)CC(=O)[C@@H](C)/C=C(/C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)/C=C/C=C/C=C\\1C",
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4",
                "CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCS)O"};
        int[][] expectedUnfiltered = {
                {0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 7, 2, 0, 11, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 17, 27, 14, 47, 10, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 1, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 18, 29, 3, 29, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 6, 0, 0, 1, 1, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 13, 5, 14, 2, 3, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0}};
        int[][] expectedInterFiltered = {
                {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 17, 2, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 15, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0},
                {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        this.assertCounts(smilesArray, expectedUnfiltered, expectedInterFiltered, new BiosynfoniFingerprinter(false, false),
                new BiosynfoniFingerprinter(false, true));
    }

    /**
     * Checks the count fingerprints of the given molecules generated with the given fingerprinters.
     *
     * @param smilesArray           SMILES codes of the molecules
     * @param expectedUnfiltered    expected counts of the first fingerprinter, per molecule
     * @param expectedInterFiltered expected counts of the second fingerprinter, per molecule
     * @param unfilteredBfp         first fingerprinter, without overlap filters
     * @param interFilteredBfp      second fingerprinter, with the inter-substructure overlap filter
     */
    private void assertCounts(String[] smilesArray, int[][] expectedUnfiltered, int[][] expectedInterFiltered,
                              BiosynfoniFingerprinter unfilteredBfp, BiosynfoniFingerprinter interFilteredBfp)
            throws CDKException {
        for (int m = 0; m < smilesArray.length; m++) {
            ICountFingerprint unfilteredCfp = unfilteredBfp.getCountFingerprint(this.smilesParser.parseSmiles(smilesArray[m]));
            ICountFingerprint interFilteredCfp = interFilteredBfp.getCountFingerprint(this.smilesParser.parseSmiles(smilesArray[m]));
            int[] unfilteredCounts = new int[unfilteredBfp.getSize()];
            int[] interFilteredCounts = new int[interFilteredBfp.getSize()];
            for (int i = 0; i < unfilteredBfp.getSize(); i++) {
                unfilteredCounts[i] = unfilteredCfp.getCountForHash(i);
                interFilteredCounts[i] = interFilteredCfp.getCountForHash(i);
            }
            Assertions.assertArrayEquals(expectedUnfiltered[m], unfilteredCounts, smilesArray[m]);
            Assertions.assertArrayEquals(expectedInterFiltered[m], interFilteredCounts, smilesArray[m]);
        }
    }

    /**
     * Verifies that the fingerprint does not depend on whether the hydrogen atoms of a molecule are explicit or
     * implicit, for all four combinations of the overlap filters and for the refined and the original fingerprint. Several patterns use the SMARTS primitives D and h,
     * which count explicit hydrogen atoms as connections and only implicit ones as hydrogen atoms. Without suppressing
     * the explicit hydrogen atoms, the methyl group of toluene (key 20) would not be matched and the oxygen atom of
     * ethanol would be matched as ether (key 29) instead of hydroxyl (key 30).
     */
    @Test
    void testExplicitHydrogens() throws CDKException {
        String[][] pairs = {
                {"CCO", "[H]OC([H])([H])C([H])([H])[H]"},
                {"Cc1ccccc1", "[H]C([H])([H])c1ccccc1"},
                {"NCC(=O)O", "[H]N([H])C([H])([H])C(=O)O[H]"}};
        boolean[] toggles = {false, true};
        for (boolean isOriginalFingerprint : toggles) {
            for (boolean intraSubOverlapToggle : toggles) {
                for (boolean interSubOverlapToggle : toggles) {
                    BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(intraSubOverlapToggle, interSubOverlapToggle,
                            isOriginalFingerprint);
                    for (String[] pair : pairs) {
                        ICountFingerprint implicitCfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(pair[0]));
                        ICountFingerprint explicitCfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(pair[1]));
                        for (int i = 0; i < bfp.getSize(); i++) {
                            Assertions.assertEquals(implicitCfp.getCountForHash(i), explicitCfp.getCountForHash(i),
                                    "key " + i + ", intra " + intraSubOverlapToggle + ", inter " + interSubOverlapToggle
                                            + ", original " + isOriginalFingerprint + ": " + pair[1]);
                        }
                    }
                }
            }
        }
        ICountFingerprint ethanolCfp = this.biosynfoniFingerprinter.getCountFingerprint(this.smilesParser.parseSmiles(pairs[0][1]));
        Assertions.assertEquals(1, ethanolCfp.getCountForHash(30));
        Assertions.assertEquals(0, ethanolCfp.getCountForHash(29));
        ICountFingerprint tolueneCfp = this.biosynfoniFingerprinter.getCountFingerprint(this.smilesParser.parseSmiles(pairs[1][1]));
        Assertions.assertEquals(1, tolueneCfp.getCountForHash(20));
    }

    /**
     * Verifies that the given molecule is not modified, i.e. that neither its explicit hydrogen atoms are removed nor
     * aromaticity flags are set by the SMARTS preparation.
     */
    @Test
    void testInputNotModified() throws CDKException {
        IAtomContainer mol = this.smilesParser.parseSmiles("[H]OC1=CC=CC=C1");
        int atomCount = mol.getAtomCount();
        this.biosynfoniFingerprinter.getCountFingerprint(mol);
        new BiosynfoniFingerprinter(true, true).getBitFingerprint(mol);
        this.biosynfoniFingerprinter.getRawFingerprint(mol);
        Assertions.assertEquals(atomCount, mol.getAtomCount());
        for (int i = 0; i < mol.getAtomCount(); i++) {
            Assertions.assertFalse(mol.getAtom(i).isAromatic(), "atom " + i);
        }
    }

    /**
     * Verifies that nitro groups, nitrate esters, and N-oxides give the same fingerprint in their pentavalent
     * ({@code N(=O)=O}, {@code N=O}) and charge-separated ({@code [N+](=O)[O-]}, {@code [N+][O-]}) form, for all four
     * combinations of the overlap filters, and that their oxygen atoms are not matched as hydroxyl group (key 30),
     * while hydroxylamines (also deprotonated), oximes, hydroxamic acids, and deprotonated hydroxyl groups still are.
     */
    @Test
    void testNitroGroupAndNOxideRepresentations() throws CDKException {
        String[][] pairs = {
                {"CN(=O)=O", "C[N+](=O)[O-]"},
                {"c1ccccc1N(=O)=O", "c1ccccc1[N+](=O)[O-]"},
                {"CON(=O)=O", "CO[N+](=O)[O-]"},
                {"OCc1ccc(cc1)N(=O)=O", "OCc1ccc(cc1)[N+](=O)[O-]"},
                {"CN(C)(C)=O", "C[N+](C)(C)[O-]"},
                {"O=n1ccccc1", "[O-][n+]1ccccc1"}};
        int[] expectedHydroxylCounts = {0, 0, 0, 1, 0, 0};
        boolean[] toggles = {false, true};
        for (boolean intraSubOverlapToggle : toggles) {
            for (boolean interSubOverlapToggle : toggles) {
                BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(intraSubOverlapToggle, interSubOverlapToggle);
                for (String[] pair : pairs) {
                    ICountFingerprint pentavalentCfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(pair[0]));
                    ICountFingerprint chargeSeparatedCfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(pair[1]));
                    for (int i = 0; i < bfp.getSize(); i++) {
                        Assertions.assertEquals(pentavalentCfp.getCountForHash(i), chargeSeparatedCfp.getCountForHash(i),
                                "key " + i + ", intra " + intraSubOverlapToggle + ", inter " + interSubOverlapToggle + ": " + pair[1]);
                    }
                }
            }
        }
        for (int m = 0; m < pairs.length; m++) {
            for (String smiles : pairs[m]) {
                Assertions.assertEquals(expectedHydroxylCounts[m],
                        this.biosynfoniFingerprinter.getCountFingerprint(this.smilesParser.parseSmiles(smiles)).getCountForHash(30), smiles);
            }
        }
        //N-methylhydroxylamine, deprotonated N,N-dimethylhydroxylamine (uncharged nitrogen atom), acetaldoxime,
        // acetohydroxamic acid, and phenolate keep their hydroxyl group
        for (String smiles : new String[]{"CNO", "CN(C)[O-]", "CC=NO", "CC(=O)NO", "[O-]c1ccccc1"}) {
            Assertions.assertEquals(1,
                    this.biosynfoniFingerprinter.getCountFingerprint(this.smilesParser.parseSmiles(smiles)).getCountForHash(30), smiles);
        }
    }

    /**
     * Verifies that the side chain carbon atoms of the phenyl-C3 (key 13) and phenyl-C2 (key 14) units must not be part
     * of any six-membered ring, regardless of the bond orders in this ring and of further fused rings like the lactone
     * ring of aryltetralin lignans, while open side chains are still matched, for all four combinations of the overlap
     * filters.
     */
    @Test
    void testPhenylSideChainRing() throws CDKException {
        //side chain in a six-membered ring: phenylcyclohexane, 4-phenylcyclohexene (non-aromatic double bond in the
        // ring), aryltetralin lignan skeletons with and without lactone ring, and aryldihydronaphthalene lignan
        // skeletons (double bond in the ring) with and without lactone ring
        String[] ringSmilesArray = {"c1ccccc1C1CCCCC1", "c1ccccc1C1CC=CCC1",
                "COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", "OCC1Cc2cc(OC)c(O)cc2C(c2ccc(O)c(OC)c2)C1CO",
                "COc1cc(cc(OC)c1OC)C1=C2C(COC2=O)Cc2cc3OCOc3cc12", "OCC1=Cc2cc(OC)c(O)cc2C(c2ccc(O)c(OC)c2)C1CO"};
        //open side chains: propylbenzene, cinnamyl alcohol, and cannabidiol, whose only phenyl-C3 unit is in the pentyl
        // chain, while its aryl group is attached to a cyclohexene ring
        String[] openSmilesArray = {"CCCc1ccccc1", "OC/C=C/c1ccccc1", "CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1"};
        boolean[] toggles = {false, true};
        for (boolean intraSubOverlapToggle : toggles) {
            for (boolean interSubOverlapToggle : toggles) {
                BiosynfoniFingerprinter bfp = new BiosynfoniFingerprinter(intraSubOverlapToggle, interSubOverlapToggle);
                String filters = ", intra " + intraSubOverlapToggle + ", inter " + interSubOverlapToggle;
                for (String smiles : ringSmilesArray) {
                    ICountFingerprint cfp = bfp.getCountFingerprint(this.smilesParser.parseSmiles(smiles));
                    Assertions.assertEquals(0, cfp.getCountForHash(13), smiles + filters);
                    Assertions.assertEquals(0, cfp.getCountForHash(14), smiles + filters);
                }
                for (String smiles : openSmilesArray) {
                    Assertions.assertEquals(1,
                            bfp.getCountFingerprint(this.smilesParser.parseSmiles(smiles)).getCountForHash(13), smiles + filters);
                }
            }
        }
    }

    /**
     * Verifies that the refined fingerprint is the default and that the original fingerprint is only used if it is
     * explicitly requested in the constructor.
     */
    @Test
    void testOriginalFingerprintSetting() {
        Assertions.assertFalse(BiosynfoniFingerprinter.IS_ORIGINAL_FINGERPRINT_DEFAULT);
        Assertions.assertFalse(new BiosynfoniFingerprinter().isOriginalFingerprint());
        Assertions.assertFalse(new BiosynfoniFingerprinter(true, true).isOriginalFingerprint());
        Assertions.assertFalse(new BiosynfoniFingerprinter(false, false, false).isOriginalFingerprint());
        Assertions.assertTrue(new BiosynfoniFingerprinter(false, false, true).isOriginalFingerprint());
        Assertions.assertEquals(39, new BiosynfoniFingerprinter(false, false, true).getSize());
        int refinedKeys = 0;
        for (BiosynfoniFingerprinter.BiosynfoniKey key : BiosynfoniFingerprinter.BiosynfoniKey.values()) {
            if (!key.getOriginalSmarts().equals(key.getRefinedSmarts())) {
                refinedKeys++;
            }
        }
        //15 refined SMARTS patterns, plus the 8 ring keys that are matched by ring perception
        Assertions.assertEquals(15, refinedKeys);
        Assertions.assertEquals(8, Arrays.stream(BiosynfoniFingerprinter.BiosynfoniKey.values())
                .filter(key -> key.getRingSize() > 0).count());
    }

    /**
     * Verifies every refinement by a molecule whose count of the respective key differs between the original and the
     * refined fingerprint, without overlap filters, and by molecules that are matched by both, see the documentation of
     * the keys in {@link BiosynfoniFingerprinter.BiosynfoniKey}. The counts were cross-checked against the original
     * and the variant pattern set of the Biosynfoni fragmenter of MORTAR, which give the same counts.
     */
    @Test
    void testOriginalVersusRefined() throws CDKException {
        //SMILES, key position, original count, refined count
        Object[][] cases = {
                //NADPH is no NADH (key 1), its three phosphate groups are counted once each (key 21)
                {"C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O", 1, 1, 0},
                {"C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O", 2, 1, 1},
                {"C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O", 21, 12, 3},
                //alanine methyl ester is a standard amino acid (key 3), not a non-standard one (key 4)
                {"COC(=O)[C@@H](C)N", 3, 0, 1},
                {"COC(=O)[C@@H](C)N", 4, 1, 0},
                //ornithine and 2-aminoisobutyric acid are non-standard amino acids (key 4)
                {"NCCC[C@H](N)C(=O)O", 4, 0, 1},
                {"CC(C)(N)C(=O)O", 4, 0, 1},
                //glucopyranose is no open-chain hexose (key 5) or pentose (key 6), its anomeric hydroxyl group is counted
                // (key 30); open-chain glucose is matched by both
                {"OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", 5, 1, 0},
                {"OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", 6, 2, 0},
                {"OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", 30, 4, 5},
                {"OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O", 5, 1, 1},
                //the pyrrole ring of tryptamine is no phenethylamine unit (key 10), tetrahydroisoquinoline and
                // phenethylamine are
                {"NCCc1c[nH]c2ccccc12", 10, 1, 0},
                {"C1Cc2ccccc2CN1", 10, 1, 1},
                {"NCCc1ccccc1", 10, 1, 1},
                //phenyl units (keys 13 to 15) need an aromatic ring and a side chain outside of six-membered rings
                {"CCCC1CCCCC1", 13, 1, 0},
                {"CCCC1CCCCC1", 15, 1, 0},
                {"CCCc1ccccc1", 13, 1, 1},
                {"COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", 13, 2, 0},
                {"COc1cc(cc(OC)c1OC)C1C2C(COC2=O)Cc2cc3OCOc3cc12", 14, 1, 0},
                //the methyl group of a methylmalonyl unit (key 18) is attached to a branching point
                {"CCCC", 18, 2, 0},
                {"CCCc1ccccc1", 18, 1, 0},
                //one count per phosphorus (key 21) and sulfur atom (key 22) bonded to oxygen
                {"COP(=O)(O)O", 21, 4, 1},
                {"COS(=O)(=O)O", 22, 4, 1},
                //nitriles and imines have no terminal nitrogen atom (key 27), amines have
                {"CC#N", 27, 1, 0},
                {"CC=N", 27, 1, 0},
                {"CCN", 27, 1, 1},
                //peroxide oxygen atoms are no ethers (key 29)
                {"COOC", 29, 2, 0},
                //hemiacetal hydroxyl groups are hydroxyl groups (key 30), the oxygen atoms of nitro groups are not
                {"CC(O)OC", 30, 0, 1},
                {"C[N+](=O)[O-]", 30, 1, 0},
                //ring keys (31 to 38) do not match the perimeter of fused and bridged ring systems, but every face of
                // cubane
                {"c1ccc2ccccc2c1", 38, 1, 0},
                {"c1ccc2ccccc2c1", 34, 2, 2},
                {"C1CC2CC2C1", 34, 1, 0},
                {"C1CC2CC2C1", 31, 1, 1},
                {"C12C3C4C1C5C2C3C45", 32, 6, 6},
                {"C12C3C4C1C5C2C3C45", 34, 16, 0},
                {"C12C3C4C1C5C2C3C45", 36, 1, 0},
                {"CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O", 37, 1, 0},
                {"CC12CCC3C(CCC4CC(O)CCC34C)C1CCC2O", 38, 2, 0}};
        BiosynfoniFingerprinter originalBfp = new BiosynfoniFingerprinter(false, false, true);
        for (Object[] testCase : cases) {
            IAtomContainer mol = this.smilesParser.parseSmiles((String) testCase[0]);
            int key = (int) testCase[1];
            String message = "key " + key + ": " + testCase[0];
            Assertions.assertEquals((int) testCase[2], originalBfp.getCountFingerprint(mol).getCountForHash(key),
                    "original, " + message);
            Assertions.assertEquals((int) testCase[3], this.biosynfoniFingerprinter.getCountFingerprint(mol).getCountForHash(key),
                    "refined, " + message);
        }
    }
}
