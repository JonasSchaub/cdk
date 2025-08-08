/*
 * Copyright (c) 2025 Jonas Schaub <jonas.schaub@uni-jena.de>
 *                    Achim Zielesny <achim.zielesny@w-hs.de>
 *                    Christoph Steinbeck <christoph.steinbeck@uni-jena.de>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package org.openscience.cdk.tools;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * JUnit test class for testing the functionalities of the SugarDetectionUtility
 * class.
 * <p>
 *     Identifiers starting with 'CHEMBL' refer to molecules in the
 *     <a href="https://www.ebi.ac.uk/chembl/">ChEMBL database</a>.
 *     <br>Identifiers starting with 'CNP' refer to molecules in the
 *     <a href="https://coconut.naturalproducts.net">COCONUT database</a>.
 *     <br>Identifiers starting with 'CID' refer to molecules in the
 *     <a href="https://pubchem.ncbi.nlm.nih.gov">PubChem database</a>.
 * </p>
 *
 * @author Jonas Schaub
 */
class SugarDetectionUtilityTest {
    /**
     * TODO: remove this before starting a PR
     *
     * A test for quickly inspecting the sugar extraction results of a batch of molecules using CDK Depict.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionTestPrinting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String[] glycosidicNP = new String[] {
                "CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@@H](OP(=O)(O)OP(=O)(O)OCCC(C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)O[C@@H]2CO)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@@H](O)[C@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O)[C@@H]3O)[C@@H]2O)[C@@H]1O",
                "C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                //has issues in this batch processing
                "CC(N)C(=O)NC(CCC(N)=O)C(=O)NOC1OC(O)C(O)C(O)C1O",
                //has issues in this batch processing
                "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O",
                //has issues in this batch processing
                "OC1OC(O)C(O)C1OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1OC1C(O)C(O)C(O)OC(O)C1O",
                "[H]OC1([H])C([H])(OC2=C3C(OC(=O)C4=C3C([H])([H])C([H])([H])C4([H])[H])=C([H])C(=C2[H])C([H])([H])[H])OC([H])(C([H])(O[H])C1([H])O[H])C([H])([H])O[H]",
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4",
                "O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O",
                "O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CN=C(N)C=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C=%11C=CN=C%12NC(NC)CC(C%12%11)CC%10)CO)NCC",
                "O=CCC12OC(OC1=O)C3(C)C(C)CCC3(C)C2",
                "O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O",
                "O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=NCC5)C(C)C)C(O)C1O",
                "O=C1OC(C2=COC=C2)CC3(C)C1CCC4(C)C3C5OC(=O)C4(O)C=C5",
                "O=C1OC2CC3(OC4(O)C(CC5(OC45C(=O)OC)CCCCCCCCCCCCCCCC)C2(O3)C1)CCCCCCCCCCCCCCCC",
                "O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO",
                "O=C(O)CC(OC1OC(CO)C(O)C(O)C1O)(C)CC(=O)OCC=CC2=CC(OC)=C(OC3OC(CO)C(O)C(O)C3O)C(OC)=C2",
                "O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(OC2OC(CO)C(O)C(O)C2O)C=C1",
                "O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(COC(=O)C(CC(=O)O)C(OCC3C(=C)CCC4C(C)(C)CCCC34C)C(=O)OC)CCCC12C)C(=O)OC",
                "O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2OC%10OC(CO)C(O)C(O)C%10O)CO)C(O)C1O",
                "O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC4C5=CCC6C7(C)CCC(O)C(C(=O)OC8OC(CO)C(O)C(O)C8O)(C)C7CCC6(C)C5(C)CC3)C(O)C(OC9OC(CO)C(O)C(O)C9O)C2O)C(OC%10OC(CO)C(O)C(O)C%10O)C(O)C1O",
                "O=C(O)CC(O)(C)CC(=O)OC1C(O)C(OC2C3=C(O)C(=CC=C3OC2C(=C)CO)C(=O)C)OC(CO)C1O",
                "O=C(O)CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O",
                "O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2OC3OCC(O)C(OC4OCC(O)C(O)C4O)C3O)OC5C(OC(=O)C67CCC(C)(C)CC7C8=CCC9C%10(C)CC(O)C(OC%11OC(CO)C(O)C(O)C%11O)C(C(=O)O)(C)C%10CCC9(C)C8(CO)CC6)OC(C)C(OC(=O)C=CC%12=CC(OC)=C(OC)C(OC)=C%12)C5OC%13OC(C)C(O)C(O)C%13O)C1O",
                "O=C(O)C1OC(O)C(O)C(O)C1O",
                "O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(C)CCCC12C)C(=O)O",
                "O=C1C=C(OC=2C=C3OC(C)(CCC4CNC(=O)C4)C(OOCC(O)C(O)C(O)C(O)CO)CC3=CC12)C",
                "O=C1C=C(OC=2C1=CC3=C(OC(C)(C)C(OOCC(O)C(O)C(O)C(O)CO)C3)C2[N+]=4C=C5N=CC=C5C4CC)C",
                "O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4",
                "O=C(O)CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O",
                //has issues in this batch processing
                "O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC",
                "O=CC1(C)C(OC2OC(C(=O)O)C(O)C(OC3OCC(O)C(O)C3O)C2OC4OC(CO)C(O)C(O)C4O)CCC5(C)C6CC=C7C8CC(C)(C)CCC8(C(=O)OC9OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(OC%10OC(CO)C(O)C%10O)C(C)CC)C(C)CC)C(O)C9OC%11OC(C)C(OC%12OCC(O)C(O)C%12O)C(O)C%11O)C(O)CC7(C)C6(C)CCC15",
                "O=C(OCC)CCC1=CC=2C=COC2C=3OCCNCC4(SSCC5(OC(OC31)C(O)C(O)C5O)CO)CCCC4",
                "O=C1OC2CC3(OC(CCC3)CC(O)CC)OC(CC4(O)OC(C)(C)CC4C=CCCCCCC(O)(C)C(O)C(OC5OC(C)C(N(C)C)CC5)C(O)C(C)C(O)C(O)(C=C1)C)C2C",
                "O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C",
                "O=C(O)CC(C)(O)CC(=O)OCc1ccc(cc1)OC1C(CO)OC(OC(C(=O)OCc2ccc(OC3OC(CO)CC(O)C3O)cc2)C(O)(CC(C)C)C(=O)OC2CCc3cc4cc(O)c(C)c(O)c4c(O)c3C2=O)C(O)C1O",
                "OC1(OCCC21OC3(O)CCOC3(O2)C)C",
                "OCC1OC2(OC3C(O)C(OC3(OC2)CO)CO)C(O)C1O",
                "OCC1OC2(OCC3OC(OC2)(CO)C(O)C3O)C(O)C1O",
                //might create issues in this batch processing
                "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
                //might create issues in this batch processing
                "OCC(O)C(O)C(O)C(O)C(O)C1OC(O)C(O)C(O)C1N",
                "OCC1OC2OC3C(O)C(O)C(OC3CO)OC4C(O)C(O)C(OC4CO)OC5C(O)C(O)C(OC5CO)OC6C(O)C(O)C(OC6CO)OC7C(O)C(O)C(OC7CO)OC1C(O)C2O",
                "O(*)C1OC(CN)CCC1N",
                "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
                "O=C1N=C2C(=NC=3C=C(C(=CC3N2CC(O)C(O)C(O)COC4OC(CO)C(O)C(O)C4O)C)C)C(=O)N1",
                "O=C(OCC(O)C(O)C(O)C(O)COC1=C(OC2=CC(O)=CC(OCC(O)C(O)C(O)C(O)CO)=C2C1)C=3C=C(O)C(O)=C(O)C3)C=CC4=CC=C(O)C=C4",
                "O=C(C1(CC2=C3C=CC4C5(CCC(C(CO)(C5CCC4(C3(CC(C2(CC1)CO)O)C)C)C)OC6OC(C(C(C6O)OC7OC(C(C(C7O)O)O)CO)O)C)C)C)OCC(C(C(CO)O)O)O",
                "CC(CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)CCOP(=O)(O)OP(=O)(O)OC1C(C(C(C(O1)CO)OC2C(C(C(C(O2)CO)OC3C(C(C(C(O3)COC4C(C(C(C(O4)COC5C(C(C(C(O5)CO)OC6C(C(C(C(O6)CO)O)O)O)O)O)O)OC7C(C(C(C(O7)CO)OC8C(C(C(C(O8)CO)O)O)O)O)O)O)O)OC9C(C(C(C(O9)CO)O)O)OC1C(C(C(C(O1)CO)O)O)OC1C(C(C(C(O1)CO)O)OC1C(C(C(C(O1)CO)O)OC1C(C(C(C(O1)CO)O)O)OC1C(C(C(C(O1)CO)O)O)O)O)O)O)O)NC(=O)C)O)NC(=O)C",
                "O=CC(O)C(O)C(O)C(O)COC(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O",
                "O=C(O)C=1C(=N)C(=O)CC(O)(C1O)C2OC(COC(=O)C)C(OC(=O)C3NC(=S)SC3C)C(OC4OC(C)C(O)(C(OC(=O)C(C)C)C)C(OC)C4)C2O",
                "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)(C2O)C3OC(CO)C(O)C(O)C3O)C(O)C4OCC(O)C(O)C4O",
                "O=C(O)C=1C(=O)C(O)(CC(=O)C1N)C2OC(COC(=O)C)C(OC(=O)C(N=C(S)SCC(N=C(O)C)C(=O)O)C(SCC(N=C(O)C)C(=O)O)C)C(OC3OC(C)C(O)(C(OC(=O)C(C)CC)C)C(OC)C3)C2O",
                "O=C(O)CC1C(=CC2CCC=CC1C2OC3OC(CO)C4(OC(=C5N=C([CH-]C5=C4)C(C)C)CCCO)C(O)C3(O)CC)C(=O)OC",
                "O=C(OC1C2=C(O)C=3C(=O)C=4C=CC=C(O)C4C(=O)C3C(O)=C2C(OC5OC(C)C(OC6OC(C)C(OC7OC(C(=O)CC7O)C)C(O)C6)C(N(C)C)C5)CC1(O)CC)C",
                //CNP0194094.4
                "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.O=S(=O)(O)O.[Cl-].[Cl-].[K+].[K+]",
                //CNP0194094.5
                "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.O=S(=O)(O)O.[Cl-].[Cl-].[Na+].[Na+]",
                //CNP0491043.0
                "CCC(O)COCC1OC(OC2C(COCC(O)CC)OC(OCC(O)CC)C(OCC(O)CC)C2OCC(O)CC)C(OCC(O)CC)C(OCC(O)CC)C1OCC(O)CC.COCC1OC(OC2C(COC)OC(OC)C(OC)C2OC)C(OC)C(OC)C1OC",
                //CNP0585724.0
                "COC1=CC(C(=O)OC2C3C=COC(C)C3C3(CO)OC23)=CC=C1O.OCC1OC(O)C(O)C(O)C1O",
                //CNP0570079.0
                "CC(=O)OCC1OC(OC2C(COC(C)=O)OC(OC(C)=O)C(OC(C)=O)C2OC(C)=O)C(OC(C)=O)C(OC(C)=O)C1OC(C)=O.CCC(=O)OCC1OC(OC2C(COC(=O)CC)OC(OC(=O)CC)C(OC(=O)CC)C2OC(=O)CC)C(OC(=O)CC)C(OC(=O)CC)C1OC(=O)CC.OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O",
                //CNP0577894.0
                "O=C([O-])C1OC(O)C(O)C(O)C1O.[Na+]",
                //CNP0495754.0
                "CC1OC(O)C(O)C(N)C1O.CNC1=CC=C(C(C)=O)C=C1",
                //CNP0305679.8
                "CCCCCC(=O)O.OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
                //CNP0305679.9
                "CCCCCCCC(=O)O.OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
        };
        for (String smiles : glycosidicNP) {
            List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, true, true);
            if (candidates.size() < 2) {
                continue;
            }
            System.out.println(smiles + " input mol");
            for (IAtomContainer candidate : candidates) {
                if (candidate.isEmpty() && candidates.indexOf(candidate) == 0) {
                    System.out.println("[empty aglycone]");
                    continue;
                }
                System.out.println(smiGen.create(candidate));
            }
        }
    }

    /**
     * Tests the correct circular sugar extraction on some examples.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionTestAsserting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        // Map of input SMILES codes to list of expected aglycone and sugar SMILES codes
        Map<String, List<String>> testCases = new HashMap<>((int) ((60.0/0.75) + 3.0));
        testCases.put(
                "CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@@H](OP(=O)(O)OP(=O)(O)OCCC(C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)O[C@@H]2CO)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@@H](O)[C@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O)[C@@H]3O)[C@@H]2O)[C@@H]1O",
                Arrays.asList(
                        "OP(=O)(O)OP(=O)(O)OCCC(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C",
                        "CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@H](O[C@@H]2CO)O)O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O[C@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]6O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]7O[C@H](CO)[C@@H](O[C@H]8O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]8O)[C@H](O)[C@@H]7O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]9O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]9O[C@H]%10O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]%10O[C@H]%11O[C@H](CO)[C@@H](O)[C@H](O[C@H]%12O[C@H](CO)[C@@H](O)[C@H](O[C@H]%13O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%13O[C@H]%14O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%14O)[C@H]%12O)[C@@H]%11O)[C@@H]3O)[C@@H]1O"
                )
        );
        testCases.put(
                "C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                Arrays.asList(
                        "C=CC1C(C[C@@H]2NCCC3=C2NC4=CC=CC=C34)C(C(=O)O)=CO[C@H]1O",
                        "[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O"
                )
        );
        testCases.put(
                "CC(N)C(=O)NC(CCC(N)=O)C(=O)NOC1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "CC(N)C(=O)NC(CCC(N)=O)C(=O)NO",
                        "C1(OC(O)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1O",
                        "C1(OC(C)C(O)C(O)C1OC2OC(O)C(O)C(O)C2O)O"
                )
        );
        testCases.put(
                "OC1OC(O)C(O)C1OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1OC1C(O)C(O)C(O)OC(O)C1O",
                Arrays.asList(
                        "OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1O",
                        "OC1OC(O)C(O)C1O",
                        "C1(C(O)C(O)C(O)OC(O)C1O)O"
                )
        );
        testCases.put(
                "[H]OC1([H])C([H])(OC2=C3C(OC(=O)C4=C3C([H])([H])C([H])([H])C4([H])[H])=C([H])C(=C2[H])C([H])([H])[H])OC([H])(C([H])(O[H])C1([H])O[H])C([H])([H])O[H]",
                Arrays.asList(
                        "OC1=C2C(OC(=O)C3=C2C([H])([H])C([H])([H])C3([H])[H])=C([H])C(=C1[H])C([H])([H])[H]",
                        "[H]OC1([H])C([H])(OC([H])(C([H])(O[H])C1([H])O[H])C([H])([H])O[H])O"
                )
        );
        testCases.put(
                "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4",
                Collections.singletonList(
                        "O=C(OC1C(OCC2=COC(OC(=O)CC(C)C)C3C2CC(O)C3(O)COC(=O)C)OC(CO)C(O)C1O)C=CC4=CC=C(O)C=C4"
                )
        );
        testCases.put(
                "O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O",
                Collections.singletonList(
                        "O=P(O)(O)OCC1OC(OP(=O)(O)O)C(O)C1O"
                )
        );
        testCases.put(
                "O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CN=C(N)C=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C=%11C=CN=C%12NC(NC)CC(C%12%11)CC%10)CO)NCC",
                Collections.singletonList(
                        "O=C1OC2C(CCO)CCC3(C=C4C=CCC5C(C=CC(C45)C23)CCCC(C)(CC6=CN=C(N)C=C6)CC=7C=CC=C8C(=O)C9(OC19C(=O)C87)CC(=C(C)CC%10C=%11C=CN=C%12NC(NC)CC(C%12%11)CC%10)CO)NCC"
                )
        );
        testCases.put(
                "O=CCC12OC(OC1=O)C3(C)C(C)CCC3(C)C2",
                Collections.singletonList(
                        "O=CCC12OC(OC1=O)C3(C)C(C)CCC3(C)C2"
                )
        );
        testCases.put(
                "O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O",
                Collections.singletonList(
                        "O=C(O)C12OC(OC3=CC=4OCC5C6=C(OC5C4C(=C3)C7=CC=CC(O)=C7)C(OC)=C(OC)C=C6CNC)(CO)C(O)C(O)(NCC1NC)C2O"
                )
        );
        testCases.put(
                "O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=NCC5)C(C)C)C(O)C1O",
                Collections.singletonList(
                        "O=C(O)C1OC(OC=2C=CC=3C(=O)[C-](C=[O+]C3C2)C4=CC=C(O)C=C4)C(O)(CNCC(CC=5C=NCC5)C(C)C)C(O)C1O"
                )
        );
        testCases.put(
                "O=C1OC(C2=COC=C2)CC3(C)C1CCC4(C)C3C5OC(=O)C4(O)C=C5",
                Collections.singletonList(
                        "O=C1OC(C2=COC=C2)CC3(C)C1CCC4(C)C3C5OC(=O)C4(O)C=C5"
                )
        );
        testCases.put(
                "O=C1OC2CC3(OC4(O)C(CC5(OC45C(=O)OC)CCCCCCCCCCCCCCCC)C2(O3)C1)CCCCCCCCCCCCCCCC",
                Collections.singletonList(
                        "O=C1OC2CC3(OC4(O)C(CC5(OC45C(=O)OC)CCCCCCCCCCCCCCCC)C2(O3)C1)CCCCCCCCCCCCCCCC"
                )
        );
        testCases.put(
                "O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO",
                Collections.singletonList(
                        "O=C1C2=CC=CC3=C2CN1CC(=O)C4=C(O)C5=C6OC7OC(COC(C=CC6=C(OC)C8=C5C=9C(=CC%10CCCC%10C49)CC8)C%11=CNC=%12C=CC(=CC%12%11)CNC)C(O)C(OC#CC3)C7(O)CO"
                )
        );
        testCases.put(
                "O=C(O)CC(OC1OC(CO)C(O)C(O)C1O)(C)CC(=O)OCC=CC2=CC(OC)=C(OC3OC(CO)C(O)C(O)C3O)C(OC)=C2",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC=CC1=CC(OC)=C(O)C(OC)=C1",
                        "C1(OC(CO)C(O)C(O)C1O)O",
                        "C1(OC(CO)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(OC2OC(CO)C(O)C(O)C2O)C=C1",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(O)C=C1",
                        "C1(OC(CO)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(COC(=O)C(CC(=O)O)C(OCC3C(=C)CCC4C(C)(C)CCCC34C)C(=O)OC)CCCC12C)C(=O)OC",
                Collections.singletonList(
                        "O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(COC(=O)C(CC(=O)O)C(OCC3C(=C)CCC4C(C)(C)CCCC34C)C(=O)OC)CCCC12C)C(=O)OC"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2OC%10OC(CO)C(O)C(O)C%10O)CO)C(O)C1O",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2O)CO)C(O)C1O",
                        "C1(OC(CO)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC4C5=CCC6C7(C)CCC(O)C(C(=O)OC8OC(CO)C(O)C(O)C8O)(C)C7CCC6(C)C5(C)CC3)C(O)C(OC9OC(CO)C(O)C(O)C9O)C2O)C(OC%10OC(CO)C(O)C(O)C%10O)C(O)C1O",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC3C5=CCC6C7(C)CCC(O)C(C(=O)O)(C)C7CCC6(C)C5(C)CC4)C(O)C(O)C2O)C(O)C(O)C1O",
                        "C1(OC(CO)C(O)C(O)C1O)O",
                        "C1(OC(CO)C(O)C(O)C1O)O",
                        "C1(OC(CO)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OC1C(O)C(OC2C3=C(O)C(=CC=C3OC2C(=C)CO)C(=O)C)OC(CO)C1O",
                Collections.singletonList(
                        "O=C(O)CC(O)(C)CC(=O)OC1C(O)C(OC2C3=C(O)C(=CC=C3OC2C(=C)CO)C(=O)C)OC(CO)C1O"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O",
                Collections.singletonList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1OC(C=2C(O)=CC(O)=C3C(=O)C=C(OC32)C=4C=CC(O)=C(O)C4)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2OC3OCC(O)C(OC4OCC(O)C(O)C4O)C3O)OC5C(OC(=O)C67CCC(C)(C)CC7C8=CCC9C%10(C)CC(O)C(OC%11OC(CO)C(O)C(O)C%11O)C(C(=O)O)(C)C%10CCC9(C)C8(CO)CC6)OC(C)C(OC(=O)C=CC%12=CC(OC)=C(OC)C(OC)=C%12)C5OC%13OC(C)C(O)C(O)C%13O)C1O",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2O)OC3C(OC(=O)C45CCC(C)(C)CC4C6=CCC7C8(C)CC(O)C(O)C(C(=O)O)(C)C8CCC7(C)C6(CO)CC5)OC(C)C(OC(=O)C=CC9=CC(OC)=C(OC)C(OC)=C9)C3O)C1O",
                        "C1(OCC(O)C(OC2OCC(O)C(O)C2O)C1O)O",
                        "C1(OC(CO)C(O)C(O)C1O)O",
                        "C1(OC(C)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(O)C1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "",
                        "O=C(O)C1OC(O)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(C)CCCC12C)C(=O)O",
                Collections.singletonList(
                        "O=C(O)CC(C(=O)O)C(OCC1C(=C)CCC2C(C)(C)CCCC12C)C(=O)O"
                )
        );
        testCases.put(
                "O=C1C=C(OC=2C=C3OC(C)(CCC4CNC(=O)C4)C(OOCC(O)C(O)C(O)C(O)CO)CC3=CC12)C",
                Collections.singletonList(
                        "O=C1C=C(OC=2C=C3OC(C)(CCC4CNC(=O)C4)C(OOCC(O)C(O)C(O)C(O)CO)CC3=CC12)C"
                )
        );
        testCases.put(
                "O=C1C=C(OC=2C1=CC3=C(OC(C)(C)C(OOCC(O)C(O)C(O)C(O)CO)C3)C2[N+]=4C=C5N=CC=C5C4CC)C",
                Collections.singletonList(
                        "O=C1C=C(OC=2C1=CC3=C(OC(C)(C)C(OOCC(O)C(O)C(O)C(O)CO)C3)C2[N+]=4C=C5N=CC=C5C4CC)C"
                )
        );
        testCases.put(
                "O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4",
                Arrays.asList(
                        "O=C1C=C(OC2=CC(OC=O)=C(O)C(O)=C12)C=3C=CC(O)=CC3",
                        "C1OC(O)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC",
                Collections.singletonList(
                        "O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC"
                )
        );
        testCases.put(
                "O=CC1(C)C(OC2OC(C(=O)O)C(O)C(OC3OCC(O)C(O)C3O)C2OC4OC(CO)C(O)C(O)C4O)CCC5(C)C6CC=C7C8CC(C)(C)CCC8(C(=O)OC9OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(OC%10OC(CO)C(O)C%10O)C(C)CC)C(C)CC)C(O)C9OC%11OC(C)C(OC%12OCC(O)C(O)C%12O)C(O)C%11O)C(O)CC7(C)C6(C)CCC15",
                Arrays.asList(
                        "O=CC1(C)C(O)CCC2(C)C3CC=C4C5CC(C)(C)CCC5(C(=O)OC6OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(O)C(C)CC)C(C)CC)C(O)C6O)C(O)CC4(C)C3(C)CCC12",
                        "C1(OC(C(=O)O)C(O)C(OC2OCC(O)C(O)C2O)C1OC3OC(CO)C(O)C(O)C3O)O",
                        "C1(OC(CO)C(O)C1O)O",
                        "C1(OC(C)C(OC2OCC(O)C(O)C2O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C",
                Arrays.asList(
                        "",
                        "O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C"
                )
        );
        testCases.put(
                "O=C(O)CC(C)(O)CC(=O)OCc1ccc(cc1)OC1C(CO)OC(OC(C(=O)OCc2ccc(OC3OC(CO)CC(O)C3O)cc2)C(O)(CC(C)C)C(=O)OC2CCc3cc4cc(O)c(C)c(O)c4c(O)c3C2=O)C(O)C1O",
                Arrays.asList(
                        "O=C(O)CC(C)(O)CC(=O)OCC1=CC=C(C=C1)OC2C(CO)OC(OC(C(=O)OCC3=CC=C(O)C=C3)C(O)(CC(C)C)C(=O)OC4CCC5=CC6=CC(O)=C(C)C(O)=C6C(O)=C5C4=O)C(O)C2O",
                        "C1(OC(CO)CC(O)C1O)O"
                )
        );
        testCases.put(
                "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
                Arrays.asList(
                        "OCC(O)C(O)C(O)CO",
                        "C1OC(CO)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "OCC(O)C(O)C(O)C(O)C(O)C1OC(O)C(O)C(O)C1N",
                Arrays.asList(
                        "OCC(O)C(O)C(O)C(O)CO",
                        "C1OC(O)C(O)C(O)C1N"
                )
        );
        testCases.put(
                "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
                Arrays.asList(
                        "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(O)=CC3C2C4=CC=C(O)C(O)=C4)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
                        "C1(OC(CO)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C1N=C2C(=NC=3C=C(C(=CC3N2CC(O)C(O)C(O)COC4OC(CO)C(O)C(O)C4O)C)C)C(=O)N1",
                Arrays.asList(
                        "O=C1N=C2C(=NC=3C=C(C(=CC3N2CC(O)C(O)C(O)CO)C)C)C(=O)N1",
                        "C1(OC(CO)C(O)C(O)C1O)O"
                )
        );
        testCases.put(
                "O=C(C1(CC2=C3C=CC4C5(CCC(C(CO)(C5CCC4(C3(CC(C2(CC1)CO)O)C)C)C)OC6OC(C(C(C6O)OC7OC(C(C(C7O)O)O)CO)O)C)C)C)OCC(C(C(CO)O)O)O",
                Arrays.asList(
                        "O=C(C1(CC2=C3C=CC4C5(CCC(C(CO)(C5CCC4(C3(CC(C2(CC1)CO)O)C)C)C)O)C)C)OCC(C(C(CO)O)O)O",
                        "C1(OC(C(C(C1O)OC2OC(C(C(C2O)O)O)CO)O)C)O"
                )
        );
        testCases.put(
                "CC(CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)CCOP(=O)(O)OP(=O)(O)OC1C(C(C(C(O1)CO)OC2C(C(C(C(O2)CO)OC3C(C(C(C(O3)COC4C(C(C(C(O4)COC5C(C(C(C(O5)CO)OC6C(C(C(C(O6)CO)O)O)O)O)O)O)OC7C(C(C(C(O7)CO)OC8C(C(C(C(O8)CO)O)O)O)O)O)O)O)OC9C(C(C(C(O9)CO)O)O)OC1C(C(C(C(O1)CO)O)O)OC1C(C(C(C(O1)CO)O)OC1C(C(C(C(O1)CO)O)OC1C(C(C(C(O1)CO)O)O)OC1C(C(C(C(O1)CO)O)O)O)O)O)O)O)NC(=O)C)O)NC(=O)C",
                Arrays.asList(
                        "CC(CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)CCOP(=O)(O)OP(=O)(O)O",
                        "C1(C(C(C(C(O1)CO)OC2C(C(C(C(O2)CO)OC3C(C(C(C(O3)COC4C(C(C(C(O4)COC5C(C(C(C(O5)CO)OC6C(C(C(C(O6)CO)O)O)O)O)O)O)OC7C(C(C(C(O7)CO)OC8C(C(C(C(O8)CO)O)O)O)O)O)O)O)OC9C(C(C(C(O9)CO)O)O)OC%10C(C(C(C(O%10)CO)O)O)OC%11C(C(C(C(O%11)CO)O)OC%12C(C(C(C(O%12)CO)O)OC%13C(C(C(C(O%13)CO)O)O)OC%14C(C(C(C(O%14)CO)O)O)O)O)O)O)O)NC(=O)C)O)NC(=O)C)O"
                )
        );
        testCases.put(
                "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)(C2O)C3OC(CO)C(O)C(O)C3O)C(O)C4OCC(O)C(O)C4O",
                Arrays.asList(
                        "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)C2O)CO",
                        "C1OC(CO)C(O)C(O)C1O",
                        "C1OCC(O)C(O)C1O"
                )
        );
        // CNP0194094.4
        testCases.put(
                "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.O=S(=O)(O)O.[Cl-].[Cl-].[K+].[K+]",
                Arrays.asList(
                        "O=S(=O)(O)O.[Cl-].[Cl-].[K+].[K+]",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"
                )
        );
        // CNP0194094.5
        testCases.put(
                "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.O=S(=O)(O)O.[Cl-].[Cl-].[Na+].[Na+]",
                Arrays.asList(
                        "O=S(=O)(O)O.[Cl-].[Cl-].[Na+].[Na+]",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"
                )
        );
        // CNP0491043.0
        testCases.put(
                "CCC(O)COCC1OC(OC2C(COCC(O)CC)OC(OCC(O)CC)C(OCC(O)CC)C2OCC(O)CC)C(OCC(O)CC)C(OCC(O)CC)C1OCC(O)CC.COCC1OC(OC2C(COC)OC(OC)C(OC)C2OC)C(OC)C(OC)C1OC",
                Arrays.asList(
                        "CCC(O)COCC1OC(OC2C(COCC(O)CC)OC(OCC(O)CC)C(OCC(O)CC)C2OCC(O)CC)C(OCC(O)CC)C(OCC(O)CC)C1OCC(O)CC",
                        "COCC1OC(OC2C(COC)OC(OC)C(OC)C2OC)C(OC)C(OC)C1OC"
                )
        );
        // CNP0585724.0
        testCases.put(
                "COC1=CC(C(=O)OC2C3C=COC(C)C3C3(CO)OC23)=CC=C1O.OCC1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "COC1=CC(C(=O)OC2C3C=COC(C)C3C4(CO)OC24)=CC=C1O",
                        "OCC1OC(O)C(O)C(O)C1O"
                )
        );
        // CNP0570079.0
        testCases.put(
                "CC(=O)OCC1OC(OC2C(COC(C)=O)OC(OC(C)=O)C(OC(C)=O)C2OC(C)=O)C(OC(C)=O)C(OC(C)=O)C1OC(C)=O.CCC(=O)OCC1OC(OC2C(COC(=O)CC)OC(OC(=O)CC)C(OC(=O)CC)C2OC(=O)CC)C(OC(=O)CC)C(OC(=O)CC)C1OC(=O)CC.OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O",
                Arrays.asList(
                        "CC(=O)OCC1OC(OC2C(COC(C)=O)OC(OC(C)=O)C(OC(C)=O)C2OC(C)=O)C(OC(C)=O)C(OC(C)=O)C1OC(C)=O.CCC(=O)OCC1OC(OC2C(COC(=O)CC)OC(OC(=O)CC)C(OC(=O)CC)C2OC(=O)CC)C(OC(=O)CC)C(OC(=O)CC)C1OC(=O)CC",
                        "OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O"
                )
        );
        // CNP0577894.0
        testCases.put(
                "O=C([O-])C1OC(O)C(O)C(O)C1O.[Na+]",
                Arrays.asList(
                        "[Na+]",
                        "O=C([O-])C1OC(O)C(O)C(O)C1O"
                )
        );
        // CNP0495754.0
        testCases.put(
                "CC1OC(O)C(O)C(N)C1O.CNC1=CC=C(C(C)=O)C=C1",
                Arrays.asList(
                        "CNC1=CC=C(C(C)=O)C=C1",
                        "CC1OC(O)C(O)C(N)C1O"
                )
        );
        // CNP0305679.8
        testCases.put(
                "CCCCCC(=O)O.OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
                Arrays.asList(
                        "CCCCCC(=O)O",
                        "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
                )
        );
        // CNP0305679.9
        testCases.put(
                "CCCCCCCC(=O)O.OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
                Arrays.asList(
                        "CCCCCCCC(=O)O",
                        "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
                )
        );
        // Process each test case
        for (Map.Entry<String, List<String>> entry : testCases.entrySet()) {
            String inputSmiles = entry.getKey();
            List<String> expectedSmilesList = entry.getValue();
            List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(inputSmiles), true, false, false);
            List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
            Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
        }
    }

    /**
     * Tests the correct circular sugar extraction with R saturation on some examples.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionWithRTestAsserting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        // Map of input SMILES codes to list of expected aglycone and sugar SMILES codes
        Map<String, List<String>> testCases = new HashMap<>((int) ((60.0/0.75) + 3.0));
        // Add these test cases to the testCases map in sugarExtractionWithRTestAsserting()
        testCases.put(
                "CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@@H](OP(=O)(O)OP(=O)(O)OCCC(C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)O[C@@H]2CO)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@@H](O)[C@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O)[C@@H]3O)[C@@H]2O)[C@@H]1O",
                Arrays.asList(
                        "O(P(=O)(O)OP(=O)(O)OCCC(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)*",
                        "CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@H](O[C@@H]2CO)O*)O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O[C@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]6O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]7O[C@H](CO)[C@@H](O[C@H]8O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]8O)[C@H](O)[C@@H]7O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]9O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]9O[C@H]%10O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]%10O[C@H]%11O[C@H](CO)[C@@H](O)[C@H](O[C@H]%12O[C@H](CO)[C@@H](O)[C@H](O[C@H]%13O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%13O[C@H]%14O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%14O)[C@H]%12O)[C@@H]%11O)[C@@H]3O)[C@@H]1O"
                )
        );
        testCases.put(
                "C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
                Arrays.asList(
                        "C=CC1C(C[C@@H]2NCCC3=C2NC4=CC=CC=C34)C(C(=O)O)=CO[C@H]1O*",
                        "[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O*"
                )
        );
        testCases.put(
                "CC(N)C(=O)NC(CCC(N)=O)C(=O)NOC1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "CC(N)C(=O)NC(CCC(N)=O)C(=O)NO*",
                        "C1(OC(O)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1O*",
                        "C1(OC(C)C(O)C(O)C1OC2OC(O)C(O)C(O)C2O)O*"
                )
        );
        testCases.put(
                "OC1OC(O)C(O)C1OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1OC1C(O)C(O)C(O)OC(O)C1O",
                Arrays.asList(
                        "O(C1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1O*)*",
                        "OC1OC(O)C(O)C1O*",
                        "C1(C(O)C(O)C(O)OC(O)C1O)O*"
                )
        );
        testCases.put(
                "[H]OC1([H])C([H])(OC2=C3C(OC(=O)C4=C3C([H])([H])C([H])([H])C4([H])[H])=C([H])C(=C2[H])C([H])([H])[H])OC([H])(C([H])(O[H])C1([H])O[H])C([H])([H])O[H]",
                Arrays.asList(
                        "O(C1=C2C(OC(=O)C3=C2C([H])([H])C([H])([H])C3([H])[H])=C([H])C(=C1[H])C([H])([H])[H])*",
                        "[H]OC1([H])C([H])(OC([H])(C([H])(O[H])C1([H])O[H])C([H])([H])O[H])O*"
                )
        );
        testCases.put(
                "O=C(O)CC(OC1OC(CO)C(O)C(O)C1O)(C)CC(=O)OCC=CC2=CC(OC)=C(OC3OC(CO)C(O)C(O)C3O)C(OC)=C2",
                Arrays.asList(
                        "O=C(O)CC(O*)(C)CC(=O)OCC=CC1=CC(OC)=C(O*)C(OC)=C1",
                        "C1(OC(CO)C(O)C(O)C1O)O*",
                        "C1(OC(CO)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(OC2OC(CO)C(O)C(O)C2O)C=C1",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1=CC=C(O*)C=C1",
                        "C1(OC(CO)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2OC%10OC(CO)C(O)C(O)C%10O)CO)C(O)C1O",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OC1COC(OC2C(O)C(OC(OC3C(O)C(O)C(OC4CC5CCC6C(CCC7(C)C6CC8OC9(OCC(C)CC9)C(C)C87)C5(C)CC4O)OC3CO)C2O*)CO)C(O)C1O",
                        "C1(OC(CO)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC4C5=CCC6C7(C)CCC(O)C(C(=O)OC8OC(CO)C(O)C(O)C8O)(C)C7CCC6(C)C5(C)CC3)C(O)C(OC9OC(CO)C(O)C(O)C9O)C2O)C(OC%10OC(CO)C(O)C(O)C%10O)C(O)C1O",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1OC(OCC2OC(OC(=O)C34CCC(C)(C)CC3C5=CCC6C7(C)CCC(O)C(C(=O)O*)(C)C7CCC6(C)C5(C)CC4)C(O)C(O*)C2O)C(O*)C(O)C1O",
                        "C1(OC(CO)C(O)C(O)C1O)O*",
                        "C1(OC(CO)C(O)C(O)C1O)O*",
                        "C1(OC(CO)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2OC3OCC(O)C(OC4OCC(O)C(O)C4O)C3O)OC5C(OC(=O)C67CCC(C)(C)CC7C8=CCC9C%10(C)CC(O)C(OC%11OC(CO)C(O)C(O)C%11O)C(C(=O)O)(C)C%10CCC9(C)C8(CO)CC6)OC(C)C(OC(=O)C=CC%12=CC(OC)=C(OC)C(OC)=C%12)C5OC%13OC(C)C(O)C(O)C%13O)C1O",
                Arrays.asList(
                        "O=C(O)CC(O)(C)CC(=O)OCC1(O)COC(OC2C(O)C(OC(C)C2O*)OC3C(OC(=O)C45CCC(C)(C)CC4C6=CCC7C8(C)CC(O)C(O*)C(C(=O)O)(C)C8CCC7(C)C6(CO)CC5)OC(C)C(OC(=O)C=CC9=CC(OC)=C(OC)C(OC)=C9)C3O*)C1O",
                        "C1(OCC(O)C(OC2OCC(O)C(O)C2O)C1O)O*",
                        "C1(OC(CO)C(O)C(O)C1O)O*",
                        "C1(OC(C)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(O)C1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "",
                        "O=C(O)C1OC(O)C(O)C(O)C1O")
        );
        testCases.put(
                "O=C1C=C(OC2=CC(OC(=O)C3OC(O)C(O)C(O)C3O)=C(O)C(O)=C12)C=4C=CC(O)=CC4",
                Arrays.asList(
                        "O=C1C=C(OC2=CC(OC(=O)*)=C(O)C(O)=C12)C=3C=CC(O)=CC3",
                        "C1(OC(O)C(O)C(O)C1O)*"
                )
        );
        testCases.put(
                "O=CC1(C)C(OC2OC(C(=O)O)C(O)C(OC3OCC(O)C(O)C3O)C2OC4OC(CO)C(O)C(O)C4O)CCC5(C)C6CC=C7C8CC(C)(C)CCC8(C(=O)OC9OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(OC%10OC(CO)C(O)C%10O)C(C)CC)C(C)CC)C(O)C9OC%11OC(C)C(OC%12OCC(O)C(O)C%12O)C(O)C%11O)C(O)CC7(C)C6(C)CCC15",
                Arrays.asList(
                        "O=CC1(C)C(O*)CCC2(C)C3CC=C4C5CC(C)(C)CCC5(C(=O)OC6OC(C)C(OC(=O)CC(O)CC(OC(=O)CC(O)CC(O*)C(C)CC)C(C)CC)C(O)C6O*)C(O)CC4(C)C3(C)CCC12",
                        "C1(OC(C(=O)O)C(O)C(OC2OCC(O)C(O)C2O)C1OC3OC(CO)C(O)C(O)C3O)O*",
                        "C1(OC(CO)C(O)C1O)O*",
                        "C1(OC(C)C(OC2OCC(O)C(O)C2O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C",
                Arrays.asList(
                        "",
                        "O=C(NC1C(O)OC(CO)C(O)C1OC2OC(CO)C(OC)C(O)C2OC3OC(C)C(O)C(O)C3OC)C")
        );
        testCases.put(
                "O=C(O)CC(C)(O)CC(=O)OCc1ccc(cc1)OC1C(CO)OC(OC(C(=O)OCc2ccc(OC3OC(CO)CC(O)C3O)cc2)C(O)(CC(C)C)C(=O)OC2CCc3cc4cc(O)c(C)c(O)c4c(O)c3C2=O)C(O)C1O",
                Arrays.asList(
                        "O=C(O)CC(C)(O)CC(=O)OCC1=CC=C(C=C1)OC2C(CO)OC(OC(C(=O)OCC3=CC=C(O*)C=C3)C(O)(CC(C)C)C(=O)OC4CCC5=CC6=CC(O)=C(C)C(O)=C6C(O)=C5C4=O)C(O)C2O",
                        "C1(OC(CO)CC(O)C1O)O*"
                )
        );
        testCases.put(
                "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O",
                Arrays.asList(
                        "OCC(O)C(O)C(O)C(O)*",
                        "C1(OC(CO)C(O)C(O)C1O)*"
                )
        );
        testCases.put(
                "OCC(O)C(O)C(O)C(O)C(O)C1OC(O)C(O)C(O)C1N",
                Arrays.asList(
                        "OCC(O)C(O)C(O)C(O)C(O)*",
                        "C1(OC(O)C(O)C(O)C1N)*"
                )
        );
        testCases.put(
                "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(OC4OC(CO)C(O)C(O)C4O)=CC3C2C5=CC=C(O)C(O)=C5)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
                Arrays.asList(
                        "O=C(O)C1=CC(O)C(O)C(OC(=O)C2C(=CC=3C=C(O)C(O*)=CC3C2C4=CC=C(O)C(O)=C4)C(=O)OCC(O)C(O)C(O)C(O)C(O)CO)C1",
                        "C1(OC(CO)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C1N=C2C(=NC=3C=C(C(=CC3N2CC(O)C(O)C(O)COC4OC(CO)C(O)C(O)C4O)C)C)C(=O)N1",
                Arrays.asList(
                        "O=C1N=C2C(=NC=3C=C(C(=CC3N2CC(O)C(O)C(O)CO*)C)C)C(=O)N1",
                        "C1(OC(CO)C(O)C(O)C1O)O*"
                )
        );
        testCases.put(
                "O=C(C1(CC2=C3C=CC4C5(CCC(C(CO)(C5CCC4(C3(CC(C2(CC1)CO)O)C)C)C)OC6OC(C(C(C6O)OC7OC(C(C(C7O)O)O)CO)O)C)C)C)OCC(C(C(CO)O)O)O",
                Arrays.asList(
                        "O=C(C1(CC2=C3C=CC4C5(CCC(C(CO)(C5CCC4(C3(CC(C2(CC1)CO)O)C)C)C)O*)C)C)OCC(C(C(CO)O)O)O",
                        "C1(OC(C(C(C1O)OC2OC(C(C(C2O)O)O)CO)O)C)O*"
                )
        );
        testCases.put(
                "CC(CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)CCOP(=O)(O)OP(=O)(O)OC1C(C(C(C(O1)CO)OC2C(C(C(C(O2)CO)OC3C(C(C(C(O3)COC4C(C(C(C(O4)COC5C(C(C(C(O5)CO)OC6C(C(C(C(O6)CO)O)O)O)O)O)O)OC7C(C(C(C(O7)CO)OC8C(C(C(C(O8)CO)O)O)O)O)O)O)O)OC9C(C(C(C(O9)CO)O)O)OC1C(C(C(C(O1)CO)O)O)OC1C(C(C(C(O1)CO)O)OC1C(C(C(C(O1)CO)O)OC1C(C(C(C(O1)CO)O)O)OC1C(C(C(C(O1)CO)O)O)O)O)O)O)O)NC(=O)C)O)NC(=O)C",
                Arrays.asList(
                        "CC(CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C)CCOP(=O)(O)OP(=O)(O)O*",
                        "C1(C(C(C(C(O1)CO)OC2C(C(C(C(O2)CO)OC3C(C(C(C(O3)COC4C(C(C(C(O4)COC5C(C(C(C(O5)CO)OC6C(C(C(C(O6)CO)O)O)O)O)O)O)OC7C(C(C(C(O7)CO)OC8C(C(C(C(O8)CO)O)O)O)O)O)O)O)OC9C(C(C(C(O9)CO)O)O)OC%10C(C(C(C(O%10)CO)O)O)OC%11C(C(C(C(O%11)CO)O)OC%12C(C(C(C(O%12)CO)O)OC%13C(C(C(C(O%13)CO)O)O)OC%14C(C(C(C(O%14)CO)O)O)O)O)O)O)O)NC(=O)C)O)NC(=O)C)O*"
                )
        );
        testCases.put(
                "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)(C2O)C3OC(CO)C(O)C(O)C3O)C(O)C4OCC(O)C(O)C4O",
                Arrays.asList(
                        "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)(C2O)*)C(O)*",
                        "C1(OC(CO)C(O)C(O)C1O)*",
                        "C1(OCC(O)C(O)C1O)*"
                )
        );
        testCases.put(
                "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.O=S(=O)(O)O.[Cl-].[Cl-].[K+].[K+]",
                Arrays.asList(
                        "O=S(=O)(O)O.[Cl-].[Cl-].[K+].[K+]",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"
                )
        );
        testCases.put(
                "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O.O=S(=O)(O)O.[Cl-].[Cl-].[Na+].[Na+]",
                Arrays.asList(
                        "O=S(=O)(O)O.[Cl-].[Cl-].[Na+].[Na+]",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
                        "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O"
                )
        );
        testCases.put(
                "CCC(O)COCC1OC(OC2C(COCC(O)CC)OC(OCC(O)CC)C(OCC(O)CC)C2OCC(O)CC)C(OCC(O)CC)C(OCC(O)CC)C1OCC(O)CC.COCC1OC(OC2C(COC)OC(OC)C(OC)C2OC)C(OC)C(OC)C1OC",
                Arrays.asList(
                        "CCC(O)COCC1OC(OC2C(COCC(O)CC)OC(OCC(O)CC)C(OCC(O)CC)C2OCC(O)CC)C(OCC(O)CC)C(OCC(O)CC)C1OCC(O)CC",
                        "COCC1OC(OC2C(COC)OC(OC)C(OC)C2OC)C(OC)C(OC)C1OC"
                )
        );
        testCases.put(
                "COC1=CC(C(=O)OC2C3C=COC(C)C3C3(CO)OC23)=CC=C1O.OCC1OC(O)C(O)C(O)C1O",
                Arrays.asList(
                        "COC1=CC(C(=O)OC2C3C=COC(C)C3C4(CO)OC24)=CC=C1O",
                        "OCC1OC(O)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "CC(=O)OCC1OC(OC2C(COC(C)=O)OC(OC(C)=O)C(OC(C)=O)C2OC(C)=O)C(OC(C)=O)C(OC(C)=O)C1OC(C)=O.CCC(=O)OCC1OC(OC2C(COC(=O)CC)OC(OC(=O)CC)C(OC(=O)CC)C2OC(=O)CC)C(OC(=O)CC)C(OC(=O)CC)C1OC(=O)CC.OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O",
                Arrays.asList(
                        "CC(=O)OCC1OC(OC2C(COC(C)=O)OC(OC(C)=O)C(OC(C)=O)C2OC(C)=O)C(OC(C)=O)C(OC(C)=O)C1OC(C)=O.CCC(=O)OCC1OC(OC2C(COC(=O)CC)OC(OC(=O)CC)C(OC(=O)CC)C2OC(=O)CC)C(OC(=O)CC)C(OC(=O)CC)C1OC(=O)CC",
                        "OCC1OC(OC2C(CO)OC(O)C(O)C2O)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "O=C([O-])C1OC(O)C(O)C(O)C1O.[Na+]",
                Arrays.asList(
                        "[Na+]",
                        "O=C([O-])C1OC(O)C(O)C(O)C1O"
                )
        );
        testCases.put(
                "CC1OC(O)C(O)C(N)C1O.CNC1=CC=C(C(C)=O)C=C1",
                Arrays.asList(
                        "CNC1=CC=C(C(C)=O)C=C1",
                        "CC1OC(O)C(O)C(N)C1O"
                )
        );
        testCases.put(
                "CCCCCC(=O)O.OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
                Arrays.asList(
                        "CCCCCC(=O)O",
                        "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
                )
        );
        testCases.put(
                "CCCCCCCC(=O)O.OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
                Arrays.asList(
                        "CCCCCCCC(=O)O",
                        "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
                )
        );
        // Process each test case
        for (Map.Entry<String, List<String>> entry : testCases.entrySet()) {
            String inputSmiles = entry.getKey();
            List<String> expectedSmilesList = entry.getValue();
            List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(inputSmiles), true, false, true);
            List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
            Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
        }
    }

    /**
     * Test for correct circular sugar extraction from strictosidinic acid.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionStrictosidinicAcidTest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0225072.2
        IAtomContainer mol = smiPar.parseSmiles("C=CC1C(C[C@@H]2NCCC3=C2NC2=CC=CC=C32)C(C(=O)O)=CO[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //aglycone
        Assertions.assertEquals("C=CC1C(C[C@@H]2NCCC3=C2NC4=CC=CC=C34)C(C(=O)O)=CO[C@H]1O", smiGen.create(candidates.get(0)));
        //beta-D-glucose
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct circular sugar extraction from natural product CNP0580945.1.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionCNP0580945_1Test() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0580945.1
        IAtomContainer mol = smiPar.parseSmiles("CC(=O)N[C@H]1[C@H](O[C@@H]2C3CO[C@H](O3)[C@H](NC(C)=O)[C@H]2O[C@H](C)C(=O)N[C@@H](C)C(=O)N[C@H](CCC(=O)N[C@@H](CCC[C@@H](N)C(=O)O)C(=O)N[C@H](C)C(=O)O)C(=O)O)O[C@H](CO)[C@@H](O)[C@@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //aglycone
        Assertions.assertEquals("O[C@@H]1C2CO[C@H](O2)[C@H](NC(C)=O)[C@H]1O[C@H](C)C(=O)N[C@@H](C)C(=O)N[C@H](CCC(=O)N[C@@H](CCC[C@@H](N)C(=O)O)C(=O)N[C@H](C)C(=O)O)C(=O)O", smiGen.create(candidates.get(0)));
        //N-Acetyl-beta-D-glucosamine
        Assertions.assertEquals("CC(=O)N[C@H]1[C@@H](O[C@H](CO)[C@@H](O)[C@@H]1O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from Fusacandin B (CNP0295326.2).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionFusacandin_BTest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0295326.2
        IAtomContainer mol = smiPar.parseSmiles("CCCCC/C=C/C=C/C(O)C/C=C/C=C/C(=O)O[C@@H]1[C@@H](O)[C@H](C2=C(O)C=C(O)C=C2CO)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //aglycone
        Assertions.assertEquals("CCCCCC=CC=CC(O)CC=CC=CC(=O)O[C@@H]1[C@@H](O)[C@H](C2=C(O)C=C(O)C=C2CO)O[C@H](CO)[C@H]1O", smiGen.create(candidates.get(0)));
        //2-O-beta-D-galactopyranosyl-beta-D-galactopyranose
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)O", smiGen.create(candidates.get(1)));
        sdu.setRemoveOnlyTerminalSugarsSetting(false);
        candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //disconnected aglycone
        Assertions.assertEquals("CCCCCC=CC=CC(O)CC=CC=CC(=O)O.C1=C(O)C=C(O)C=C1CO", smiGen.create(candidates.get(0)));
        //triple sugar
        Assertions.assertEquals("[C@H]1([C@@H](O)CO[C@H](CO)[C@H]1O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from CID 131999265.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionCID131999265Test() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CID	131999265
        IAtomContainer mol = smiPar.parseSmiles("CCCCCCCCCCCCOCC(COCCOCC(COCCCCCCCCCCCC)O[C@H]1[C@H](C([C@@H](C(O1)CO)O)O)O)O[C@H]2[C@H](C([C@@H](C(O2)CO)O)O)O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(3, candidates.size());
        //aglycone
        Assertions.assertEquals("CCCCCCCCCCCCOCC(COCCOCC(COCCCCCCCCCCCC)O)O", smiGen.create(candidates.get(0)));
        //sugar 1, CID	58265196 (underdefined hexopyranose)
        Assertions.assertEquals("[C@@H]1([C@H](C([C@@H](C(O1)CO)O)O)O)O", smiGen.create(candidates.get(1)));
        //sugar 2 (the same)
        Assertions.assertEquals("[C@@H]1([C@H](C([C@@H](C(O1)CO)O)O)O)O", smiGen.create(candidates.get(2)));
    }

    /**
     * Test for correct sugar extraction from CNP0381981.2.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionCNP0381981_2Test() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0381981.2
        IAtomContainer mol = smiPar.parseSmiles("CC1=C(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)C=CC2=C1OC(=O)C1=C2CCCC1");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //aglycone
        Assertions.assertEquals("CC1=C(O)C=CC2=C1OC(=O)C3=C2CCCC3", smiGen.create(candidates.get(0)));
        //beta-D-glucose
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from CNP0151033.2 containing a non-terminal sugar.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionCNP0151033_2NonTerminalTest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0151033.2
        IAtomContainer mol = smiPar.parseSmiles("CC(=O)OC[C@]1(O)[C@H]2[C@H](OC(=O)CC(C)C)OC=C(CO[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3OC(=O)/C=C/C3=CC=C(O)C=C3)[C@H]2C[C@@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        sdu.setRemoveOnlyTerminalSugarsSetting(false);
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //disconnected(!) aglycone
        Assertions.assertEquals("CC(=O)OC[C@]1(O)[C@H]2[C@H](OC(=O)CC(C)C)OC=C(CO)[C@H]2C[C@@H]1O.OC(=O)C=CC1=CC=C(O)C=C1", smiGen.create(candidates.get(0)));
        //beta-D-glucose
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from CNP0125332.1, ribose bisphosphate, which should be recognized as a sugar completely.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionRiboseBisPhosphateTest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0125332.1
        IAtomContainer mol = smiPar.parseSmiles("O=P(O)(O)OC[C@H]1O[C@H](OP(=O)(O)O)[C@H](O)[C@@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        sdu.setPreservationModeThresholdSetting(7);
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //empty aglycone
        Assertions.assertEquals("", smiGen.create(candidates.get(0)));
        //Ribose 1,5-bisphosphate
        Assertions.assertEquals("O=P(O)(O)OC[C@H]1O[C@H](OP(=O)(O)O)[C@H](O)[C@@H]1O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from CNP0104886.1, which contains a sugar acid.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionCNP0104886_1Test() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0104886.1
        IAtomContainer mol = smiPar.parseSmiles("CNC(=N)NC[C@H](CNC[C@]1(O)[C@H](OC2=CC=C3C(=O)C(C4=CC=C(O)C=C4)=COC3=C2)O[C@H](C(=O)O)[C@@H](O)[C@@H]1O)C(C)C");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        sdu.setRemoveOnlyTerminalSugarsSetting(false);
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //disconnected aglycone
        Assertions.assertEquals("CNC(=N)NC[C@H](CNC)C(C)C.OC1=CC=C2C(=O)C(C3=CC=C(O)C=C3)=COC2=C1", smiGen.create(candidates.get(0)));
        //sugar acid
        Assertions.assertEquals("[C@H]1(O)[C@@H](O[C@H](C(=O)O)[C@@H](O)[C@@H]1O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from Tangshenoside I (CNP0218440.3), which contains two beta-D-glucose moieties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionTangshenoside_ITest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0218440.3
        IAtomContainer mol = smiPar.parseSmiles("COC1=CC(/C=C/COC(=O)C[C@](C)(CC(=O)O)O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)=CC(OC)=C1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(3, candidates.size());
        //aglycone
        Assertions.assertEquals("COC1=CC(C=CCOC(=O)C[C@](C)(CC(=O)O)O)=CC(OC)=C1O", smiGen.create(candidates.get(0)));
        //beta-D-glucose
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O", smiGen.create(candidates.get(1)));
        //beta-D-glucose
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)O", smiGen.create(candidates.get(2)));
    }

    /**
     * Test for correct sugar extraction from Gitonin (CNP0209335.4), which contains a complex sugar structure.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionGitoninTest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0209335.4
        IAtomContainer mol = smiPar.parseSmiles("C[C@@H]1CC[C@@]2(OC1)O[C@H]1C[C@H]3[C@@H]4CC[C@H]5C[C@@H](O[C@@H]6O[C@H](CO)[C@H](O[C@@H]7O[C@H](CO)[C@@H](O)[C@H](O[C@@H]8OC[C@@H](O)[C@H](O)[C@H]8O)[C@H]7O[C@@H]7O[C@H](CO)[C@H](O)[C@H](O)[C@H]7O)[C@H](O)[C@H]6O)[C@H](O)C[C@]5(C)[C@H]4CC[C@]3(C)[C@H]1[C@@H]2C");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //aglycone
        Assertions.assertEquals("C[C@@H]1CC[C@@]2(OC1)O[C@H]3C[C@H]4[C@@H]5CC[C@H]6C[C@@H](O)[C@H](O)C[C@]6(C)[C@H]5CC[C@]4(C)[C@H]3[C@@H]2C", smiGen.create(candidates.get(0)));
        //connected sugars
        Assertions.assertEquals("[C@H]1(O[C@H](CO)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O[C@@H]3OC[C@@H](O)[C@H](O)[C@H]3O)[C@H]2O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]1O)O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from alpha-l-mannuronic acid (CNP0171089.22), which should be recognized as a sugar.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionAlphaMannuronicAcidTest() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0171089.22
        IAtomContainer mol = smiPar.parseSmiles("O=C(O)[C@@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //empty aglycone
        Assertions.assertEquals("", smiGen.create(candidates.get(0)));
        //alpha-l-mannuronic acid
        Assertions.assertEquals("O=C(O)[C@@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from glycan G00008, which contains a complex structure with multiple sugars.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionglycanG00008Test() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0083402.1
        IAtomContainer mol = smiPar.parseSmiles("CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@@H](OP(=O)(O)OP(=O)(O)OCCC(C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)O[C@@H]2CO)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO[C@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@@H](O)[C@H](O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]3O[C@H]3O[C@H](CO)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O)[C@@H]3O)[C@@H]2O)[C@@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        //aglycone
        Assertions.assertEquals("OP(=O)(O)OP(=O)(O)OCCC(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C", smiGen.create(candidates.get(0)));
        //sugars
        Assertions.assertEquals("CC(=O)N[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](NC(C)=O)[C@H](O[C@@H]2CO)O)O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O[C@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]6O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]7O[C@H](CO)[C@@H](O[C@H]8O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]8O)[C@H](O)[C@@H]7O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]9O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]9O[C@H]%10O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]%10O[C@H]%11O[C@H](CO)[C@@H](O)[C@H](O[C@H]%12O[C@H](CO)[C@@H](O)[C@H](O[C@H]%13O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%13O[C@H]%14O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%14O)[C@H]%12O)[C@@H]%11O)[C@@H]3O)[C@@H]1O", smiGen.create(candidates.get(1)));
    }

    /**
     * Test for correct sugar extraction from a glycosidic natural product with an ester bond between aglycone and sugar.
     * TODO: The extraction result should be improved in a future version, to produce an alcohol and a carboxy acid moiety on
     * either side of the former ester bond.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionTestEster() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        //CNP0214620.1
        IAtomContainer mol = smiPar.parseSmiles("O=C(OC1=CC2=C(C(O)=C1O)C(=O)C=C(C1=CC=C(O)C(O)=C1)O2)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O");
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(mol, true, false, false);
        Assertions.assertEquals(2, candidates.size());
        // aglycone
        Assertions.assertEquals("O=COC1=CC2=C(C(O)=C1O)C(=O)C=C(C3=CC=C(O)C(O)=C3)O2", smiGen.create(candidates.get(0)));
        //sugar
        Assertions.assertEquals("C1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", smiGen.create(candidates.get(1)));
    }

    /**
     * The following tests had the purpose to single out molecules that had issues in batch processing in the past.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest1() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String glycosidicNP = "OCC(O)C(O)C(O)C(O)C(O)C1OC(O)C(O)C(O)C1N";
        List<String> expectedSmilesList = Arrays.asList(
                "OCC(O)C(O)C(O)C(O)CO",
                "C1OC(O)C(O)C(O)C1N"
        );
        IAtomContainer molecule = smiPar.parseSmiles(glycosidicNP);
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(molecule, true, false, false);
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * See above.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest2() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, false, false);
        List<String> expectedSmilesList = Arrays.asList(
                "O=C(O)CC(O)(C(=O)O)C(C(=O)O)CCCCCCCCCCCCCC"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * See above.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest3() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "OCC(O)C(O)C(O)C(O)C1OC(CO)C(O)C(O)C1O";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, false, false);
        List<String> expectedSmilesList = Arrays.asList(
                "OCC(O)C(O)C(O)CO",
                "C1OC(CO)C(O)C(O)C1O"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * See above.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest4() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)(C2O)C3OC(CO)C(O)C(O)C3O)C(O)C4OCC(O)C(O)C4O";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, false, false);
        List<String> expectedSmilesList = Arrays.asList(
                "O=C(C=CC1=CC=C(O)C=C1)C=2C(=O)C(C(=O)C(O)C2O)CO",
                "C1OC(CO)C(O)C(O)C1O",
                //TODO this sugar should include the C-OH connecting it to the aglycone
                "C1OCC(O)C(O)C1O"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * See above.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest5() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "CC(N)C(=O)NC(CCC(N)=O)C(=O)NOC1OC(O)C(O)C(O)C1O";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, false, false);
        List<String> expectedSmilesList = Arrays.asList(
                "CC(N)C(=O)NC(CCC(N)=O)C(=O)NO",
                "C1(OC(O)C(O)C(O)C1O)O"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * See above.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest6() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "OC1OC(O)C(O)C1OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1OC1C(O)C(O)C(O)OC(O)C1O";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, false, false);
        List<String> expectedSmilesList = Arrays.asList(
                "OC1C(OCCCCCCCCCCCCCCCCC)OC(OCCCCCCCCCCC)C(O)C1O",
                "OC1OC(O)C(O)C1O",
                "C1(C(O)C(O)C(O)OC(O)C1O)O"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * See above.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest7() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, false, false);
        List<String> expectedSmilesList = Arrays.asList(
                "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1O",
                "C1(OC(C)C(O)C(O)C1OC2OC(O)C(O)C(O)C2O)O"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * This test is for a specific case where the sugar extraction should yield an empty aglycone but did not in the past.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest8() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0336316.1
        String smiles = "NC(=O)N[C@@H](C=O)[C@@H](O)[C@H](O)[C@H](O)CO";
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, true, false);
        List<String> expectedSmilesList = Arrays.asList(
                "",
                "NC(=O)N[C@@H](C=O)[C@@H](O)[C@H](O)[C@H](O)CO"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    /**
     * This test is for a specific case where the sugar extraction should yield an empty aglycone but did not in the past.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void sugarExtractionIndividualTest9() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Stereo);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0595604.0
        String smiles = "CC(=O)NC(C=O)C(O)C(OC1OC(CO)C(OC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(OC3OC(CO)C(O)C(O)C3O)C2O)C(O)C1NC(C)=O)C(O)CO";
        List<IAtomContainer> candidates =sdu.copyAndExtractAglyconeAndSugars(smiPar.parseSmiles(smiles), true, true, false);
        List<String> expectedSmilesList = Arrays.asList(
                "",
                "CC(=O)NC(C=O)C(O)C(OC1OC(CO)C(OC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(OC4OC(CO)C(O)C(O)C4O)C2O)C(O)C1NC(C)=O)C(O)CO"
        );
        List<String> generatedSmilesList = this.generateSmilesList(candidates, smiGen);
        Assertions.assertLinesMatch(expectedSmilesList, generatedSmilesList);
    }

    //TODO: remove this before starting the PR
    /**
     * Test for processing the COCONUT database. This method reads the COCONUT SDF file and processes each molecule to
     * extract aglycone and sugar candidates.
     *
     * @throws Exception if anything goes wrong
     */
    //@Test
    void testCOCONUT() throws Exception {
        IteratingSDFReader reader = new IteratingSDFReader(
                new FileReader("C:\\Users\\jonas\\Research\\Projects\\Project_COCONUT_Curation\\COCONUT_versions\\coconut_sdf_2d_lite-03-2025.sdf"),
                SilentChemObjectBuilder.getInstance());
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        int exceptionsCounter = 0;
        List<String> coconutIds = new ArrayList<>(30);
        while (reader.hasNext()) {
            IAtomContainer molecule = reader.next();
            try {
                List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(molecule, true, true, true);
//            if (candidates.size() > 1) {
//                System.out.println("COCONUT ID: " + molecule.getProperty("COCONUT_ID"));
//                System.out.println("Aglycone: " + candidates.get(0));
//                System.out.println("Sugar: " + candidates.get(1));
//            }
            } catch (Exception e) {
                System.err.println("Error processing molecule with COCONUT ID: " + molecule.getProperty("identifier"));
                exceptionsCounter++;
                coconutIds.add(molecule.getProperty("identifier"));
                e.printStackTrace();
            }
        }
        System.out.println("Total exceptions encountered: " + exceptionsCounter);
        for (String id : coconutIds) {
            System.out.println(id);
        }
    }

    /**
     * Tests the postprocessing split of ester groups connecting sugar moieties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testEsterSplitting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0138295
        String smiles = "O=CC(O)C(O)C(O)C(O)COC(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O";
        IAtomContainer molecule = smiPar.parseSmiles(smiles);
        sdu.splitEsters(molecule, false);
        Assertions.assertEquals("O=CC(O)C(O)C(O)C(O)COC(O)(C(O)CO)C(O)C(O)C=O.O=C(O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C",
                smiGen.create(molecule));
        molecule = smiPar.parseSmiles(smiles);
        sdu.splitEsters(molecule, true);
        Assertions.assertEquals("*OC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C.*OCC(O)C(O)(OCC(O)C(O)C(O)C(O)C=O)C(O)C(O)C=O",
                smiGen.create(molecule));
    }

    /**
     * Tests the postprocessing split of ether groups cross-linking sugar moieties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testEtherCrosslinkingSplitting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0138295
        String smiles = "O=CC(O)C(O)C(O)C(O)COC(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O";
        IAtomContainer molecule = smiPar.parseSmiles(smiles);
        sdu.splitEthersCrosslinking(molecule, false);
        Assertions.assertEquals("O=CC(O)C(O)C(O)C(O)CO.O=CC(O)C(O)C(O)C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C",
                smiGen.create(molecule));
        molecule = smiPar.parseSmiles(smiles);
        sdu.splitEthersCrosslinking(molecule, true);
        Assertions.assertEquals("*OCC(O)C(O)C(O)C(O)C=O.*C(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O",
                smiGen.create(molecule));
    }

    /**
     * Tests the postprocessing split of ether groups connecting sugar moieties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testEtherSplitting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0138295
        String smiles = "O=CC(O)C(O)C(O)C(O)COC(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O";
        IAtomContainer molecule = smiPar.parseSmiles(smiles);
        sdu.splitEthers(molecule, false);
        //note that this shows why the ether splitting should be done last because it also matches esters and crosslinking ehthers that are treated differently
        Assertions.assertEquals("O=CC(O)C(O)C(O)C(O)CO.O=CC(O)C(O)C(O)(O)C(O)CO.O=C(O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C",
                smiGen.create(molecule));
        molecule = smiPar.parseSmiles(smiles);
        sdu.splitEthers(molecule, true);
        Assertions.assertEquals("*OC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C.*OCC(O)C(O)C(O)C(O)C=O.*OCC(O)C(O)(O*)C(O)C(O)C=O",
                smiGen.create(molecule));
    }

    /**
     * Tests the postprocessing split of peroxide groups connecting sugar moieties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testPeroxideSplitting() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0206763.1
        String smiles = "COC1=C2C[C@@H](C(C)(C)O)OC2=CC2=C1C(=O)C=C(COOC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO)O2";
        IAtomContainer molecule = smiPar.parseSmiles(smiles);
        sdu.splitPeroxides(molecule, false);
        //note that this shows why the ether splitting should be done last because it also matches esters and crosslinking ehthers that are treated differently
        Assertions.assertEquals("O=C1C=C(OC=2C=C3OC(CC3=C(OC)C12)C(O)(C)C)CO.OCC(O)C(O)C(O)C(O)CO",
                smiGen.create(molecule));
        molecule = smiPar.parseSmiles(smiles);
        sdu.splitPeroxides(molecule, true);
        Assertions.assertEquals("*OCC=1OC=2C=C3OC(CC3=C(OC)C2C(=O)C1)C(O)(C)C.*OCC(O)C(O)C(O)C(O)CO",
                smiGen.create(molecule));
    }

    /**
     * Test for splitting O-glycosidic bonds in a molecule, using the splitOGlycosidicBonds method used for postprocessing
     * extracted circular sugar moieties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testSplitOGlycosidicBonds() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0595604.0
        String smiles = "CC(=O)NC(C=O)C(O)C(OC1OC(CO)C(OC2OC(COC3OC(CO)C(O)C(O)C3O)C(O)C(OC3OC(CO)C(O)C(O)C3O)C2O)C(O)C1NC(C)=O)C(O)CO";
        IAtomContainer molecule = smiPar.parseSmiles(smiles);
        sdu.splitOGlycosidicBonds(molecule, false);
        Assertions.assertEquals("O=CC(NC(=O)C)C(O)C(O)C(O)CO.O=C(NC1C(O)OC(CO)C(O)C1O)C.OCC1OC(O)C(O)C(O)C1O.OCC1OC(O)C(O)C(O)C1O.OCC1OC(O)C(O)C(O)C1O",
                smiGen.create(molecule));
        molecule = smiPar.parseSmiles(smiles);
        sdu.splitOGlycosidicBonds(molecule, true);
        Assertions.assertEquals("*OCC1OC(O*)C(O)C(O*)C1O.*OC1OC(CO)C(O)C(O)C1O.*OC1OC(CO)C(O)C(O)C1O.*OC1OC(CO)C(O*)C(O)C1NC(=O)C.*OC(C(O)CO)C(O)C(C=O)NC(=O)C",
                smiGen.create(molecule));
    }

    /**
     * Test for splitting ether, ester, and peroxide bonds in a molecule, using the postprocessing method for extracted
     * sugar moeties.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testSplitEtherEsterPeroxideBondsPostprocessing() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SmilesGenerator smiGen = new SmilesGenerator(SmiFlavor.Canonical);
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0138295
        String smiles = "O=CC(O)C(O)C(O)C(O)COC(O)(C(O)COC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C)C(O)C(O)C=O";
        IAtomContainer molecule = smiPar.parseSmiles(smiles);
        sdu.splitEtherEsterAndPeroxideBondsPostProcessing(molecule, false);
        Assertions.assertEquals("O=CC(O)C(O)C(O)C(O)CO.O=CC(O)C(O)C(O)C(O)CO.O=C(O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C",
                smiGen.create(molecule));
        molecule = smiPar.parseSmiles(smiles);
        sdu.splitEtherEsterAndPeroxideBondsPostProcessing(molecule, true);
        Assertions.assertEquals("*OC(=O)C(O)C(O)C(O)C(O)COC1=CC=CC=2C(=O)C3=CC(=CC(O)=C3C(=O)C12)C.*OCC(O)C(O)C(O)C(O)C=O.*OCC(O)C(*)(O)C(O)C(O)C=O",
                smiGen.create(molecule));
    }

    /**
     * Tests the retrieval of atom indices for aglycone and sugar moieties from a molecule.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testRetrievalOfAtomIndices() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        String smiles = "CCCCCC=CC=CC(O)CC=CC=CC(=O)OC1C(O)C(C2=C(O)C=C(O)C=C2CO)OC(CO)C1OC1OC(C)C(O)C(O)C1OC1OC(O)C(O)C(O)C1O";
        IAtomContainer mol = smiPar.parseSmiles(smiles);
        Map<IAtom, IAtom> inputAtomToAglyconeAtomMap = new HashMap<>((int) ((mol.getAtomCount() / 0.75f) + 2), 0.75f);
        Map<IAtom, IAtom> inputAtomToSugarAtomMap = new HashMap<>((int) ((mol.getAtomCount() / 0.75f) + 2), 0.75f);
        List<IAtomContainer> candidates = sdu.copyAndExtractAglyconeAndSugars(
                mol,
                true,
                false,
                false,
                true,
                inputAtomToAglyconeAtomMap,
                new HashMap<>((int) ((mol.getAtomCount() / 0.75f) + 2), 0.75f),
                inputAtomToSugarAtomMap,
                new HashMap<>((int) ((mol.getAtomCount() / 0.75f) + 2), 0.75f));
        int[] aglyconeAtomIndices = sdu.getAtomIndicesOfGroup(mol, candidates.get(0), inputAtomToAglyconeAtomMap);
        Assertions.assertArrayEquals(new int[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38}, aglyconeAtomIndices);
        Assertions.assertEquals(3, candidates.size());
        int[] sugarOneAtomIndices = sdu.getAtomIndicesOfGroup(mol, candidates.get(1), inputAtomToSugarAtomMap);
        int[] sugarTwoAtomIndices = sdu.getAtomIndicesOfGroup(mol, candidates.get(2), inputAtomToSugarAtomMap);
        Assertions.assertArrayEquals(new int[] {38, 39, 40, 41, 42, 43, 44, 45, 46, 47}, sugarOneAtomIndices);
        Assertions.assertArrayEquals(new int[] {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58}, sugarTwoAtomIndices);
    }

    /**
     * Tests the retrieval of group indices for all atoms in a molecule, which is useful for visualisation.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testGetGroupIndicesForAllAtoms() throws Exception {
        SmilesParser smiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        SugarDetectionUtility sdu = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
        //CNP0295326.4
        String smiles = "CCCCC/C=C/C=C/[C@@H](O)C/C=C/C=C/C(=O)OC1C(O)[C@H](C2=C(O)C=C(O)C=C2CO)O[C@H](CO)[C@H]1O[C@@H]1OC(CO)[C@H](O)[C@H](O)C1O[C@@H]1OC(CO)[C@H](O)[C@H](O)C1O";
        IAtomContainer mol = smiPar.parseSmiles(smiles);
        Map<IAtom, IAtom> inputAtomToAglyconeAtomMap = new HashMap<>((int) ((mol.getAtomCount() / 0.75f) + 2), 0.75f);
        Map<IAtom, IAtom> inputAtomToSugarAtomMap = new HashMap<>((int) ((mol.getAtomCount() / 0.75f) + 2), 0.75f);
        List<IAtomContainer> aglyconeAndSugarsList = sdu.copyAndExtractAglyconeAndSugars(
                mol,
                true,
                false,
                false,
                true,
                inputAtomToAglyconeAtomMap,
                new HashMap<>((int) ((mol.getBondCount() / 0.75f) + 2), 0.75f),
                inputAtomToSugarAtomMap,
                new HashMap<>((int) ((mol.getBondCount() / 0.75f) + 2), 0.75f));
        int[] groupIndices = sdu.getGroupIndicesForAllAtoms(mol, aglyconeAndSugarsList, inputAtomToAglyconeAtomMap, inputAtomToSugarAtomMap);
        for (IAtom atom : mol.atoms()) {
            atom.setMapIdx(groupIndices[atom.getIndex()] + 1);
        }
        String smi = new SmilesGenerator(SmiFlavor.Isomeric | SmiFlavor.AtomAtomMap).create(mol);
        Assertions.assertEquals("[CH3:1][CH2:1][CH2:1][CH2:1][CH2:1]/[CH:1]=[CH:1]/[CH:1]=[CH:1]/[C@@H:1]([OH:1])[CH2:1]/[CH:1]=[CH:1]/[CH:1]=[CH:1]/[C:1](=[O:1])[O:1][CH:1]1[CH:1]([OH:1])[C@H:1]([C:1]2=[C:1]([OH:1])[CH:1]=[C:1]([OH:1])[CH:1]=[C:1]2[CH2:1][OH:1])[O:1][C@H:1]([CH2:1][OH:1])[C@H:1]1[O:2][C@@H:2]3[O:2][CH:2]([CH2:2][OH:2])[C@H:2]([OH:2])[C@H:2]([OH:2])[CH:2]3[O:3][C@@H:3]4[O:3][CH:3]([CH2:3][OH:3])[C@H:3]([OH:3])[C@H:3]([OH:3])[CH:3]4[OH:3]",
                smi);
    }

    /**
     * Helper method to convert a list of IAtomContainer to a list of SMILES strings.
     *
     * @param candidates List of IAtomContainer molecules to convert
     * @param smiGen SMILES generator to use for conversion
     * @return List of SMILES strings
     * @throws CDKException if SMILES generation fails
     */
    protected List<String> generateSmilesList(List<IAtomContainer> candidates, SmilesGenerator smiGen) throws CDKException {
        List<String> result = new ArrayList<>(candidates.size());
        for (IAtomContainer mol : candidates) {
            result.add(smiGen.create(mol));
        }
        return result;
    }
}