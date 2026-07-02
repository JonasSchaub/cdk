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

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;

/**
 * Implementation of the Biosynfoni molecular fingerprint (<a href="https://doi.org/10.1186/s13321-025-01081-6">
 * Biosynfoni paper (Journal for Chemoinformatics, 2025)</a>), a lightweight,
 * biosynthesis-informed, and interpretable fingerprint specifically designed
 * for natural products and natural product-inspired molecules.
 * Look at the <a href="https://github.com/lucinamay/biosynfoni">
 * Original Python implementation</a> for more context
 * <p>
 * The fingerprint consists of a predefined set of SMARTS patterns describing
 * biosynthetically relevant structural motifs. Each SMARTS pattern corresponds
 * to one fingerprint position (feature key). There are 39 SMARTS used.
 * </p>
 *
 * <p>
 * Two fingerprint types are supported:
 * </p>
 * <ul>
 *   <li><b>Bit fingerprint</b>: a bit is set if at least one match for the
 *       corresponding SMARTS pattern exists.</li>
 *   <li><b>Count fingerprint</b>: stores the number of accepted matches for
 *       each SMARTS pattern. This the default fingerprint representation described in the reference publication. </li>
 * </ul>
 *
 *
 * <p>
 * Like the original Biosynfoni implementation, this version can optionally
 * apply overlap filtering:
 * <ul>
 *   <li><b>Intra-pattern filtering</b>: prevents multiple matches of the same
 *       SMARTS pattern from reusing atoms.</li>
 *   <li><b>Inter-pattern filtering</b>: prevents matches of later SMARTS
 *       patterns from reusing atoms already assigned to earlier patterns.</li>
 * </ul>
 * </p>
 *
 * <p>
 * Prior to matching, molecules are preprocessed by assigning atom types,
 * adding implicit hydrogens, detecting ring membership and aromaticity.
 * The supplied molecule may therefore be modified during fingerprint generation.
 * </p>
 *
 *   <p><b>Known differences to the reference implementation:</b></p>
 *   <ul>
 *     <li>
 *       The SMARTS pattern used to identify carbon atoms not inside a 6 ring  was replaced
 *       from {@code [#6;!$([r6])]} to {@code [#6!$(*1*****1)]}. The original
 *       pattern excludes only carbon atoms in six-membered rings, whereas the
 *       replacement excludes carbon atoms that are part of any six-membered
 *       ring. This change was introduced because the SMARTS matching behaviour
 *       in CDK did not reproduce the intended semantics of the original
 *       Biosynfoni implementation. The affected fingerprint keys represent phenyl-derived substructures
 *       originating from the shikimate pathway (fingerprint keys 13, 14, and 15)
 *       See the <a href="https://github.com/cdk/cdk/issues/1292#issue-4641968710">  Git Issue </a>
 *       for more information.
 *     </li>
 *     <li>
 *       Matching can always differ due to the use of different toolkits and therefore the use of different aromaticity models for example:
 *       The aromaticity assignment of positively charged sulfur atoms may differ
 *       between toolkits, resulting in different SMARTS matching behaviour.
 *     </li>
 *   </ul>
 *
 *   <li> Count fingerprints generated using inter-substructure overlap filtering
 *   differ from the reference implementation(<a href="https://github.com/lucinamay/biosynfoni">
 *   Original Python implementation</a>) because this implementation
 *   resolves overlapping matches using canonical atom numbering rather than
 *   the input atom order. This removes path dependency and ensures
 *   deterministic fingerprints regardless of atom indexing.  </li>
 * </ul>
 * </p>
 *
 * @author Marlon Raffelt (MaRa1778)
 */
public class BiosynfoniFingerprinter extends AbstractFingerprinter implements IFingerprinter {

    /**
     * Default Biosynfoni feature definitions.
     * <p>
     * Each enum constant defines a SMARTS pattern and its corresponding label.
     * Together, these feature definitions form the default Biosynfoni fingerprint.
     * </p>
     */
    public enum DefaultBiosynfoniKey {
        //Cofactors
        /**
         * Coenzyme A (CoA).
         * SMARTS Matches the CoA backbone including the pantetheine chain,
         * diphosphate linker and adenosine moiety.
         */
        CO_COA("co_coa", "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O"),
        /**
         * Nicotinamide adenine dinucleotide (NADH).
         * SMARTS Matches the complete reduced NADH cofactor.
         */
        CO_NADH("co_nadh", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O"),// 1: co_nadh
        /**
         * Nicotinamide adenine dinucleotide phosphate (NADPH).
         * SMARTS matches the complete reduced NADPH cofactor.
         */
        CO_NADPH("co_nadph", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O"),
        //Aminoacids
        /**
         * Standard proteinogenic amino acids.
         * SMARTS matches the common amino acid backbone together with one of the
         * twenty canonical side chains.
         */
        ALL_STD_AMINOS("allstnd_aminos", "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]"),
        //#todo
        NON_STD_AMINOS("nonstnd_aminos", "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]"),
        //Sugar deviates
        /**
         * Open-chain hexose.
         * SMARTS matches a linear six-carbon sugar containing a hydroxyl group
         * on every carbon atom.
         */
        S_OPENPYR_C6O6("s_openpyr_C6O6", "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])"),
        /**
         * Open-chain pentose.
         * SMARTS matches a linear five-carbon sugar containing hydroxyl groups.
         */
        S_OPENFUR_C5O5("s_openfur_C5O5", "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])"),
        /**
         * Pyranose sugar.
         * SMARTS matches a six-membered sugar ring containing one oxygen atom.
         */
        S_PYRANOSE_C5O4("s_pyranose_C5O4", "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1"),
        /**
         * Furanose sugar.
         * SMARTS matches a five-membered sugar ring containing one oxygen atom.
         */
        S_FURANOSE_C4O3("s_furanose_C4O3", "C~1~[#8]~C~C(~[#8])~C(~[#8])~1"),
        /**
         * Indole ring system.
         * SMARTS matches an indole scaffold with a two-carbon side chain terminating
         * in a nitrogen atom.
         */
        D_INDOLE("d_indoleC2N_12", "c1cccc2c1c(~[#6]~[#6]~[#7])cn2"),
        /**
         * Phenethylamine motif.
         * SMARTS matches a benzene ring connected to a two-carbon aliphatic chain
         * ending in a nitrogen atom.
         */
        D_PHENYL_C2N("d_phenylC2N_9", "c1ccccc1[#6][#6][#7]"),
        /**
         * Six-membered nitrogen-containing ring.
         * SMARTS matches a six-membered ring consisting of five carbon atoms and one
         * nitrogen atom, regardless of aromaticity.
         */
        D_C5N("d_c5n_6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1"),
        /**
         * Five-membered nitrogen-containing ring.
         * SMARTS matches a five-membered ring consisting of four carbon atoms and one
         * nitrogen atom, regardless of aromaticity.
         */
        D_C4N("d_c4n_5", "[#6]~1~[#6]~[#6]~[#6]~[#7]~1"),

        //Phenyls from Shikimate pathway

        /**
         * Phenylpropyl motif.
         * SMARTS matches a phenyl ring attached to a three-carbon side chain.
         */
        D_PHENYL_C3("d_phenylC3_9_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]~[#6!$(*1*****1)]"),
        /**
         * Phenethyl motif.
         * SMARTS matches a phenyl ring attached to a two-carbon side chain.
         */
        D_PHENYL_C2("d_phenylC2_8_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]"),
        /**
         * Benzyl motif.
         * SMARTS matches a phenyl ring attached to a single carbon atom.
         */
        D_PHENYL_C1("d_phenylC1_7_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]"),

        //Carbon chains

        /**
         * Isoprene unit.
         *
         * <p>SMARTS matches a branched five-carbon motif corresponding to an isoprene
         * building block. Since atomic numbers are used, both aromatic and
         * aliphatic carbon atoms may be matched.</p>
         */
        D_ISOPRENE("d_isoprene_5", "[#6]~[#6](~[#6])~[#6]~[#6]"),
        /**
         * Acetyl group.
         *
         * <p>SMARTS matches a two-carbon chain connected to an oxygen atom. Since
         * atomic numbers are used, both aromatic and aliphatic carbon atoms
         * may be matched.</p>
         */
        D2_ACETYL("d2_acetyl_C2O1", "[#6]~[#6]~[#8]"),
        /**
         * Methylmalonyl motif.
         *
         * <p>SMARTS matches a three-carbon motif terminating in a methyl carbon
         * (degree 1, three implicit hydrogens), corresponding to the
         * methylmalonyl building block.</p>
         */
        D2_METHYLMALONYL("d2_methylmalonyl_C3", "[#6]~[#6][C;D1;h3]"),
        /**
         * Ethyl group.
         *
         * <p>SMARTS matches two carbon atoms connected by a single bond. Since atomic
         * numbers are used, both aromatic and aliphatic carbon atoms may be
         * matched.</p>
         */
        D_ETHYL("d_ethyl_2", "[#6]~[#6]"),
        /**
         * Methyl group.
         *
         * <p>SMARTS matches a terminal methyl carbon (degree 1) carrying three
         * implicit hydrogen atoms.</p>
         */
        D_METHYL("d_methyl_1", "[C;D1;h3]"),

        //Phosphorus- and sulfur-containing functional groups

        /**
         * Phosphate group.
         *
         * <p>SMARTS matches a phosphorus atom connected to an oxygen atom by any bond.
         * Since the atomic numbers ({@code P} and {@code O}) are used, both
         * aromatic and aliphatic environments are matched. The pattern identifies
         * the characteristic phosphorus-oxygen linkage found in phosphate groups,
         * but is not restricted to a specific phosphate species.</p>
         */
        PHOSPHATE("phosphate_2", "P~O"),
        /**
         * Sulfonate group.
         *
         * <p>SMARTS matches a sulfur atom connected to an oxygen atom by any bond.
         * Since the atomic numbers ({@code S} and {@code O}) are used, both
         * aromatic and aliphatic environments are matched. The pattern identifies
         * the characteristic sulfur-oxygen linkage found in sulfonate and related
         * sulfur-oxygen functional groups.</p>
         */
        SULFONATE("sulfonate_2", "S~O"),

        //halogenoids

        /**
         * Fluorine atom.
         * SMARTS matches any fluorine atom.
         */
        HAL_F("hal_f", "[#9]"),
        /**
         * Chlorine atom.
         * SMARTS matches any chlorine atom.
         */
        HAL_CL("hal_cl", "[#17]"),
        /**
         * Bromine atom.
         * SMARTS matches any bromine atom.
         */
        HAL_BR("hal_br", "[#35]"),
        /**
         * Iodine atom.
         * SMARTS matches any iodine atom.
         */
        HAL_I("hal_i", "[#53]"),

        //functional groups

        /**
         * Terminal nitrogen atom.
         *
         * <p>SMARTS matches any nitrogen atom with exactly one bonded neighbour
         * (SMARTS degree = 1), irrespective of whether it is aromatic or
         * aliphatic. This pattern is not restricted to nitrate groups despite
         * the constant name.</p>
         */
        N_NITRATE("n_nitrate_1", "[N;D1]"),
        /**
         * Epoxide group.
         * SMARTS matches an oxygen atom that is part of a three-membered ring.
         */
        O_EPOXY("o_epoxy_1", "[O;x2;r3]"),
        /**
         * Ether group.
         * SMARTS matches a non-cyclic ether oxygen connecting two atoms while
         * excluding esters, phosphates and sulfates.
         */
        O_ETHER("o_ether_1", "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]"),
        /**
         * Hydroxyl group.
         * SMARTS matches a hydroxyl (-OH) group attached to a carbon or nitrogen atom
         * while excluding carboxylic acids, phosphates and sulfates.
         */
        O_HYDROXYL("o_hydroxyl_1", "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]"),

        //rings

        /**
         * Three-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C3("r_c3", "[#6]~1~[#6]~[#6]~1"),
        /**
         * Four-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C4("r_c4", "[#6]~1~[#6]~[#6]~[#6]~1"),
        /**
         * five-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C5("r_c5", "[#6]~1~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Six-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C6("r_c6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Seven-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C7("r_c7", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Eight-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C8("r_c8", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Nine-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C9("r_c9", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Ten-membered carbon ring.
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         */
        R_C10("r_c10", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1");

        /**
         * Human-readable name of the fingerprint feature.
         */
        private final String label;

        /**
         * SMARTS pattern used to identify the corresponding structural motif.
         */
        private final String smarts;

        /**
         * Creates a Biosynfoni fingerprint key definition.
         *
         * <p>Each key consists of a human-readable label and a SMARTS pattern
         * describing a biosynthetically relevant structural motif. During
         * fingerprint generation, the SMARTS pattern is matched against the
         * molecule to determine the presence or frequency of the corresponding
         * feature.</p>
         *
         * @param label  descriptive name of the fingerprint feature
         * @param smarts SMARTS pattern used to identify the structural motif
         */
        DefaultBiosynfoniKey(String label, String smarts) {
            this.label = label;
            this.smarts = smarts;
        }

        /**
         * Returns the human-readable label of this fingerprint feature.
         *
         * @return the feature label
         */
        public String getLabel() {
            return this.label;
        }

        /**
         * Returns the SMARTS pattern defining this fingerprint feature.
         *
         * @return the SMARTS pattern
         */
        public String getSmarts() {
            return this.smarts;
        }
    }

    /**
     * Enables filtering of overlapping matches within the same substructure.
     */
    private final boolean intraSubOverlapToggle;

    /**
     * Enables filtering of overlapping matches between different substructures.
     */
    private final boolean interSubOverLapToggle;

    /**
     * SMARTS patterns used for fingerprint generation.
     */
    private String[] smartsList;

    /**
     * Number of SMARTS patterns.
     */
    private int fingerprintSize;

    /**
     * Creates a Biosynfoni fingerprinter using the default SMARTS patterns with
     * both intra- and inter-substructure overlap filtering disabled. This
     * constructor provides the default configuration of the fingerprinter.
     */
    public BiosynfoniFingerprinter() {
        this(false, false, null);
    }

    /**
     * Creates a Biosynfoni fingerprinter using the default SMARTS patterns and
     * the specified overlap filtering configuration.
     *
     * @param intraSubOverlapToggle whether intra-substructure overlap filtering is enabled
     *                              (default: {@code false})
     * @param interSubOverLapToggle whether inter-substructure overlap filtering is enabled
     *                              (default: {@code false})
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle) {
        this(intraSubOverlapToggle, interSubOverLapToggle, null);
    }

    /**
     * Creates a Biosynfoni fingerprinter with the specified overlap filtering
     * configuration and SMARTS patterns.
     *
     * <p>If {@code smarts} is {@code null}, the default Biosynfoni SMARTS
     * patterns are used. Otherwise, the supplied SMARTS patterns define the
     * fingerprint features used during fingerprint generation.</p>
     *
     * @param intraSubOverlapToggle whether intra-substructure overlap filtering is enabled
     * @param interSubOverLapToggle whether inter-substructure overlap filtering is enabled
     * @param smarts                SMARTS patterns defining the fingerprint. If {@code null}, the
     *                              default Biosynfoni SMARTS patterns are used.
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle, String[] smarts) {
        this.interSubOverLapToggle = interSubOverLapToggle;
        this.intraSubOverlapToggle = intraSubOverlapToggle;
        this.smartsList = smarts;
        if (smarts == null) {
            this.fingerprintSize = DefaultBiosynfoniKey.values().length;
        } else {
            this.fingerprintSize = smarts.length;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {

        if (container == null) {
            throw new NullPointerException("container must not be null");
        }
        BitSet BIOSYNFingerprintBIT = new BitSet(this.fingerprintSize);

        if (container.isEmpty()) {
            return new BitSetFingerprint(BIOSYNFingerprintBIT);
        }

        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);
        for (int i = 0; i < filteredMatches.size(); i++) {
            if (!filteredMatches.get(i).isEmpty()) {
                BIOSYNFingerprintBIT.set(i);
            }

        }
        return new BitSetFingerprint(BIOSYNFingerprintBIT);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {

        if (container == null) {
            throw new NullPointerException("container must not be null");
        }
        if (container.isEmpty()) {
            throw new CDKException("container must not be empty");
        }
        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);

        return new CountFingerprint(filteredMatches);
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {
        int[] counts = new int[this.fingerprintSize];

        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);
        for (int i = 0; i < filteredMatches.size(); i++) {
            counts[i] = filteredMatches.get(i).size();
        }
        Map<String, Integer> rawFingerprint = new HashMap<>();
        int index = 0;
        for (DefaultBiosynfoniKey key : DefaultBiosynfoniKey.values()){
            rawFingerprint.put(key.getLabel(),counts[index]);
            index++;
        }

        return rawFingerprint;
    }

    @Override
    public int getSize() {
        return this.fingerprintSize;
    }

    /**
     * Identify and return SMARTS matches for the molecule, grouped by SMARTS pattern.
     * <p>
     * Purpose:
     * - For each SMARTS pattern (either the default set from {@link DefaultBiosynfoniKey}
     * or a custom {@code smartsList}), find all unique atom-index matches in the
     * supplied molecule and collect them into per-pattern lists.
     * <p>
     * Behavior / Algorithm:
     * - Ensures the molecule is prepared for SMARTS matching (calls
     * {@link SmartsPattern#prepare(IAtomContainer)}) and runs
     * {@link #preprocessMolecule(IAtomContainer)} once to detect atom types,
     * implicit hydrogens, rings and aromaticity.
     * - For each SMARTS, constructs a {@link SmartsPattern} and delegates to
     * {@link #getSubMatches(SmartsPattern, IAtomContainer, List)} to obtain the
     * list of atom-index matches for that pattern. If overlap filters are
     * enabled via constructor flags, {@link #getSubMatches} will apply them.
     * <p>
     * Return value:
     * - A {@code List<List<int[]>>} where the outer list index corresponds to the
     * SMARTS index (default order or the order in {@code smartsList}) and each
     * inner list contains zero or more {@code int[]} arrays with atom indices
     * matching the respective SMARTS.
     * <p>
     * Side effects:
     * - The supplied {@code IAtomContainer} is mutated by
     * {@link #preprocessMolecule(IAtomContainer)} (atom typing, hydrogens,
     * aromaticity flags). This method does not copy the molecule.
     *
     * @param aMolecule the molecule to search for SMARTS matches (mutated)
     * @return grouped SMARTS matches (outer list = SMARTS order, inner lists = matches)
     */
    private List<List<int[]>> getFilteredMatches(IAtomContainer aMolecule) {
        SmartsPattern.prepare(aMolecule);
        List<List<int[]>> filteredMatches = new ArrayList<>(this.fingerprintSize);
        // prepare aromaticity and hydrogen's once
        IAtomContainer preparedMol = this.preprocessMolecule(aMolecule);
        if (this.smartsList == null) {
            for (DefaultBiosynfoniKey key : DefaultBiosynfoniKey.values()) {

                SmartsPattern pattern = SmartsPattern.create(key.getSmarts());
                List<int[]> subMatches = this.getSubMatches(pattern, preparedMol, filteredMatches);
                filteredMatches.add(subMatches);
            }

        } else {
            for (String smarts : this.smartsList) {
                SmartsPattern pattern = SmartsPattern.create(smarts);
                List<int[]> subMatches = this.getSubMatches(pattern, preparedMol, filteredMatches);
                filteredMatches.add(subMatches);
            }
        }
        return filteredMatches;
    }

    /**
     * Find unique atom-index matches for a single SMARTS pattern and apply optional overlap filters.
     * <p>
     * Purpose:
     * - Return all unique atom-index matches for the given {@code pattern} in
     * {@code aMolecule} and apply intra-/inter-pattern overlap filtering when
     * enabled.
     * <p>
     * Algorithm / Notes:
     * - Uses SMARTS matching API to obtain unique atom matches: the underlying
     * call returns an {@code int[][]} where each row contains atom indices for
     * one match. These are converted into a {@code List<int[]>} for easier
     * processing.
     * - If {@link #intraSubOverlapToggle} is {@code true}, calls
     * {@link #intraSubOverlap(List)} to make matches atom-disjoint within the
     * same SMARTS pattern.
     * - If {@link #interSubOverLapToggle} is {@code true}, calls
     * {@link #interSubOverlap(List, List)} to prevent reuse of atoms already
     * accepted by previously processed SMARTS patterns (order-dependent).
     * <p>
     * Edge cases:
     * - If no matches are found the returned list is empty.
     * - The method assumes {@code aMolecule} was preprocessed (aromaticity,
     * hydrogens) by {@link #getFilteredMatches}.
     *
     * @param pattern         a compiled SMARTS pattern
     * @param aMolecule       the (preprocessed) molecule to match against
     * @param filteredMatches previously accepted matches for inter-pattern filtering
     * @return a List of matches where each match is an int[] of atom indices
     */
    private List<int[]> getSubMatches(SmartsPattern pattern, IAtomContainer aMolecule, List<List<int[]>> filteredMatches) {
        List<int[]> subMatches = new ArrayList<>();
        int[][] uniqueMatches = pattern.matchAll(aMolecule).uniqueAtoms().toArray();

        subMatches.addAll(Arrays.asList(uniqueMatches));
        if (this.intraSubOverlapToggle) {
            subMatches = this.intraSubOverlap(subMatches);
        }
        if (this.interSubOverLapToggle) {
            subMatches = this.interSubOverlap(subMatches, filteredMatches);
        }
        return subMatches;
    }

    /**
     * Filters overlapping Matches from the same Structure
     * Matches are processed in sorted order. A match is accepted only if none of its
     * atom indices have already been assigned to a previously accepted match.
     * Accepted matches block all of their atoms from being reused in later
     * matches. This ensures that the returned matches are atom-disjoint.
     * Description
     * - Sorts the provided atom-index matches by their first atom index and,
     * when equal, by their second atom index to ensure deterministic
     * processing order. <p>
     * - Iterates through the sorted matches while maintaining a set of
     * already accepted atom indices ({@code blockedAtoms}). <p>
     * - Uses {@link #hasOverlap(int[], Set)} to determine whether a match
     * shares atoms with previously accepted matches. <p>
     * - Accepts only atom-disjoint matches and records their atom indices via
     * {@link #addBlockedAtoms(int[], Set)}. <p>
     * - Returns a filtered list in which no two matches reuse the same atoms
     * within a single SMARTS pattern.<p>
     *
     * @param subMatches the atom-index matches ({@code int[]}) obtained for a
     *                   single SMARTS pattern
     * @return a list of atom-disjoint matches for the SMARTS pattern
     *
     */
    private List<int[]> intraSubOverlap(List<int[]> subMatches) {
        List<int[]> filteredSubMatches = new ArrayList<>(this.fingerprintSize);

        try{subMatches.sort((a, b) -> {

            if (a[0] != b[0]) {
                return Integer.compare(a[0], b[0]);
            }
            return Integer.compare(a[1], b[1]);
        });}catch(ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException){
            System.err.println("Failed to compare Indices. This should never happen");
        } //#todo Fix the problem with an comparison

        Set<Integer> blockedAtoms = new HashSet<>();
        for (int[] aMatch : subMatches) {

            if (!this.hasOverlap(aMatch, blockedAtoms)) {
                filteredSubMatches.add(aMatch);
                this.addBlockedAtoms(aMatch, blockedAtoms);
            }

        }
        return filteredSubMatches;
    }

    /**
     * Filters matches to prevent atom reuse across different SMARTS patterns.
     * - Collects atom indices from all previously accepted SMARTS matches into
     * a {@code blockedAtoms} set.
     * - Iterates over the current SMARTS matches and checks whether a match
     * reuses atoms that were already assigned to earlier SMARTS patterns.
     * - Uses {@link #hasOverlap(int[], Set)} to determine whether a match
     * overlaps with previously accepted atoms.
     * - Retains only matches that are atom-disjoint from all earlier SMARTS matches.
     * - Returns a filtered list in which atoms assigned to previous SMARTS patterns cannot be reused.
     *
     * @param subMatches  the atom-index matches ({@code int[]}) of the current
     *                    SMARTS pattern
     * @param prevMatches the previously accepted matches from earlier SMARTS
     *                    patterns
     * @return a list of matches that do not overlap with atoms used by
     * previous SMARTS patterns
     */
    private List<int[]> interSubOverlap(List<int[]> subMatches, List<List<int[]>> prevMatches) {
        List<int[]> filteredSubMatches = new ArrayList<>(this.fingerprintSize);
        Set<Integer> blockedAtoms = new HashSet<>();

        for (List<int[]> matches : prevMatches) {
            for (int[] aMatch : matches) {
                this.addBlockedAtoms(aMatch, blockedAtoms);
            }
        }

        for (int[] aMatch : subMatches) {
            if (!this.hasOverlap(aMatch, blockedAtoms)) {
                filteredSubMatches.add(aMatch);

            }
        }
        return filteredSubMatches;
    }

    /**
     * Checks whether a match overlaps with blocked atoms.
     * - Iterates over the atom indices of the provided match. <p>
     * - Checks whether any atom index is contained in the
     * {@code blockedAtoms} set. <p>
     * - Returns {@code true} as soon as an overlapping atom is found;
     * otherwise returns {@code false}. <p>
     *
     * @param aMatch       the atom-index match ({@code int[]}) to be checked
     * @param blockedAtoms the set of atom indices already assigned to
     *                     accepted matches
     * @return {@code true} if the match overlaps with blocked atoms,
     * otherwise {@code false}
     */

    private boolean hasOverlap(int[] aMatch, Set<Integer> blockedAtoms) {
        for (int atom : aMatch) {
            if (blockedAtoms.contains(atom)) {
                return true;
            }
        }
        return false;
    }

    /**
     * - Iterates over the atom indices contained in the provided match. <p>
     * - Adds each atom index to the supplied {@code blockedAtoms} set. <p>
     * - Updates the set in place, ensuring that duplicate atom indices are
     * stored only once. <p>
     *
     * @param match        the atom-index match ({@code int[]}) whose atoms should be
     *                     added to the blocked set
     * @param blockedAtoms the set of blocked atom indices to be updated
     */
    private void addBlockedAtoms(int[] match, Set<Integer> blockedAtoms) {
        for (int atom : match) {
            blockedAtoms.add(atom);
        }
    }

    /**
     * Prepares a molecule for SMARTS matching by assigning atom types,
     * adding implicit hydrogens, detecting ring membership, and
     * perceiving aromaticity.
     * The method uses the recommended Ring finding method for pattern matching({@code Cycles.sssr()})
     * {@code AtomContainerManipulator.perceiveAtomTypesAndConfigureAtoms(IAtomContainer)}
     * and adds implicit hydrogens via {@link CDKHydrogenAdder}.
     * - Resets all aromaticity and ring-membership flags on atoms and bonds
     * before recomputing these properties. <p>
     * - Detects ring membership using the recommended ring-finding method for
     * SMARTS matching, {@code Cycles.sssr(IAtomContainer)}, and marks atoms
     * and bonds belonging to rings.<p>
     * - Creates an {@link Aromaticity} model using
     * {@code ElectronDonation.cdk()} and
     * {@code Cycles.cdkAromaticSet()}.
     * - Applies aromaticity perception to the molecule and updates atom and
     * bond aromaticity flags in place. <p>
     * -Applies unique numbering from {@link #canonicalIndex(IAtomContainer)} if {@link #interSubOverlap(List, List)} shall be applied
     * - Returns the same molecule instance with updated atom typing,
     * hydrogen counts, ring membership, and aromaticity information. <p>
     *
     * @param aMolecule the molecule whose atom information should be improved for later usage
     * @return the same molecule with updated aromaticity
     */
    private IAtomContainer preprocessMolecule(IAtomContainer aMolecule) {
        IAtomContainer clonedMolecule;
        if (this.interSubOverLapToggle) {
            try {
                clonedMolecule = this.canonicalIndex(aMolecule);
            } catch (Exception e) {
                clonedMolecule = aMolecule;
            }
        } else {
            clonedMolecule = aMolecule;
        }
        try {

            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clonedMolecule);
            CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(clonedMolecule.getBuilder());

            hydrogenAdder.addImplicitHydrogens(clonedMolecule);

            for (IAtom atom : clonedMolecule.atoms()) {
                atom.setIsAromatic(false);
                atom.setIsInRing(false);
            }
            for (IBond bond : clonedMolecule.bonds()) {
                bond.setIsAromatic(false);
                bond.setIsInRing(false);
            }
            Cycles cycles = Cycles.sssr(clonedMolecule);
            IRingSet rings = cycles.toRingSet();

            for (IAtomContainer molecule : rings.atomContainers()) {
                for (IAtom atom : molecule.atoms()) {
                    atom.setIsInRing(true);
                }
                for (IBond bond : molecule.bonds()) {
                    bond.setIsInRing(true);
                }
            }
            Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            aromaticity.apply(clonedMolecule);
            return clonedMolecule;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Creates a canonical atom ordering for the given molecule.
     * <p>
     * Uses CDK's canonical labeling algorithm ({@link Canon#label(IAtomContainer, int[][])})
     * to generate a deterministic atom order that is independent of the original
     * atom numbering in the input molecule.
     * </p>
     *
     * <p>
     * The method:
     * <ul>
     *   <li>Computes canonical labels for all atoms.</li>
     *   <li>Sorts atoms according to their canonical labels.</li>
     *   <li>Creates a new molecule containing cloned atoms in canonical order.</li>
     *   <li>Reconstructs all bonds using the new atom indices.</li>
     * </ul>
     * </p>
     *
     * <p>
     * This canonicalization is used to eliminate atom-order dependencies during
     * overlap filtering, especially when
     * {@link #interSubOverlap(List, List)} is enabled.
     * Molecules that are structurally identical but differ only in atom numbering
     * will therefore produce identical atom indices for matching operations.
     * </p>
     *
     * @param aMolecule the molecule to canonicalize
     * @return a new {@link IAtomContainer} with atoms arranged in canonical order
     * @throws CloneNotSupportedException if an atom cannot be cloned
     */
    private IAtomContainer canonicalIndex(IAtomContainer aMolecule) throws CloneNotSupportedException {

        int[][] g = GraphUtil.toAdjList(aMolecule);
        long[] labels = Canon.label(aMolecule, g);

        Integer[] indices = new Integer[labels.length];
        for (int i = 0; i < labels.length; i++) {
            indices[i] = i;
        }
        Arrays.sort(indices, (i, j) -> {
            int cmp = Long.compare(labels[i], labels[j]);
            if (cmp != 0) return cmp;
            return Integer.compare(i, j);
        });


        int[] oldToNew = new int[aMolecule.getAtomCount()];
        IAtomContainer canonical = aMolecule.getBuilder().newInstance(IAtomContainer.class);

        for (int newIdx = 0; newIdx < indices.length; newIdx++) {
            int oldIdx = indices[newIdx];
            IAtom atom = (IAtom) aMolecule.getAtom(oldIdx).clone();
            canonical.addAtom(atom);
            oldToNew[oldIdx] = newIdx;
        }

        for (IBond bond : aMolecule.bonds()) {
            int oldA = aMolecule.indexOf(bond.getBegin());
            int oldB = aMolecule.indexOf(bond.getEnd());

            IAtom newA = canonical.getAtom(oldToNew[oldA]);
            IAtom newB = canonical.getAtom(oldToNew[oldB]);

            IBond newBond = aMolecule.getBuilder().newInstance(IBond.class, newA, newB, bond.getOrder(), bond.getStereo());

            canonical.addBond(newBond);
        }

        return canonical;
    }


}

/**
 * Immutable implementation of {@link ICountFingerprint} storing the occurrence
 * count of each Biosynfoni substructure.
 * <p>
 * The fingerprint is created from the filtered SMARTS matches generated during
 * fingerprint calculation. Each fingerprint position corresponds to one
 * Biosynfoni substructure key, and the stored value represents the number of
 * non-overlapping matches for that key.
 * <p>
 * This implementation provides read-only access to the fingerprint counts.
 * Methods intended for mutable fingerprints are currently not supported.
 */
final class CountFingerprint implements ICountFingerprint {
    private int[] counts;
    private boolean behaveAsBitFingerprint;

    public CountFingerprint(List<List<int[]>> filteredMatches) {
        this.counts = new int[filteredMatches.size()];

        for (int i = 0; i < filteredMatches.size(); i++) {
            counts[i] = filteredMatches.get(i).size();
        }
    }

    @Override
    public long size() {
        return counts.length;
    }

    @Override
    public int numOfPopulatedbins() {
        int populated = 0;
        for (int count : counts) {
            if (count > 0) {
                populated++;
            }
        }
        return populated;
    }

    @Override
    public int getCount(int index) {
        return counts[index];
    }

    @Override
    public int getHash(int index) {
        return index;
    }

    @Override
    public void merge(ICountFingerprint fp) {
        if(size() == fp.size()) {
            for(int i =0; i < fp.size(); i++) {
                counts[i] += fp.getCount(i);
            }
        }else {
            throw new IllegalArgumentException("Fingerprints are not the same size");
        }
    }

    @Override
    public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {
        this.behaveAsBitFingerprint = behaveAsBitFingerprint;
        convertCountToBitFingerprint();
    }

    @Override
    public boolean hasHash(int hash) {
        return hash >= 0 && hash < counts.length;
    }

    @Override
    public int getCountForHash(int hash) {
        return hasHash(hash) ? counts[hash] : 0;
    }
    private void convertCountToBitFingerprint() {

    }
}