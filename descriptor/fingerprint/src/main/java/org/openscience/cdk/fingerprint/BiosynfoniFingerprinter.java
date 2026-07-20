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


import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Implementation of the Biosynfoni molecular fingerprint (<a href="https://doi.org/10.1186/s13321-025-01081-6">
 * Biosynfoni paper (Journal of Chemoinformatics, 2025)</a>), a lightweight,
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
 * Optionally, overlap filtering can be applied, as in the reference
 * implementation:
 * <ul>
 *   <li><b>Intra-pattern filtering</b>: prevents multiple matches of the same
 *       SMARTS pattern from reusing atoms.</li>
 *   <li><b>Inter-pattern filtering</b>: prevents matches of later SMARTS
 *       patterns from reusing atoms already assigned to earlier patterns.</li>
 * </ul>
 * </p>
 *
 * <p>
 * The given molecules will be updated, using {@link SmartsPattern#prepare(IAtomContainer)} to ensure that all properties
 * required are available. This is an Augmentation of the given molecule.
 * </p>
 *
 *   <p><b>Known differences to the reference implementation:</b></p>
 *   <ul>
 *     <li>
 *      The original pattern does not exclude carbon atoms that are part of all six-membered rings because, in SMARTS,
 *      the r{@literal <}x{@literal >} primitive refers to the size of the smallest ring containing an atom.
 *      Consequently, carbon atoms that belong to both a six-membered ring and a smaller fused ring are not matched as
 *      r{@literal <}6{@literal >} and are therefore not excluded.
 *      in CDK did not reproduce the intended semantics of the original
 *      Biosynfoni implementation. The affected fingerprint keys represent phenyl-derived substructures
 *      originating from the shikimate pathway (fingerprint keys 13, 14, and 15)
 *      See the <a href="https://github.com/cdk/cdk/issues/1292#issue-4641968710">  Git Issue </a>
 *      for more information.
 *     </li>
 *   <li>
 *      Fingerprints generated using inter-substructure overlap filtering
 *      differ from the reference implementation(<a href="https://github.com/lucinamay/biosynfoni">
 *      Original Python implementation</a>) because this implementation
 *      resolves overlapping matches using canonical atom numbering rather than
 *      the input atom order. This removes path dependency and ensures
 *      deterministic fingerprints regardless of atom indexing.
 *   </li>
 *   <li>
 *       Unlike the reference Python implementation, count fingerprints are returned as an {@link IntArrayCountFingerprint}.
 *       Consequently, features with a count of zero are omitted from the underlying representation
 *       instead of being stored explicitly. So this Fingerprint is viable with {@link org.openscience.cdk.similarity.Tanimoto}
 *   </li>
 *   <li>
 *       Matching can always differ due to the use of different toolkits and therefore the use of different aromaticity models for example:
 *       The aromaticity assignment of positively charged sulfur atoms may differ
 *       between toolkits, resulting in different SMARTS matching behaviour.
 *   </li>
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
    enum BiosynfoniKey {
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
         * The pattern matches only the amino atom.
         */
        ALL_STD_AMINOS("allstnd_aminos", "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]"),

        /**
         * The SMARTS represented by this pattern describes the common α-amino acid
         * backbone together with a predefined set of allowed side chains. Therefore,
         * it does not match every possible α-amino acid, but only those whose side
         * chain is explicitly included in the SMARTS definition.
         * allowed side chains:
         * <ul>
         *   <li>Alanine</li>
         *   <li>Arginine</li>
         *   <li>Asparagine</li>
         *   <li>Aspartate (Aspartic acid)</li>
         *   <li>Cysteine</li>
         *   <li>Glutamate (Glutamic acid)</li>
         *   <li>Histidine</li>
         *   <li>Isoleucine</li>
         *   <li>Leucine</li>
         *   <li>Lysine</li>
         *   <li>Methionine</li>
         *   <li>Phenylalanine</li>
         *   <li>Serine</li>
         *   <li>Threonine</li>
         *   <li>Tryptophan</li>
         *   <li>Tyrosine</li>
         *   <li>Valine</li>
         *   </ul>
         *  This Pattern will not Match when {@link #ALL_STD_AMINOS} matched before. The pattern matches only the amino atom.
         */
        NON_STD_AMINOS("nonstnd_aminos", "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]"),
        //Sugar derivatives
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
         * SMARTS matches a phosphorus atom connected to an oxygen atom by any bond.
         */
        PHOSPHATE("phosphate_2", "P~O"),
        /**
         * Sulfonate group.
         * SMARTS matches a sulfur atom connected to an oxygen atom by any bond.
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
         * The x2 constraint ensures that the oxygen has exactly two ring bonds,
         * restricting the match to oxygen atoms incorporated into the ring rather than exocyclic oxygen atoms.
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
        BiosynfoniKey(String label, String smarts) {
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
    private final String[] smartsList;

    /**
     * Array holding prepared Patterns
     */
    private final SmartsPattern[] smartsPatterns;

    /**
     * Number of SMARTS patterns.
     */
    private final int fingerprintSize;

    /**
     * logging Tool for possible occurring Errors
     */
    private static final ILoggingTool logger = LoggingToolFactory.createLoggingTool(BiosynfoniFingerprinter.class);

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
     *                              default Biosynfoni SMARTS patterns are used. Still in Code to simplify debugging
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle, String[] smarts) {
        this.interSubOverLapToggle = interSubOverLapToggle;
        this.intraSubOverlapToggle = intraSubOverlapToggle;
        this.smartsList = smarts;
        if (smarts == null) {

            this.fingerprintSize = BiosynfoniKey.values().length;

        } else {
            this.fingerprintSize = smarts.length;
        }

        this.smartsPatterns = new SmartsPattern[this.fingerprintSize];

        int i = 0;
        for (BiosynfoniKey key : BiosynfoniKey.values()) {
            SmartsPattern pattern = SmartsPattern.create(key.getSmarts());
            this.smartsPatterns[i] = pattern;
            i++;
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
        BitSet BiosynfoiBitSet = new BitSet(this.fingerprintSize);

        if (container.isEmpty()) {
            return new BitSetFingerprint(BiosynfoiBitSet);
        }

        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);
        for (int i = 0; i < filteredMatches.size(); i++) {
            if (!filteredMatches.get(i).isEmpty()) {
                BiosynfoiBitSet.set(i);
            }

        }
        return new BitSetFingerprint(BiosynfoiBitSet);
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

        return new BiosynfoniCountFingerprint(filteredMatches);
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {

        if (container == null) {
            throw new NullPointerException("container must not be null");
        }
        if (container.isEmpty()) {
            throw new CDKException("container must not be empty");
        }

        int[] counts = new int[this.fingerprintSize];

        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);
        Map<String, Integer> rawFingerprint = new HashMap<>();
        int index = 0;
        for (BiosynfoniKey key : BiosynfoniKey.values()) {
            counts[index] = filteredMatches.get(index).size();
            rawFingerprint.put(key.getLabel(), counts[index]);
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
     * - For each SMARTS pattern (either the default set from {@link BiosynfoniKey}
     * or a custom {@code smartsList}), find all unique atom-index matches in the
     * supplied molecule and collect them into per-pattern lists.
     * <p>
     * Behavior / Algorithm:
     * - Ensures the molecule is prepared for SMARTS matching (calls
     * {@link SmartsPattern#prepare(IAtomContainer)}) and runs
     * <p>
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
     * (atom typing, hydrogens,
     * aromaticity flags). This method does not copy the molecule.
     *
     * @param aMolecule the molecule to search for SMARTS matches (mutated)
     * @return grouped SMARTS matches (outer list = SMARTS order, inner lists = matches)
     */
    private List<List<int[]>> getFilteredMatches(IAtomContainer aMolecule) {
        SmartsPattern.prepare(aMolecule);
        List<List<int[]>> filteredMatches = new ArrayList<>(this.fingerprintSize);

        if (this.smartsList == null) {
            for (SmartsPattern pattern : this.smartsPatterns) {
                List<int[]> subMatches = this.getSubMatches(pattern, aMolecule, filteredMatches);
                filteredMatches.add(subMatches);
            }

        } else {
            for (String smarts : this.smartsList) {
                SmartsPattern pattern = SmartsPattern.create(smarts);
                List<int[]> subMatches = this.getSubMatches(pattern, aMolecule, filteredMatches);
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
     * {@link #intraSubOverlap(List, int[])} to make matches atom-disjoint within the
     * same SMARTS pattern.
     * - If {@link #interSubOverLapToggle} is {@code true}, calls
     * {@link #interSubOverlap(List, List, int)} to prevent reuse of atoms already
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
        int[][] uniqueMatches = pattern.matchAll(aMolecule).uniqueAtoms().toArray();
        int moleculeAtomCount = aMolecule.getAtomCount();
        List<int[]> subMatches = new ArrayList<>(Arrays.asList(uniqueMatches));
        if (this.intraSubOverlapToggle) {
            subMatches = this.intraSubOverlap(subMatches, getCanonicalIndexMap(aMolecule));
        }
        if (this.interSubOverLapToggle) {
            subMatches = this.interSubOverlap(subMatches, filteredMatches, moleculeAtomCount);
        }
        return subMatches;
    }

    /**
     * Filters overlapping matches within a single SMARTS pattern.
     *
     * <p>Matches are first ordered according to their canonical atom indices to eliminate dependence on the input atom order.
     * The ordered matches are then processed sequentially. A match is accepted
     * only if none of its atoms have already been assigned to a previously
     * accepted match of the same SMARTS pattern.</p>
     *
     * <p>Accepted matches block all of their atoms from being reused by later
     * matches, resulting in a set of atom-disjoint matches for the SMARTS
     * pattern.</p>
     *
     * @param subMatches the matches of a single SMARTS pattern
     * @param newIndex   mapping from the original atom indices to their canonical
     *                   ordering
     * @return the filtered list of non-overlapping matches
     */
    private List<int[]> intraSubOverlap(List<int[]> subMatches, int[] newIndex) {
        List<int[]> filteredSubMatches = new ArrayList<>(this.fingerprintSize);
        for (int[] subMatch : subMatches) {
            for (int i = 0; i < subMatch.length; i++) {
                subMatch[i] = newIndex[subMatch[i]];
            }
        }


        try {
            subMatches.sort((a, b) -> {

                if (a[0] != b[0]) {
                    return Integer.compare(a[0], b[0]);
                }
                return Integer.compare(a[1], b[1]);
            });
        } catch (ArrayIndexOutOfBoundsException arrayIndexOutOfBoundsException) {
            logger.error("Failed to compare indices. This should never happen.", arrayIndexOutOfBoundsException);
        }

        Set<Integer> blockedAtoms = new HashSet<>();
        for (int[] aMatch : subMatches) {

            if (!hasOverlap(aMatch, blockedAtoms)) {
                filteredSubMatches.add(aMatch);
                for (int atomIndex : aMatch) {
                    blockedAtoms.add(atomIndex);
                }
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
    private List<int[]> interSubOverlap(List<int[]> subMatches, List<List<int[]>> prevMatches, int moleculeAtomCount) {
        List<int[]> filteredSubMatches = new ArrayList<>(this.fingerprintSize);
        Set<Integer> blockedAtoms = new HashSet<>(moleculeAtomCount);

        for (List<int[]> matches : prevMatches) {
            for (int[] aMatch : matches) {
                for (int atomIndex : aMatch) {
                    blockedAtoms.add(atomIndex);
                }
            }
        }

        for (int[] aMatch : subMatches) {
            if (!hasOverlap(aMatch, blockedAtoms)) {
                filteredSubMatches.add(aMatch);

            }
        }
        return filteredSubMatches;
    }

    /**
     * Checks whether a SMARTS match reuses atoms that have already been assigned
     * to previously accepted matches.
     *
     * <p>The method iterates over all atom indices in the supplied match and
     * returns {@code true} as soon as one of them is contained in the
     * {@code blockedAtoms} set. Otherwise, {@code false} is returned.</p>
     *
     * @param aMatch       the atom indices of the SMARTS match to check
     * @param blockedAtoms the set of atom indices already assigned to previously
     *                     accepted matches
     * @return {@code true} if the match shares at least one atom with
     * {@code blockedAtoms}; otherwise {@code false}
     */
    private static boolean hasOverlap(int[] aMatch, Set<Integer> blockedAtoms) {
        for (int atomIndex : aMatch) {
            if (blockedAtoms.contains(atomIndex)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Computes a mapping from the original atom indices of a molecule to their
     * canonical indices.
     *
     * <p>The canonical indices are determined using the CDK canonical labelling
     * algorithm. The returned array is indexed by the original atom index, and
     * each value corresponds to the atom's position in the canonical ordering.
     *
     * <p>For an atom with original index {@code i}, its canonical index is given
     * by {@code canonicalIndex[i]}.
     *
     * @param aMolecule the molecule for which the canonical atom index mapping is
     *                  computed
     * @return an array mapping original atom indices to canonical atom indices
     */
    private int[] getCanonicalIndexMap(IAtomContainer aMolecule) {

        int[][] g = GraphUtil.toAdjList(aMolecule);
        long[] labels = Canon.label(aMolecule, g);

        Integer[] indices = new Integer[labels.length];
        for (int i = 0; i < labels.length; i++) {
            indices[i] = i;
        }
        Arrays.sort(indices, Comparator.comparingLong((Integer i) -> labels[i]).thenComparingInt(i -> i));

        int[] canonicalIndex = new int[indices.length];
        for (int canon = 0; canon < indices.length; canon++) {
            canonicalIndex[indices[canon]] = canon;
        }
        return canonicalIndex;
    }


}

