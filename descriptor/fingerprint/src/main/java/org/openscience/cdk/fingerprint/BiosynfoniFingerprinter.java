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
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Implementation of the Biosynfoni molecular fingerprint (<a href="https://doi.org/10.1186/s13321-025-01081-6">
 * Biosynfoni paper (Journal of Chemoinformatics, 2025)</a>), a lightweight,
 * biosynthesis-informed, and interpretable fingerprint specifically designed
 * for natural products and natural product-inspired molecules.
 * Look at the <a href="https://github.com/lucinamay/biosynfoni">
 * Original Python implementation</a> for more context.
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
 * The given {@link IAtomContainer} is not modified. The SMARTS patterns are matched against a copy, in which
 * explicit hydrogen atoms are suppressed (i.e. converted to implicit hydrogen counts, see
 * {@link AtomContainerManipulator#copyAndSuppressedHydrogens(IAtomContainer)}) and on which ring perception and
 * aromaticity detection are applied ({@link SmartsPattern#prepare(IAtomContainer)}). Suppressing the explicit
 * hydrogen atoms is necessary because several patterns use the SMARTS primitives D (number of explicit connections,
 * which includes explicit hydrogen atoms) and h (number of implicit hydrogen atoms), which were written for
 * heavy-atom graphs with implicit hydrogen atoms, as the reference implementation receives them from RDKit. With
 * explicit hydrogen atoms, e.g. a methyl group would not be matched and the oxygen atom of a hydroxyl group would be
 * matched as ether. Like this, the fingerprint of a molecule does not depend on whether its hydrogen atoms are
 * explicit or implicit.
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
 *       The oxygen atoms of nitro groups, nitrate esters, and N-oxides are not matched as hydroxyl groups (see
 *       {@link BiosynfoniKey#O_HYDROXYL}). RDKit converts these groups into their charge-separated form when reading
 *       a molecule, so the reference implementation matches the negatively charged oxygen atom as hydroxyl group.
 *   </li>
 *   <li>
 *      Fingerprints generated using intra-substructure overlap filtering can differ from the
 *      <a href="https://github.com/lucinamay/biosynfoni">original Python implementation</a>. The reference
 *      implementation accepts the matches of a SMARTS pattern greedily in the order in which RDKit finds them, i.e.
 *      along the input atom order, so that its counts can depend on the atom order of the input molecule. This
 *      implementation processes the matches in a canonical atom order instead (a depth-first traversal along the
 *      bonds based on the canonical labels of {@link Canon}), so that the fingerprint does not depend on the input
 *      atom order. If several selections of non-overlapping matches are possible, the selected one can therefore
 *      differ from the one of the reference implementation. The inter-substructure overlap filter does not need a
 *      canonical atom order and gives the same results as the reference implementation, because it blocks the atoms
 *      of all accepted matches of the previous SMARTS patterns and checks each match of the current SMARTS pattern
 *      independently of the other ones.
 *   </li>
 *   <li>
 *       Unlike the reference Python implementation, count fingerprints are returned as an {@link IntArrayCountFingerprint}.
 *       Consequently, features with a count of zero are omitted from the underlying representation
 *       instead of being stored explicitly. This way this Fingerprint is can be used with {@link org.openscience.cdk.similarity.Tanimoto}
 *   </li>
 *   <li>
 *       Matching can always differ due to the use of different toolkits and therefore the use of different aromaticity models, for example:
 *       The aromaticity assignment of positively charged sulfur atoms may differ
 *       between toolkits, resulting in different SMARTS matching behaviour.
 *   </li>
 * </ul>
 * </p>
 *
 * @author Marlon Raffelt (MaRa1778)
 */
public class BiosynfoniFingerprinter extends AbstractFingerprinter {

    /**
     * Default Biosynfoni feature definitions.
     * <p>
     * Each enum constant defines a SMARTS pattern and its corresponding label.
     * Together, these feature definitions form the default Biosynfoni fingerprint.
     * //ToDo add Biosynfonif version and our changes to it
     * </p>
     */
    public enum BiosynfoniKey {
        //Cofactors
        /**
         * Coenzyme A (CoA).
         * <p>
         * SMARTS matches the CoA backbone including the pantetheine chain,
         * diphosphate linker and adenosine moiety.
         * </p>
         */
        CO_COA("co_coa", "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O"),
        /**
         * Nicotinamide adenine dinucleotide (NADH).
         * <p>
         * SMARTS matches the complete reduced NADH cofactor.
         * </p>
         */
        CO_NADH("co_nadh", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O"),// 1: co_nadh
        /**
         * Nicotinamide adenine dinucleotide phosphate (NADPH).
         * <p>
         * SMARTS matches the complete reduced NADPH cofactor.
         * </p>
         */
        CO_NADPH("co_nadph", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O"),
        //Amino acids
        /**
         * Standard proteinogenic amino acids.
         * <p>
         * SMARTS matches the common amino acid backbone together with one of the
         * twenty canonical side chains.
         * The pattern matches only the amino atom.
         * </p>
         */
        ALL_STD_AMINOS("allstnd_aminos", "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]"),
        /**
         * Looser proteinogenic amino acids.
         * <p>
         * SMARTS matches the common amino acid backbone together with one of the
         * twenty canonical side chains.
         * The pattern matches only the amino atom.
         * This key is looser to the bonds connected to it and allows it to be inside a structure
         * </p>
         */
        NON_STD_AMINOS("nonstnd_aminos", "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]"),
        //Sugar derivatives
        /**
         * Open-chain hexose.
         * <p>
         * SMARTS matches a linear six-carbon sugar containing a hydroxyl group
         * on every carbon atom.
         * </p>
         */
        S_OPENPYR_C6O6("s_openpyr_C6O6", "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])"),
        /**
         * Open-chain pentose.
         * <p>
         * SMARTS matches a linear five-carbon sugar containing hydroxyl groups.
         * </p>
         */
        S_OPENFUR_C5O5("s_openfur_C5O5", "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])"),
        /**
         * Pyranose sugar.
         * <p>
         * SMARTS matches a six-membered sugar ring containing one oxygen atom.
         * </p>
         */
        S_PYRANOSE_C5O4("s_pyranose_C5O4", "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1"),
        /**
         * Furanose sugar.
         * <p>
         * SMARTS matches a five-membered sugar ring containing one oxygen atom.
         * </p>
         */
        S_FURANOSE_C4O3("s_furanose_C4O3", "C~1~[#8]~C~C(~[#8])~C(~[#8])~1"),
        //alkaloids
        /**
         * Indole ring system.
         * <p>
         * SMARTS matches an indole scaffold with a two-carbon side chain terminating
         * in a nitrogen atom.
         * </p>
         */
        D_INDOLE("d_indoleC2N_12", "c1cccc2c1c(~[#6]~[#6]~[#7])cn2"),
        /**
         * Phenethylamine motif.
         * <p>
         * SMARTS matches a benzene ring connected to a two-carbon aliphatic chain
         * ending in a nitrogen atom.
         * </p>
         */
        D_PHENYL_C2N("d_phenylC2N_9", "c1ccccc1[#6][#6][#7]"),
        /**
         * Six-membered nitrogen-containing ring.
         * <p>
         * SMARTS matches a six-membered ring consisting of five carbon atoms and one
         * nitrogen atom, regardless of aromaticity.
         * </p>
         */
        D_C5N("d_c5n_6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1"),
        /**
         * Five-membered nitrogen-containing ring.
         * <p>
         * SMARTS matches a five-membered ring consisting of four carbon atoms and one
         * nitrogen atom, regardless of aromaticity.
         * </p>
         */
        D_C4N("d_c4n_5", "[#6]~1~[#6]~[#6]~[#6]~[#7]~1"),
        //Phenyls from Shikimate pathway
        /**
         * Phenylpropyl motif.
         * <p>
         * SMARTS matches a phenyl ring attached to a three-carbon side chain.
         * </p>
         */
        D_PHENYL_C3("d_phenylC3_9_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]~[#6!$(*1*****1)]"),
        /**
         * Phenethyl motif.
         * <p>
         * SMARTS matches a phenyl ring attached to a two-carbon side chain.
         * </p>
         */
        D_PHENYL_C2("d_phenylC2_8_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]"),
        /**
         * Benzyl motif.
         * <p>
         * SMARTS matches a phenyl ring attached to a single carbon atom.
         * </p>
         */
        D_PHENYL_C1("d_phenylC1_7_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]"),
        //Carbon chains
        /**
         * Isoprene unit.
         * <p> SMARTS matches a branched five-carbon motif corresponding to an isoprene
         * building block. Since atomic numbers are used, both aromatic and
         * aliphatic carbon atoms may be matched.</p>
         */
        D_ISOPRENE("d_isoprene_5", "[#6]~[#6](~[#6])~[#6]~[#6]"),
        /**
         * Acetyl group.
         * <p> SMARTS matches a two-carbon chain connected to an oxygen atom. Since
         * atomic numbers are used, both aromatic and aliphatic carbon atoms
         * may be matched.</p>
         */
        D2_ACETYL("d2_acetyl_C2O1", "[#6]~[#6]~[#8]"),
        /**
         * Methylmalonyl motif.
         * <p> SMARTS matches a three-carbon motif terminating in a methyl carbon
         * (degree 1, three implicit hydrogens), corresponding to the
         * methylmalonyl building block.</p>
         */
        D2_METHYLMALONYL("d2_methylmalonyl_C3", "[#6]~[#6][C;D1;h3]"),
        /**
         * Ethyl group.
         * <p> SMARTS matches two carbon atoms connected by a single bond. Since atomic
         * numbers are used, both aromatic and aliphatic carbon atoms may be
         * matched.</p>
         */
        D_ETHYL("d_ethyl_2", "[#6]~[#6]"),
        /**
         * Methyl group.
         * <p> SMARTS matches a terminal methyl carbon (degree 1) carrying three
         * implicit hydrogen atoms.</p>
         */
        D_METHYL("d_methyl_1", "[C;D1;h3]"),
        //Phosphorus- and sulfur-containing functional groups
        /**
         * Phosphate group.
         * <p> SMARTS matches a phosphorus atom connected to an oxygen atom by any bond.
         * </p>
         */
        PHOSPHATE("phosphate_2", "P~O"),
        /**
         * Sulfonate group.
         * <p> SMARTS matches a sulfur atom connected to an oxygen atom by any bond.
         * </p>
         */
        SULFONATE("sulfonate_2", "S~O"),
        //halide
        /**
         * Fluorine atom.
         * <p> SMARTS matches any fluorine atom.
         * </p>
         */
        HAL_F("hal_f", "[#9]"),
        /**
         * Chlorine atom.
         * <p> SMARTS matches any chlorine atom.
         * </p>
         */
        HAL_CL("hal_cl", "[#17]"),
        /**
         * Bromine atom.
         * <p> SMARTS matches any bromine atom.
         * </p>
         */
        HAL_BR("hal_br", "[#35]"),
        /**
         * Iodine atom.
         * <p> SMARTS matches any iodine atom.
         * </p>
         */
        HAL_I("hal_i", "[#53]"),
        //functional groups
        /**
         * Terminal nitrogen atom.
         * <p> SMARTS matches any aliphatic nitrogen atom with exactly one bonded neighbor
         * (SMARTS degree = 1), This pattern is not restricted to nitrate groups despite
         * the constant name.</p>
         */
        N_NITRATE("n_nitrate_1", "[N;D1]"),
        /**
         * Epoxide group.
         * <p> SMARTS matches an oxygen atom that is part of a three-membered ring.
         * The x2 constraint ensures that the oxygen has exactly two ring bonds,
         * restricting the match to oxygen atoms incorporated into the ring rather than exocyclic oxygen atoms.
         * </p>
         */
        O_EPOXY("o_epoxy_1", "[O;x2;r3]"),
        /**
         * Ether group.
         * <p>
         * SMARTS matches a non-cyclic ether oxygen connecting two atoms while
         * excluding esters, phosphates and sulfates.
         * </p>
         */
        O_ETHER("o_ether_1", "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]"),
        /**
         * Hydroxyl group.
         * <p>
         * SMARTS matches a hydroxyl (-OH) group or a deprotonated hydroxyl group (-O<sup>-</sup>) attached to a carbon
         * or nitrogen atom while excluding carboxylic acids, phosphates and sulfates.
         * </p>
         * <p>
         * In addition to the original pattern, oxygen atoms bonded to a nitrogen atom that carries a further oxygen
         * atom ({@code !$(*[#7]~[#8])}), i.e. the oxygen atoms of nitro groups and nitrate esters, and negatively
         * charged oxygen atoms bonded to a positively charged nitrogen atom ({@code !$([#8-][#7+])}), i.e. the oxygen
         * atoms of N-oxides, are excluded. These groups can be written in a pentavalent form ({@code N(=O)=O},
         * {@code N=O}), in which the oxygen atom is not matched because it is doubly bonded, or in a charge-separated
         * form ({@code [N+](=O)[O-]}, {@code [N+][O-]}), in which the singly bonded {@code [O-]} would otherwise be
         * matched as deprotonated hydroxyl group, so that the result depended on the representation of the same
         * molecule. In neither form, these oxygen atoms are (deprotonated) hydroxyl groups: the negative charge is a
         * formal charge and does not result from the removal of a proton. The first expression applies the rule of
         * the original pattern for carboxylic acids and carboxylates ({@code !$(*C~O)}) to nitrogen atoms.
         * Hydroxylamines, oximes, and hydroxamic acids, whose uncharged nitrogen atom carries only one oxygen atom,
         * are still matched.
         * </p>
         */
        O_HYDROXYL("o_hydroxyl_1", "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(*[#7]~[#8]);!$([#8-][#7+]);!$(P);!$(S)]"),
        //rings
        /**
         * Three-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity.
         * </p>
         */
        R_C3("r_c3", "[#6]~1~[#6]~[#6]~1"),
        /**
         * Four-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of three carbon atoms, regardless of aromaticity. </p>
         */
        R_C4("r_c4", "[#6]~1~[#6]~[#6]~[#6]~1"),
        /**
         * Five-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of four carbon atoms, regardless of aromaticity. </p>
         */
        R_C5("r_c5", "[#6]~1~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Six-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of five carbon atoms, regardless of aromaticity. </p>
         */
        R_C6("r_c6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Seven-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of six carbon atoms, regardless of aromaticity. </p>
         */
        R_C7("r_c7", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Eight-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of seven carbon atoms, regardless of aromaticity. </p>
         */
        R_C8("r_c8", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Nine-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of eight carbon atoms, regardless of aromaticity. </p>
         */
        R_C9("r_c9", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        /**
         * Ten-membered carbon ring.
         * <p>
         * SMARTS matches a ring consisting of nine carbon atoms, regardless of aromaticity. </p>
         */
        R_C10("r_c10", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1");

        /**
         * Readable name of the fingerprint feature.
         */
        private final String label;

        /**
         * SMARTS pattern used to identify the corresponding structural motif.
         */
        private final String smarts;

        /**
         * Creates a Biosynfoni fingerprint key definition.
         *
         * <p>Each key consists of a readable label and a SMARTS pattern
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
         * Returns the readable label of this fingerprint feature.
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
    private final boolean interSubOverlapToggle;

    /**
     * Array holding prepared Patterns.
     */
    private final SmartsPattern[] smartsPatterns;

    /**
     * Number of SMARTS patterns.
     */
    private static final int FINGERPRINT_SIZE = BiosynfoniKey.values().length;

    /**
     * Creates a Biosynfoni fingerprinter using the default SMARTS patterns with
     * both intra- and inter-substructure overlap filtering disabled. This
     * constructor provides the default configuration of the fingerprinter.
     */
    public BiosynfoniFingerprinter() {
        this(false, false);
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
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle) {
        this.interSubOverlapToggle = interSubOverLapToggle;
        this.intraSubOverlapToggle = intraSubOverlapToggle;

        this.smartsPatterns = new SmartsPattern[BiosynfoniFingerprinter.FINGERPRINT_SIZE];

        int i = 0;
        for (BiosynfoniKey key : BiosynfoniKey.values()) {
            //the molecule is prepared once for all patterns in getFilteredMatches(), so the automatic
            // preparation of every pattern (ring perception and aromaticity detection) is disabled
            SmartsPattern pattern = SmartsPattern.create(key.getSmarts()).setPrepare(false);
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
        BitSet biosynfoiBitSet = new BitSet(BiosynfoniFingerprinter.FINGERPRINT_SIZE);

        if (container.isEmpty()) {
            return new BitSetFingerprint(biosynfoiBitSet);
        }

        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);
        for (int i = 0; i < filteredMatches.size(); i++) {
            if (!filteredMatches.get(i).isEmpty()) {
                biosynfoiBitSet.set(i);
            }

        }
        return new BitSetFingerprint(biosynfoiBitSet);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
        if (container == null) {
            throw new NullPointerException("container must not be null");
        }

        List<List<int[]>> filteredMatches = new ArrayList<>(BiosynfoniFingerprinter.FINGERPRINT_SIZE);

        if (container.isEmpty()) {
            return new BiosynfoniCountFingerprint(filteredMatches);
        }

        filteredMatches = this.getFilteredMatches(container);

        return new BiosynfoniCountFingerprint(filteredMatches);
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {

        if (container == null) {
            throw new NullPointerException("container must not be null");
        }
        if (container.isEmpty()) {
            return new HashMap<>((BiosynfoniFingerprinter.FINGERPRINT_SIZE + 1) * 4 / 3);
        }

        List<List<int[]>> filteredMatches = this.getFilteredMatches(container);
        Map<String, Integer> rawFingerprint = new HashMap<>((BiosynfoniFingerprinter.FINGERPRINT_SIZE + 1) * 4 / 3);

        int index = 0;
        for (BiosynfoniKey key : BiosynfoniKey.values()) {
            rawFingerprint.put(key.getLabel(), filteredMatches.get(index).size());
            index++;
        }
        return rawFingerprint;
    }

    @Override
    public int getSize() {
        return BiosynfoniFingerprinter.FINGERPRINT_SIZE;
    }

    /**
     * Finds and filters the substructure matches for all SMARTS patterns
     * defined for this fingerprint.
     *
     * <p>The patterns are matched against a copy of the molecule in which explicit
     * hydrogen atoms are suppressed (see the class documentation), so that the given
     * molecule is not modified. The copy is prepared once for SMARTS pattern
     * matching (ring perception and aromaticity detection). If the
     * intra-substructure overlap filter is enabled, a canonical atom index
     * mapping is generated to make the filter independent of the input atom
     * order. Each SMARTS pattern is then matched against the copy, and the
     * resulting matches are filtered before being added to the result list. If
     * the inter-substructure overlap filter is enabled, the atoms of the accepted
     * matches are blocked for all following SMARTS patterns.</p>
     *
     * <p>The atom indices of the returned matches refer to the copy. They are identical
     * to the atom indices of the given molecule if it does not contain explicit hydrogen
     * atoms that are suppressed.</p>
     *
     * @param container the molecule to search for SMARTS pattern matches
     * @return a list containing the filtered matches for each SMARTS pattern
     */
    private List<List<int[]>> getFilteredMatches(IAtomContainer container) {
        IAtomContainer molecule = AtomContainerManipulator.copyAndSuppressedHydrogens(container);
        //marks ring atoms and bonds and applies the Daylight aromaticity model for all patterns
        SmartsPattern.prepare(molecule);
        int[] canonicalIndex = null;
        if (this.intraSubOverlapToggle) {
            canonicalIndex = this.getCanonicalIndexMap(molecule);
        }
        //atoms of the accepted matches of the previous SMARTS patterns, only used by the inter-substructure overlap filter
        boolean[] blockedByPreviousPatterns = new boolean[molecule.getAtomCount()];
        List<List<int[]>> filteredMatches = new ArrayList<>(BiosynfoniFingerprinter.FINGERPRINT_SIZE);
        for (SmartsPattern pattern : this.smartsPatterns) {
            List<int[]> subMatches = this.getSubMatches(pattern, molecule, canonicalIndex, blockedByPreviousPatterns);
            if (this.interSubOverlapToggle) {
                for (int[] subMatch : subMatches) {
                    for (int atomIndex : subMatch) {
                        blockedByPreviousPatterns[atomIndex] = true;
                    }
                }
            }
            filteredMatches.add(subMatches);
        }
        return filteredMatches;
    }

    /**
     * Finds all unique atom matches for a SMARTS pattern and applies the configured overlap filters,
     *
     * <p>First, all unique atom matches of the pattern against the molecule
     * are determined. Depending on the configured overlap options, matches
     * are then filtered for overlaps either within the same SMARTS pattern
     * or with matches identified by previously processed SMARTS patterns.</p>
     *
     * @param pattern                   the SMARTS pattern used for substructure matching
     * @param aMolecule                 the molecule to search for matches
     * @param canonicalIndex            the canonical atom index mapping used for intra-pattern overlap filtering;
     *                                  null if the intra-substructure overlap filter is disabled
     * @param blockedByPreviousPatterns flags of the atoms of the accepted matches of the previously processed SMARTS
     *                                  patterns, indexed by atom index, used for inter-pattern overlap filtering
     * @return a list of filtered atom matches for the given SMARTS pattern
     */
    private List<int[]> getSubMatches(SmartsPattern pattern, IAtomContainer aMolecule, int[] canonicalIndex,
                                      boolean[] blockedByPreviousPatterns) {
        int[][] uniqueMatches = pattern.matchAll(aMolecule).uniqueAtoms().toArray();
        List<int[]> subMatches = new ArrayList<>(Arrays.asList(uniqueMatches));
        if (this.intraSubOverlapToggle) {
            subMatches = this.intraSubOverlap(subMatches, canonicalIndex);
        }
        if (this.interSubOverlapToggle) {
            subMatches = this.interSubOverlap(subMatches, blockedByPreviousPatterns);
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
     * pattern. The canonical indices are only used to determine the processing
     * order; the returned matches contain the original atom indices, like the
     * matches of the other filters.</p>
     *
     * @param subMatches     the matches of a single SMARTS pattern
     * @param canonicalIndex mapping from the original atom indices to their canonical
     *                       ordering
     * @return the filtered list of non-overlapping matches
     */
    private List<int[]> intraSubOverlap(List<int[]> subMatches, int[] canonicalIndex) {
        List<int[]> filteredSubMatches = new ArrayList<>(subMatches.size());
        List<int[]> canonMatches = new ArrayList<>(subMatches.size());

        for (int[] subMatch : subMatches) {
            int[] canonMatch = new int[subMatch.length];

            for (int i = 0; i < subMatch.length; i++) {
                canonMatch[i] = canonicalIndex[subMatch[i]];
            }

            Arrays.sort(canonMatch);
            canonMatches.add(canonMatch);
        }
        //the positions of the matches are sorted by their canonical indices, so that the original matches are kept
        Integer[] order = new Integer[subMatches.size()];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }
        //the reference implementation processes the matches in the order in which RDKit finds them, i.e. along the
        // input atom order (its sorting by match length does not change the order of matches of the same pattern);
        // here, the canonical atom order is used instead, which also follows the bonds, see getCanonicalIndexMap()
        Arrays.sort(order, (aPosition, bPosition) -> {
            int[] a = canonMatches.get(aPosition);
            int[] b = canonMatches.get(bPosition);
            int length = Math.min(a.length, b.length);
            for (int i = 0; i < length; i++) {
                int comparison = Integer.compare(a[i], b[i]);
                if (comparison != 0) {
                    return comparison;
                }
            }
            return Integer.compare(a.length, b.length);
        });
        boolean[] blockedAtoms = new boolean[canonicalIndex.length];
        for (int position : order) {
            int[] subMatch = subMatches.get(position);
            if (!BiosynfoniFingerprinter.hasOverlap(subMatch, blockedAtoms)) {
                filteredSubMatches.add(subMatch);
                for (int atomIndex : subMatch) {
                    blockedAtoms[atomIndex] = true;
                }
            }
        }

        return filteredSubMatches;
    }


    /**
     * Filters matches to prevent atom reuse across different SMARTS patterns.
     * Retains only the matches of the current SMARTS pattern that do not share any atom with the accepted
     * matches of the previously processed SMARTS patterns. The flags of these atoms are maintained in
     * {@link #getFilteredMatches(IAtomContainer)}, so they do not have to be collected again for every pattern.
     *
     * @param subMatches                the atom-index matches ({@code int[]}) of the current
     *                                  SMARTS pattern
     * @param blockedByPreviousPatterns flags of the atoms of the accepted matches of the previously processed SMARTS
     *                                  patterns, indexed by atom index
     * @return a list of matches that do not overlap with atoms used by
     * previous SMARTS patterns
     */
    private List<int[]> interSubOverlap(List<int[]> subMatches, boolean[] blockedByPreviousPatterns) {
        List<int[]> filteredSubMatches = new ArrayList<>(subMatches.size());
        for (int[] aMatch : subMatches) {
            if (!BiosynfoniFingerprinter.hasOverlap(aMatch, blockedByPreviousPatterns)) {
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
     * returns {@code true} as soon as one of them is flagged in
     * {@code blockedAtoms}. Otherwise, {@code false} is returned.</p>
     *
     * @param aMatch       the atom indices of the SMARTS match to check
     * @param blockedAtoms flags of the atoms already assigned to previously
     *                     accepted matches, indexed by atom index
     * @return {@code true} if the match shares at least one atom with
     * {@code blockedAtoms}; otherwise {@code false}
     */
    private static boolean hasOverlap(int[] aMatch, boolean[] blockedAtoms) {
        for (int atomIndex : aMatch) {
            if (blockedAtoms[atomIndex]) {
                return true;
            }
        }
        return false;
    }

    /**
     * Computes a mapping from the original atom indices of a molecule to their
     * canonical indices, i.e. a numbering of the atoms that does not depend on
     * the input atom order.
     *
     * <p>The atoms are numbered in the order of a depth-first traversal along the
     * bonds that starts at the atom with the smallest canonical label of the CDK
     * canonical labelling ({@link Canon#label(IAtomContainer, int[][])}) and visits
     * the neighbours of an atom in the order of their canonical labels. Disconnected
     * components are traversed one after the other, in the order of their smallest
     * canonical label. This follows the same principle as the atom order of the unique
     * SMILES of the molecule, without depending on the SMILES module. In contrast to
     * numbering the atoms directly by their canonical labels, neighbouring atoms get
     * close indices, so that the intra-substructure overlap filter selects the matches
     * along chains and rings. Like this, the matches are packed about as densely as by
     * the reference implementation, which processes them along the input atom order,
     * e.g. two instead of one ethyl unit in the chain of four carbon atoms of
     * thiophene 1,1-dioxide. Stereochemistry is not considered, as in
     * the SMARTS patterns. The traversal is iterative to support large molecules.</p>
     *
     * <p>The returned array is indexed by the original atom index: for an atom with
     * original index {@code i}, its canonical index is given by {@code canonicalIndex[i]}.
     * The mapping is independent of the overlap filter settings; it is only used by the
     * intra-substructure overlap filter. The method is package-private for testing.</p>
     *
     * @param aMolecule the molecule for which the canonical atom index mapping is
     *                  computed
     * @return an array mapping original atom indices to canonical atom indices
     */
    int[] getCanonicalIndexMap(IAtomContainer aMolecule) {
        int[][] g = GraphUtil.toAdjList(aMolecule);
        long[] labels = Canon.label(aMolecule, g);

        Integer[] atomsByLabel = new Integer[labels.length];
        for (int i = 0; i < labels.length; i++) {
            atomsByLabel[i] = i;
        }
        Arrays.sort(atomsByLabel, Comparator.comparingLong((Integer i) -> labels[i]));
        //neighbours are pushed in descending label order, so that the one with the smallest label is visited first
        Comparator<Integer> descendingLabel = Comparator.comparingLong((Integer i) -> labels[i]).reversed();

        int[] canonicalIndex = new int[labels.length];
        Arrays.fill(canonicalIndex, -1);
        int nextIndex = 0;
        Deque<Integer> stack = new ArrayDeque<>();
        for (int start : atomsByLabel) {
            if (canonicalIndex[start] >= 0) {
                continue;
            }
            stack.push(start);
            while (!stack.isEmpty()) {
                int atom = stack.pop();
                if (canonicalIndex[atom] >= 0) {
                    continue;
                }
                canonicalIndex[atom] = nextIndex++;
                Integer[] neighbours = new Integer[g[atom].length];
                for (int k = 0; k < neighbours.length; k++) {
                    neighbours[k] = g[atom][k];
                }
                Arrays.sort(neighbours, descendingLabel);
                for (int neighbour : neighbours) {
                    if (canonicalIndex[neighbour] < 0) {
                        stack.push(neighbour);
                    }
                }
            }
        }
        return canonicalIndex;
    }
}

/**
 * Immutable implementation of {@link IntArrayCountFingerprint} storing the occurrence
 * count of each Biosynfoni substructure.
 * <p>
 * This class provides a count-based fingerprint representation compatible with
 * Tanimoto similarity calculations.
 * <p>
 * The fingerprint is created from the filtered SMARTS matches generated during
 * fingerprint calculation. Each fingerprint position corresponds to one Biosynfoni
 * substructure key, and the stored value represents the number of matches for that key.
 * <p>
 * This implementation provides read-only access to the fingerprint counts.
 * Methods intended for mutable fingerprints are currently not supported.
 */
final class BiosynfoniCountFingerprint extends IntArrayCountFingerprint {
    public BiosynfoniCountFingerprint(List<List<int[]>> filteredMatches) {
        int[] counts = new int[filteredMatches.size()];

        for (int i = 0; i < filteredMatches.size(); i++) {
            counts[i] = filteredMatches.get(i).size();
        }

        List<Integer> hashes = new ArrayList<>(counts.length);
        List<Integer> values = new ArrayList<>(counts.length);

        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > 0) {
                hashes.add(i);
                values.add(counts[i]);
            }
        }

        this.hitHashes = new int[hashes.size()];
        this.numOfHits = new int[values.size()];

        for (int i = 0; i < hashes.size(); i++) {
            this.hitHashes[i] = hashes.get(i);
            this.numOfHits[i] = values.get(i);
        }
    }
}