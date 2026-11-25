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
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
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
 * <p><b>Original and refined fingerprint:</b> several SMARTS patterns of the reference implementation (version
 * {@code full_1103}) match more substructures than their name and their explanation in the reference implementation
 * describe, and some patterns leave room for chemical interpretation. Therefore, two pattern sets are available, which
 * are chosen for the whole fingerprinter at once via the constructor parameter {@code isOriginalFingerprint} (see
 * {@link #isOriginalFingerprint()}):
 * </p>
 * <ul>
 *   <li><b>Refined fingerprint</b> (default, {@code isOriginalFingerprint = false}): the refined patterns follow the
 *   key names more closely. They differ from the original ones in 15 SMARTS patterns and in the matching of the 8 ring
 *   keys, which are matched by ring perception (see {@link BiosynfoniKey#getRingSize()}); the other 16 keys are
 *   identical. The fingerprint keeps its 39 positions and their order.</li>
 *   <li><b>Original fingerprint</b> ({@code isOriginalFingerprint = true}): the original patterns of the reference
 *   implementation, which reproduce its counts apart from the differences listed below.</li>
 * </ul>
 * <p>The refinements, which are described in detail at the respective keys (see {@link BiosynfoniKey}):</p>
 * <ul>
 *   <li><b>Patterns that match more than their name describes</b>: the open-chain sugar patterns also match sugar
 *   rings ({@link BiosynfoniKey#S_OPENPYR_C6O6}, {@link BiosynfoniKey#S_OPENFUR_C5O5}); the non-standard amino acid
 *   pattern lists the standard side chains instead of a generic one ({@link BiosynfoniKey#NON_STD_AMINOS}); the
 *   phenethylamine pattern matches the pyrrole ring of indoles ({@link BiosynfoniKey#D_PHENYL_C2N}); the hydroxyl
 *   pattern ignores hemiacetal hydroxyl groups and matches the oxygen atoms of nitro groups and N-oxides depending on
 *   their notation ({@link BiosynfoniKey#O_HYDROXYL}); the phosphate and sulfonate patterns count every
 *   phosphorus-oxygen and sulfur-oxygen bond instead of every group ({@link BiosynfoniKey#PHOSPHATE},
 *   {@link BiosynfoniKey#SULFONATE}); the NADH pattern also matches NADPH ({@link BiosynfoniKey#CO_NADH}); the side
 *   chain condition of the phenyl-C3 and phenyl-C2 patterns depends on the smallest ring of an atom
 *   ({@link BiosynfoniKey#D_PHENYL_C3}); the terminal nitrogen pattern also matches nitriles and imines
 *   ({@link BiosynfoniKey#N_NITRATE}); the ether pattern matches peroxides ({@link BiosynfoniKey#O_ETHER}); and the
 *   ring patterns also match the perimeter of fused and bridged ring systems, e.g. a ten-membered ring in naphthalene
 *   ({@link BiosynfoniKey#getRingSize()}).</li>
 *   <li><b>Interpretations</b>: the rings of the phenyl units from the shikimate pathway have to be aromatic
 *   ({@link BiosynfoniKey#D_PHENYL_C3}), the methyl group of a methylmalonyl unit has to be attached to a branching
 *   point ({@link BiosynfoniKey#D2_METHYLMALONYL}), and standard amino acids that are N-alkylated or esterified count
 *   as standard amino acids ({@link BiosynfoniKey#ALL_STD_AMINOS}).</li>
 * </ul>
 * <p>Deliberately broad patterns of the reference implementation, e.g. for isoprene, acetyl, and ethyl units, are
 * kept, as is the definition of the sugar rings by their number of oxygen atoms. The same pattern sets are available in
 * the Biosynfoni fragmenter of <a href="https://github.com/FelixBaensch/MORTAR">MORTAR</a>, where the refined patterns
 * are called pattern variants.</p>
 *
 *   <p><b>Known differences to the reference implementation</b>, independent of the pattern set:</p>
 *   <ul>
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
 *       between toolkits, resulting in different SMARTS matching behaviour. RDKit also converts nitro groups and
 *       N-oxides into their charge-separated form when reading a molecule, so that the original hydroxyl pattern
 *       matches their negatively charged oxygen atom there, but in the CDK only if the molecule is given in this
 *       form (see {@link BiosynfoniKey#O_HYDROXYL}).
 *   </li>
 * </ul>
 * </p>
 *
 * @author Marlon Raffelt (MaRa1778)
 */
public class BiosynfoniFingerprinter extends AbstractFingerprinter {

    /**
     * Biosynfoni substructure keys in the order of the fingerprint positions.
     * <p>
     * The keys and their order correspond to the substructure set {@code full_1103}, the default version of the
     * <a href="https://github.com/lucinamay/biosynfoni">original Python implementation</a>. The order is relevant for
     * the inter-substructure overlap filter, because the matches of earlier keys block their atoms for all later keys.
     * </p>
     * <p>
     * Each key holds its original SMARTS pattern of the reference implementation ({@link #getOriginalSmarts()}) and
     * its refined SMARTS pattern ({@link #getRefinedSmarts()}), which is identical to the original one for 16 keys.
     * The ring keys are matched by ring perception instead of their SMARTS pattern in the refined fingerprint, see
     * {@link #getRingSize()}. The documentation of each key describes its refinement. Which patterns are used is
     * set for the whole fingerprinter, see {@link BiosynfoniFingerprinter#isOriginalFingerprint()}.
     * </p>
     */
    public enum BiosynfoniKey {
        //Cofactors
        /**
         * Coenzyme A (CoA).
         * <p>
         * SMARTS matches the CoA backbone including the pantetheine chain,
         * diphosphate linker and adenosine moiety, also in acyl-CoA thioesters.
         * </p>
         */
        CO_COA("co_coa", "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O"),
        /**
         * Nicotinamide adenine dinucleotide (NADH).
         * <p>
         * SMARTS matches the complete reduced NADH cofactor, also the oxidised form NAD<sup>+</sup>.
         * </p>
         * <p>
         * Refined pattern: the 2'-oxygen atom of the adenosine ribose must not be bonded to a phosphorus atom
         * ({@code [O;!$(*P)]}). The original pattern also matches NADPH, so that NADPH counts as NADH as well, and with
         * the inter-substructure overlap filter, the following NADPH key can never match.
         * </p>
         */
        CO_NADH("co_nadh",
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O",
                "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~[O;!$(*P)])~O)~O)~O"),
        /**
         * Nicotinamide adenine dinucleotide phosphate (NADPH).
         * <p>
         * SMARTS matches the complete reduced NADPH cofactor, also the oxidised form NADP<sup>+</sup>.
         * </p>
         */
        CO_NADPH("co_nadph", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O"),
        //Amino acids
        /**
         * Standard proteinogenic amino acids (free or peptide-bound).
         * <p>
         * SMARTS matches the common amino acid backbone together with one of the
         * twenty canonical side chains (pattern taken from the Daylight SMARTS examples).
         * The pattern matches only the amino atom.
         * </p>
         * <p>
         * Refined pattern: residues with a standard side chain whose amino nitrogen atom carries a further carbon
         * substituent, e.g. N-methylated residues in peptides, or whose carboxy group is esterified, e.g. in
         * depsipeptides, are matched as well, since N-methylation and esterification are modifications of standard
         * amino acids. The original pattern does not match them.
         * </p>
         */
        ALL_STD_AMINOS("allstnd_aminos",
                "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]",
                "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H2][CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2,OX1-,N])]"),
        /**
         * Non-standard amino acids, i.e. alpha-amino acids that are not matched by {@link #ALL_STD_AMINOS}.
         * <p>
         * The pattern matches only the amino atom.
         * </p>
         * <p>
         * Refined pattern: the first part describes a generic alpha-amino acid, including alpha,alpha-disubstituted
         * ones like 2-aminoisobutyric acid, as stated in the explanation of the reference implementation ("generic
         * amino acid but not a standard amino acid"). The first part of the original pattern lists the 20 standard
         * side chains, so that only amino acids with a standard side chain but a modified amino or carboxy group are
         * matched, e.g. alanine methyl ester, and no amino acid with a non-standard side chain, e.g. ornithine.
         * </p>
         */
        NON_STD_AMINOS("nonstnd_aminos",
                "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]",
                "[$([NX3,NX4+][CX4;H1,H0]([*])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H2][CX3](=[OX1])[OX2,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3](C)(C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2,OX1-,N])])]"),
        //Sugar derivatives
        /**
         * Open-chain hexose.
         * <p>
         * SMARTS matches a linear six-carbon sugar containing an oxygen atom
         * on every carbon atom.
         * </p>
         * <p>
         * Refined pattern: the carbon-carbon bonds must not be ring bonds ({@code !@}), as the name requires. The
         * original pattern also matches cyclic hexoses, e.g. glucopyranose, so that with the inter-substructure
         * overlap filter, hexopyranoses are counted as open-chain hexoses and the following pyranose key cannot match
         * them.
         * </p>
         */
        S_OPENPYR_C6O6("s_openpyr_C6O6",
                "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
                "C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])"),
        /**
         * Open-chain pentose.
         * <p>
         * SMARTS matches a linear five-carbon sugar containing an oxygen atom on every carbon atom.
         * </p>
         * <p>
         * Refined pattern: the carbon-carbon bonds must not be ring bonds ({@code !@}), as the name and the
         * explanation of the reference implementation ("linear furanose, in open form") require. The original
         * pattern also matches furanoses and hexopyranoses.
         * </p>
         */
        S_OPENFUR_C5O5("s_openfur_C5O5",
                "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])",
                "C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])!@C(~[#8])"),
        /**
         * Pyranose sugar.
         * <p>
         * SMARTS matches a six-membered sugar ring containing one oxygen atom, with oxygen atoms at three
         * neighbouring ring carbon atoms.
         * </p>
         */
        S_PYRANOSE_C5O4("s_pyranose_C5O4", "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1"),
        /**
         * Furanose sugar.
         * <p>
         * SMARTS matches a five-membered sugar ring containing one oxygen atom, with oxygen atoms at two
         * neighbouring ring carbon atoms.
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
         * SMARTS matches a benzene ring connected to a two-carbon chain
         * ending in a nitrogen atom.
         * </p>
         * <p>
         * Refined pattern: the first chain carbon atom must not be part of a five-membered ring, and the chain is
         * written with any-bond expressions. In the original pattern, the chain is written without bond symbols, i.e.
         * with single or aromatic bonds, and may form a five-membered ring with the benzene ring, so that the pyrrole
         * ring of every indole is matched as phenethylamine unit. Isoquinolines, whose phenethylamine unit forms a
         * six-membered ring with the benzene ring, are still matched.
         * </p>
         */
        D_PHENYL_C2N("d_phenylC2N_9",
                "c1ccccc1[#6][#6][#7]",
                "c1ccccc1[#6;!$(*~1~*~*~*~*~1)]~[#6]~[#7]"),
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
         * Phenylpropyl motif, i.e. a C6-C3 unit from the shikimate pathway.
         * <p>
         * SMARTS matches a phenyl ring attached to a three-carbon side chain. The ring atoms must not be part of
         * further rings ({@code R1}), and the second and third side chain carbon atoms must not be part of a
         * six-membered ring.
         * </p>
         * <p>
         * Refined pattern: the ring has to be aromatic, so that cyclohexyl groups and monocyclic terpenes, e.g.
         * menthane derivatives, are not matched; non-aromatic shikimate derivatives like shikimic acid are therefore
         * not matched either. The side chain carbon atoms must not be part of any six-membered ring, which is written
         * as the recursive SMARTS {@code $(*~1~*~*~*~*~*~1)} with any-bond expressions, so that rings with a
         * non-aromatic double bond, e.g. cyclohexene, are recognised as well; without bond symbols, as in
         * {@code $(*1*****1)}, a bond means single or aromatic in SMARTS. The original pattern excludes these atoms
         * via {@code [#6;!$([r6])]} only if their smallest ring is six-membered, so that, e.g., carbon atoms in a
         * six-membered ring fused to the five-membered lactone ring of aryltetralin lignans are not excluded. See the
         * <a href="https://github.com/cdk/cdk/issues/1292">GitHub issue</a> for the meaning of the r primitive.
         * </p>
         */
        D_PHENYL_C3("d_phenylC3_9_strict",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]",
                "[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1~[#6]~[#6;!$(*~1~*~*~*~*~*~1)]~[#6;!$(*~1~*~*~*~*~*~1)]"),
        /**
         * Phenethyl motif, i.e. a C6-C2 unit.
         * <p>
         * SMARTS matches a phenyl ring attached to a two-carbon side chain. The second side chain carbon atom must not
         * be part of a six-membered ring. The refined pattern is described at {@link #D_PHENYL_C3}.
         * </p>
         */
        D_PHENYL_C2("d_phenylC2_8_strict",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]",
                "[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1~[#6]~[#6;!$(*~1~*~*~*~*~*~1)]"),
        /**
         * Benzyl motif, i.e. a C6-C1 unit.
         * <p>
         * SMARTS matches a phenyl ring attached to a single carbon atom. Refined pattern: the ring has to be aromatic,
         * see {@link #D_PHENYL_C3}.
         * </p>
         */
        D_PHENYL_C1("d_phenylC1_7_strict",
                "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]",
                "[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1~[#6]"),
        //Carbon chains
        /**
         * Isoprene unit.
         * <p> SMARTS matches a branched five-carbon motif corresponding to an isoprene
         * building block. Since atomic numbers are used, both aromatic and
         * aliphatic carbon atoms may be matched. According to the explanation of the reference implementation, the
         * pattern deliberately ignores double bonds to detect isoprene units within fused rings.</p>
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
         * <p>
         * Refined pattern: the methyl group must be attached to an aliphatic carbon atom with at least three
         * heavy-atom neighbours that is not a carbonyl carbon atom, i.e. to a branching point as in polyketide chains
         * built from methylmalonyl units; terpenoid methyl groups are still matched, but with the inter-substructure
         * overlap filter, they are taken by the isoprene key first. The original pattern matches every
         * C-C-CH<sub>3</sub> motif, including chain ends, methyl groups on aromatic rings, and methyl ketones, i.e.
         * typically acetate-derived structures.
         * </p>
         */
        D2_METHYLMALONYL("d2_methylmalonyl_C3",
                "[#6]~[#6][C;D1;h3]",
                "[#6]~[C;D3,D4;!$(*=[#8])]-[C;D1;h3]"),
        /**
         * Ethyl group.
         * <p> SMARTS matches two connected carbon atoms, regardless of the bond order. Since atomic
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
         * <p>
         * Refined pattern: only the phosphorus atom is matched, so that every phosphorus atom bonded to at least one
         * oxygen atom is counted once. The original pattern {@code P~O} matches every phosphorus-oxygen pair, so that
         * one phosphate group is counted four times, although the explanation of the reference implementation states
         * that additional oxygen atoms are not matched.
         * </p>
         */
        PHOSPHATE("phosphate_2", "P~O", "[P;$(*~[#8])]"),
        /**
         * Sulfonate group.
         * <p> SMARTS matches a sulfur atom connected to an oxygen atom by any bond, e.g. in sulfoxides, sulfones,
         * sulfonates, and sulfates.
         * </p>
         * <p>
         * Refined pattern: only the sulfur atom is matched, so that every sulfur atom bonded to at least one oxygen
         * atom is counted once. The original pattern {@code S~O} matches every sulfur-oxygen pair, e.g. four times
         * for one sulfate group.
         * </p>
         */
        SULFONATE("sulfonate_2", "S~O", "[S;$(*~[#8])]"),
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
         * (SMARTS degree = 1). This pattern is not restricted to nitrate groups despite
         * the constant name.</p>
         * <p>
         * Refined pattern: the nitrogen atom must be bonded by a single bond, as stated in the explanation of the
         * reference implementation ("N with one single bond"), so that nitriles and imines are not matched.
         * </p>
         */
        N_NITRATE("n_nitrate_1", "[N;D1]", "[N;D1;!$(*=*);!$(*#*)]"),
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
         * <p>
         * Refined pattern: oxygen atoms bonded to a further oxygen atom, i.e. acyclic peroxides, are excluded as well.
         * </p>
         */
        O_ETHER("o_ether_1",
                "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]",
                "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S);!$(*[#8])]"),
        /**
         * Hydroxyl group.
         * <p>
         * SMARTS matches a hydroxyl (-OH) group or a deprotonated hydroxyl group (-O<sup>-</sup>) attached to a carbon
         * or nitrogen atom while excluding carboxylic acids.
         * </p>
         * <p>
         * The original pattern excludes oxygen atoms bonded to a carbon atom that carries a further oxygen atom
         * ({@code !$(*C~O)}), which also excludes the hydroxyl groups of hemiacetals and hemiketals, e.g. the anomeric
         * hydroxyl group of reducing sugars; its parts {@code !$(P)} and {@code !$(S)} have no effect because they test
         * the oxygen atom itself.
         * </p>
         * <p>
         * Refined pattern: only oxygen atoms bonded to a carbonyl carbon atom are excluded ({@code !$(*C=O)}).
         * Additionally, oxygen atoms bonded to a nitrogen atom that carries a further oxygen atom
         * ({@code !$(*[#7]~[#8])}), i.e. the oxygen atoms of nitro groups and nitrate esters, and negatively charged
         * oxygen atoms bonded to a positively charged nitrogen atom ({@code !$([#8-][#7+])}), i.e. the oxygen atoms of
         * N-oxides, are excluded. These groups can be written in a pentavalent form ({@code N(=O)=O}, {@code N=O}), in
         * which the oxygen atom is not matched because it is doubly bonded, or in a charge-separated form
         * ({@code [N+](=O)[O-]}, {@code [N+][O-]}), in which the singly bonded {@code [O-]} would otherwise be matched
         * as deprotonated hydroxyl group, so that the result depended on the representation of the same molecule. In
         * neither form, these oxygen atoms are (deprotonated) hydroxyl groups: the negative charge is a formal charge
         * and does not result from the removal of a proton. The reference implementation matches this {@code [O-]}
         * because RDKit converts nitro groups and N-oxides into the charge-separated form when reading a molecule.
         * Hydroxylamines, oximes, and hydroxamic acids, whose uncharged nitrogen atom carries only one oxygen atom,
         * are still matched.
         * </p>
         */
        O_HYDROXYL("o_hydroxyl_1",
                "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]",
                "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C=O);!$(*[#7]~[#8]);!$([#8-][#7+])]"),
        //rings
        /**
         * Three-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C3("r_c3", "[#6]~1~[#6]~[#6]~1", 3),
        /**
         * Four-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C4("r_c4", "[#6]~1~[#6]~[#6]~[#6]~1", 4),
        /**
         * Five-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C5("r_c5", "[#6]~1~[#6]~[#6]~[#6]~[#6]~1", 5),
        /**
         * Six-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C6("r_c6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1", 6),
        /**
         * Seven-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C7("r_c7", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", 7),
        /**
         * Eight-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C8("r_c8", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", 8),
        /**
         * Nine-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C9("r_c9", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", 9),
        /**
         * Ten-membered carbon ring, regardless of aromaticity. Matched by ring perception in the refined
         * fingerprint, see {@link #getRingSize()}.
         */
        R_C10("r_c10", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1", 10);

        /**
         * Readable name of the fingerprint feature, as used in the reference implementation.
         */
        private final String label;

        /**
         * Original SMARTS pattern of the reference implementation.
         */
        private final String originalSmarts;

        /**
         * Refined SMARTS pattern, identical to the original one if the key has no refinement.
         */
        private final String refinedSmarts;

        /**
         * Number of ring atoms of a ring key, 0 for all other keys.
         */
        private final int ringSize;

        /**
         * Creates a key whose refined SMARTS pattern is identical to the original one.
         *
         * @param label  descriptive name of the fingerprint feature
         * @param smarts SMARTS pattern used to identify the structural motif
         */
        BiosynfoniKey(String label, String smarts) {
            this(label, smarts, smarts, 0);
        }

        /**
         * Creates a key with a refined SMARTS pattern.
         *
         * @param label          descriptive name of the fingerprint feature
         * @param originalSmarts original SMARTS pattern of the reference implementation
         * @param refinedSmarts  refined SMARTS pattern
         */
        BiosynfoniKey(String label, String originalSmarts, String refinedSmarts) {
            this(label, originalSmarts, refinedSmarts, 0);
        }

        /**
         * Creates a ring key, which is matched by its SMARTS pattern in the original fingerprint and by ring
         * perception in the refined fingerprint (see {@link #getRingSize()}).
         *
         * @param label    descriptive name of the fingerprint feature
         * @param smarts   SMARTS pattern of the ring
         * @param ringSize number of ring atoms
         */
        BiosynfoniKey(String label, String smarts, int ringSize) {
            this(label, smarts, smarts, ringSize);
        }

        /**
         * Creates a Biosynfoni fingerprint key definition.
         *
         * <p>Each key consists of a readable label and the SMARTS patterns
         * describing a biosynthetically relevant structural motif. During
         * fingerprint generation, the SMARTS pattern is matched against the
         * molecule to determine the presence or frequency of the corresponding
         * feature.</p>
         *
         * @param label          descriptive name of the fingerprint feature
         * @param originalSmarts original SMARTS pattern of the reference implementation
         * @param refinedSmarts  refined SMARTS pattern
         * @param ringSize       number of ring atoms of a ring key, 0 for all other keys
         */
        BiosynfoniKey(String label, String originalSmarts, String refinedSmarts, int ringSize) {
            this.label = label;
            this.originalSmarts = originalSmarts;
            this.refinedSmarts = refinedSmarts;
            this.ringSize = ringSize;
        }

        /**
         * Returns the readable label of this fingerprint feature, as used in the reference implementation.
         *
         * @return the feature label
         */
        public String getLabel() {
            return this.label;
        }

        /**
         * Returns the original SMARTS pattern of this key, as used in the reference implementation (version
         * {@code full_1103}). It is used if {@link BiosynfoniFingerprinter#isOriginalFingerprint()} is {@code true}.
         *
         * @return the original SMARTS pattern
         */
        public String getOriginalSmarts() {
            return this.originalSmarts;
        }

        /**
         * Returns the refined SMARTS pattern of this key, which is identical to the original one if the key has no
         * refinement. It is used if {@link BiosynfoniFingerprinter#isOriginalFingerprint()} is {@code false} (default),
         * except for the ring keys, which are matched by ring perception then (see {@link #getRingSize()}).
         *
         * @return the refined SMARTS pattern
         */
        public String getRefinedSmarts() {
            return this.refinedSmarts;
        }

        /**
         * Returns the number of ring atoms of a ring key, or 0 for all other keys.
         * <p>
         * A ring SMARTS pattern like {@code [#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1} matches every closed path of the given
         * size, also along the perimeter of fused or bridged ring systems, e.g. the ten carbon atoms of the perimeter
         * of naphthalene, decalin, or steroids, or the six carbon atoms of the perimeter of bicyclo[3.1.0]hexane. In
         * the refined fingerprint, ring keys are therefore matched by ring perception instead of their SMARTS
         * pattern: every relevant cycle (see {@link Cycles#relevant(IAtomContainer)}) of the given size that consists
         * only of carbon atoms is a match. The relevant cycles are the union of all minimum cycle bases, so they do
         * not depend on the input atom order, and they include all faces of symmetric ring systems like cubane.
         * </p>
         *
         * @return number of ring atoms, 0 if the key is not a ring key
         */
        public int getRingSize() {
            return this.ringSize;
        }
    }

    /**
     * Default of whether the original SMARTS patterns of the reference implementation are used instead of the refined
     * ones, see {@link #isOriginalFingerprint()}: {@code false}, i.e. the refined fingerprint is generated by default.
     */
    public static final boolean IS_ORIGINAL_FINGERPRINT_DEFAULT = false;

    /**
     * Whether the original SMARTS patterns of the reference implementation are used ({@code true}) or the refined
     * ones ({@code false}, default), see {@link #isOriginalFingerprint()}. Applies to all keys at once.
     */
    private final boolean isOriginalFingerprint;

    /**
     * Enables filtering of overlapping matches within the same substructure.
     */
    private final boolean intraSubOverlapToggle;

    /**
     * Enables filtering of overlapping matches between different substructures.
     */
    private final boolean interSubOverlapToggle;

    /**
     * Array holding prepared Patterns, indexed by the key position; null for the ring keys in the refined fingerprint,
     * which are matched by ring perception (see {@link BiosynfoniKey#getRingSize()}).
     */
    private final SmartsPattern[] smartsPatterns;

    /**
     * Number of SMARTS patterns.
     */
    private static final int FINGERPRINT_SIZE = BiosynfoniKey.values().length;

    /**
     * Creates a Biosynfoni fingerprinter using the refined SMARTS patterns (see {@link #isOriginalFingerprint()}) with
     * both intra- and inter-substructure overlap filtering disabled. This
     * constructor provides the default configuration of the fingerprinter.
     */
    public BiosynfoniFingerprinter() {
        this(false, false, BiosynfoniFingerprinter.IS_ORIGINAL_FINGERPRINT_DEFAULT);
    }

    /**
     * Creates a Biosynfoni fingerprinter using the refined SMARTS patterns (see {@link #isOriginalFingerprint()}) with
     * the specified overlap filtering configuration.
     *
     * @param intraSubOverlapToggle whether intra-substructure overlap filtering is enabled
     * @param interSubOverLapToggle whether inter-substructure overlap filtering is enabled
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle) {
        this(intraSubOverlapToggle, interSubOverLapToggle, BiosynfoniFingerprinter.IS_ORIGINAL_FINGERPRINT_DEFAULT);
    }

    /**
     * Creates a Biosynfoni fingerprinter with the specified overlap filtering configuration and pattern set.
     *
     * @param intraSubOverlapToggle whether intra-substructure overlap filtering is enabled
     * @param interSubOverLapToggle whether inter-substructure overlap filtering is enabled
     * @param isOriginalFingerprint whether the original SMARTS patterns of the reference implementation are used
     *                              instead of the refined ones (default {@code false}, see
     *                              {@link #IS_ORIGINAL_FINGERPRINT_DEFAULT} and {@link #isOriginalFingerprint()})
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle,
                                   boolean isOriginalFingerprint) {
        this.interSubOverlapToggle = interSubOverLapToggle;
        this.intraSubOverlapToggle = intraSubOverlapToggle;
        this.isOriginalFingerprint = isOriginalFingerprint;

        this.smartsPatterns = new SmartsPattern[BiosynfoniFingerprinter.FINGERPRINT_SIZE];

        int i = 0;
        for (BiosynfoniKey key : BiosynfoniKey.values()) {
            if (isOriginalFingerprint || key.getRingSize() == 0) {
                //the molecule is prepared once for all patterns in getFilteredMatches(), so the automatic
                // preparation of every pattern (ring perception and aromaticity detection) is disabled
                String smarts = isOriginalFingerprint ? key.getOriginalSmarts() : key.getRefinedSmarts();
                this.smartsPatterns[i] = SmartsPattern.create(smarts).setPrepare(false);
            }
            i++;
        }

    }

    /**
     * Returns whether the original SMARTS patterns of the reference implementation (version {@code full_1103}) are
     * used ({@code true}) or the refined ones ({@code false}, default, see {@link #IS_ORIGINAL_FINGERPRINT_DEFAULT}).
     * The setting applies to all keys at once and can only be set in the constructor
     * {@link #BiosynfoniFingerprinter(boolean, boolean, boolean)}. The refinements are described in the class
     * documentation and at the respective keys (see {@link BiosynfoniKey}). The original fingerprint reproduces the
     * reference implementation, apart from the differences listed in the class documentation; the refined fingerprint
     * keeps the 39 positions and their order.
     *
     * @return true if the original SMARTS patterns are used, false if the refined ones are used
     */
    public boolean isOriginalFingerprint() {
        return this.isOriginalFingerprint;
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
     * order. Each SMARTS pattern is then matched against the copy (the ring keys
     * of the refined fingerprint by ring perception, see {@link BiosynfoniKey#getRingSize()}), and the
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
        //ring keys are matched by ring perception in the refined fingerprint, see BiosynfoniKey.getRingSize()
        Map<Integer, List<int[]>> carbocyclesBySize = this.isOriginalFingerprint
                ? null : BiosynfoniFingerprinter.getRelevantCarbocyclesBySize(molecule);
        //atoms of the accepted matches of the previous SMARTS patterns, only used by the inter-substructure overlap filter
        boolean[] blockedByPreviousPatterns = new boolean[molecule.getAtomCount()];
        List<List<int[]>> filteredMatches = new ArrayList<>(BiosynfoniFingerprinter.FINGERPRINT_SIZE);
        BiosynfoniKey[] keys = BiosynfoniKey.values();
        for (int i = 0; i < keys.length; i++) {
            List<int[]> uniqueMatches = this.smartsPatterns[i] == null
                    ? carbocyclesBySize.getOrDefault(keys[i].getRingSize(), Collections.emptyList())
                    : Arrays.asList(this.smartsPatterns[i].matchAll(molecule).uniqueAtoms().toArray());
            List<int[]> subMatches = this.getSubMatches(uniqueMatches, canonicalIndex, blockedByPreviousPatterns);
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
     * Applies the configured overlap filters to the unique atom matches of a key.
     *
     * <p>Depending on the configured overlap options, the matches
     * are filtered for overlaps either within the same SMARTS pattern
     * or with matches identified by previously processed SMARTS patterns.</p>
     *
     * @param uniqueMatches             the unique atom matches of the key, i.e. of its SMARTS pattern or, for the ring
     *                                  keys in the refined fingerprint, the relevant carbocycles of its size; not altered
     * @param canonicalIndex            the canonical atom index mapping used for intra-pattern overlap filtering;
     *                                  null if the intra-substructure overlap filter is disabled
     * @param blockedByPreviousPatterns flags of the atoms of the accepted matches of the previously processed SMARTS
     *                                  patterns, indexed by atom index, used for inter-pattern overlap filtering
     * @return a list of filtered atom matches for the given SMARTS pattern
     */
    private List<int[]> getSubMatches(List<int[]> uniqueMatches, int[] canonicalIndex,
                                      boolean[] blockedByPreviousPatterns) {
        List<int[]> subMatches = new ArrayList<>(uniqueMatches);
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
     * Returns the relevant cycles of the given molecule that consist only of carbon atoms, grouped by their size, as
     * matches of the ring keys in the refined fingerprint (see {@link BiosynfoniKey#getRingSize()}).
     *
     * @param aMolecule the molecule, prepared for SMARTS pattern matching
     * @return atom indices of the relevant carbocycles by ring size
     */
    private static Map<Integer, List<int[]>> getRelevantCarbocyclesBySize(IAtomContainer aMolecule) {
        Map<Integer, List<int[]>> carbocyclesBySize = new HashMap<>();
        for (int[] path : Cycles.relevant(aMolecule).paths()) {
            //the path is closed, i.e. its last atom index repeats the first one
            int[] cycle = Arrays.copyOf(path, path.length - 1);
            boolean isCarbocycle = true;
            for (int atomIndex : cycle) {
                if (aMolecule.getAtom(atomIndex).getAtomicNumber() != IElement.C) {
                    isCarbocycle = false;
                    break;
                }
            }
            if (isCarbocycle) {
                carbocyclesBySize.computeIfAbsent(cycle.length, size -> new ArrayList<>()).add(cycle);
            }
        }
        return carbocyclesBySize;
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
     * SMILES of the molecule, without depending on the SMILES module, and is identical to it for molecules
     * without explicit hydrogen atoms; the Biosynfoni fragmenter of MORTAR uses the same numbering. In contrast to
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