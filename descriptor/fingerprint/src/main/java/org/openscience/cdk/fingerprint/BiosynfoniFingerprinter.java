package org.openscience.cdk.fingerprint;


import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;

/**
 * Because of the overlap filter methods it uses not the Substructure Fingerprint, only orientates at the implementation
 * The Smart expression [#6;!$([r6])] changed to [#6!$(*1*****1)] do to differences in toolkit matching methods,
 * using the changes also in the original implementation will fix all differences containing these SMARTS
 * Loaded molecules containg [S+] in aromatic Rings will match as aliphatic
 * The only diffrences detected so far are
 * <p>
 * Current Features are Fix SMARTS, Count and Bit Fingerprint
 * <p>
 * Overlap filtering:
 * - The {@code intraSubOverlapToggle} (constructor) activates filtering within the same SMARTS
 * so accepted matches are atom-disjoint (the first accepted match blocks its atoms).
 * - The {@code interSubOverLapToggle} (constructor) (experimental) prevents reuse of atoms
 * already accepted by earlier SMARTS (depends on SMARTS order).
 * <p>
 * Preparation and side effects:
 * - Before matching, {@link #preprocessMolecule(IAtomContainer)} is called:
 * - Atom types are perceived and implicit hydrogens are added.
 * - Ring and aromaticity information are computed and set on the passed {@code IAtomContainer}.
 * The method mutates the supplied molecule.
 * - Errors during aromaticity detection or hydrogen addition are wrapped as {@code RuntimeException}.
 * Other matching/fingerprint errors propagate as {@code CDKException} where applicable.
 * <p>
 * Public API (most important):
 * - {@link #BiosynfoniFingerprinter()}:
 * Default constructor (uses default SMARTS, overlap filters off).
 * - {@link #BiosynfoniFingerprinter(boolean, boolean)}:
 * Configure overlap filters.
 * - {@link #BiosynfoniFingerprinter(boolean, boolean, String[])}:
 * Provide a custom SMARTS list; {@link #getSize()} equals {@code smarts.length}.
 * - {@link #getBitFingerprint(IAtomContainer)}:
 * Returns a bit fingerprint; may throw {@code CDKException}.
 * - {@link #getCountFingerprint(IAtomContainer)}:
 * Returns a key-based count fingerprint; SMARTS position = hash position.
 * - {@link #getRawFingerprint(IAtomContainer)}:
 * Not implemented (throws {@code CDKException}).
 * - {@link #getSize()}:
 * Returns the number of SMARTS used.
 *
 */
public class BiosynfoniFingerprinter extends AbstractFingerprinter implements IFingerprinter {
    public void get() {
    }

    public enum DefaultBiosynfoniKey {

        CO_COA("co_coa", "SCCN~C(~O)CCN~C(~O)C(C(C)(C)COP(O)(~O)OP(~O)(O)OCC1C(C(C(O1)[#7]2~[#6]~[#7]~[#6]~3~[#6](~[#7]~[#6]~[#7]~[#6]~3~2)~[#7])O)OP(~O)(O)O)~O"),
        CO_NADH("co_nadh", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O)~O)~O)~O"),// 1: co_nadh
        CO_NADPH("co_nadph", "[#6]~1~[#6]~[#6]~[#7](~[#6]~[#6]~1~[#6](~O)~[#7])~[#6]~2~[#6](~[#6](~[#6](~O~2)~[#6]~O~P(~O)(~O)~O~P(~O)(~O)~O~[#6]~[#6]~3~[#6](~[#6](~[#6](~O~3)~[#7]~4~[#6]~[#7]~[#6]~5~[#6](~[#7]~[#6]~[#7]~[#6]~5~4)~[#7])~O~P(~O)(~O)~O)~O)~O)~O"),

        ALL_STD_AMINOS("allstnd_aminos", "[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]"),
        NON_STD_AMINOS("nonstnd_aminos", "[$([NX3,NX4+][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]"),

        S_OPENPYR_C6O6("s_openpyr_C6O6", "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])"),
        S_OPENFUR_C5O5("s_openfur_C5O5", "C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])~C(~[#8])"),
        S_PYRANOSE_C5O4("s_pyranose_C5O4", "C~1~[#8]~C~C(~[#8])~C(~[#8])~C(~[#8])~1"),
        S_FURANOSE_C4O3("s_furanose_C4O3", "C~1~[#8]~C~C(~[#8])~C(~[#8])~1"),

        D_INDOLE("d_indoleC2N_12", "c1cccc2c1c(~[#6]~[#6]~[#7])cn2"),
        D_PHENYL_C2N("d_phenylC2N_9", "c1ccccc1[#6][#6][#7]"),
        D_C5N("d_c5n_6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#7]~1"),
        D_C4N("d_c4n_5", "[#6]~1~[#6]~[#6]~[#6]~[#7]~1"),


        D_PHENYL_C3("d_phenylC3_9_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]~[#6!$(*1*****1)]"),
        D_PHENYL_C2("d_phenylC2_8_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6!$(*1*****1)]"),
        D_PHENYL_C1("d_phenylC1_7_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]"),
        D_ISOPRENE("d_isoprene_5", "[#6]~[#6](~[#6])~[#6]~[#6]"),

        D2_ACETYL("d2_acetyl_C2O1", "[#6]~[#6]~[#8]"),
        D2_METHYLMALONYL("d2_methylmalonyl_C3", "[#6]~[#6][C;D1;h3]"),
        D_ETHYL("d_ethyl_2", "[#6]~[#6]"),
        D_METHYL("d_methyl_1", "[C;D1;h3]"),

        PHOSPHATE("phosphate_2", "P~O"),
        SULFONATE("sulfonate_2", "S~O"),

        HAL_F("hal_f", "[#9]"),
        HAL_CL("hal_cl", "[#17]"),
        HAL_BR("hal_br", "[#35]"),
        HAL_I("hal_i", "[#53]"),

        N_NITRATE("n_nitrate_1", "[N;D1]"),
        O_EPOXY("o_epoxy_1", "[O;x2;r3]"),
        O_ETHER("o_ether_1", "[O;D2;!h;!$(*C=O);X2;!R;!$(*P);!$(*S)]"),
        O_HYDROXYL("o_hydroxyl_1", "[#8;D1;h,!v2;$(*[#6,#7]);!$(*C~O);!$(P);!$(S)]"),

        R_C3("r_c3", "[#6]~1~[#6]~[#6]~1"),
        R_C4("r_c4", "[#6]~1~[#6]~[#6]~[#6]~1"),
        R_C5("r_c5", "[#6]~1~[#6]~[#6]~[#6]~[#6]~1"),
        R_C6("r_c6", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        R_C7("r_c7", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        R_C8("r_c8", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        R_C9("r_c9", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1"),
        R_C10("r_c10", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1");


        public final String label;
        public final String smarts;

        DefaultBiosynfoniKey(String label, String smarts) {
            this.label = label;
            this.smarts = smarts;
        }
    }

    private boolean intraSubOverlapToggle = false;
    private boolean interSubOverLapToggle = false;
    private String[] smartsList = null;
    private int smartsSize = DefaultBiosynfoniKey.values().length;

    /**
     * uses default SMARTS pattern.
     * Can toggle usage of overlap filter methods
     * on true the overlap filters can be activated
     * See {@link BiosynfoniFingerprinter} for more information
     *
     * @param intraSubOverlapToggle if {@code true}, applies intra-Smarts overlap filtering so that matches within
     *                              the same SMARTS pattern cannot reuse atoms
     * @param interSubOverLapToggle if {@code true}, applies inter-Smarts overlap filtering so that atoms assigned
     *                              to earlier SMARTS patterns cannot be reused
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle) {
        this.intraSubOverlapToggle = intraSubOverlapToggle;
        this.interSubOverLapToggle = interSubOverLapToggle;
    }

    /**
     * uses default creation methods for the Fingerprints \n
     * No overlapping Structures are filtered
     * See {@link BiosynfoniFingerprinter} for more information
     */
    public BiosynfoniFingerprinter() {

    }

    /**
     * uses given SMARTS pattern.
     * Can toggle usage of overlap filter methods
     * on true the overlap filters can be activated
     * See {@link BiosynfoniFingerprinter} for more information
     *
     * @param intraSubOverlapToggle if this is true Fingerprinter will use {@code intraSubOverlap}
     * @param interSubOverLapToggle not yet implemented
     * @param smarts                if this is Null the Fingerprints will be created with the given SMARTS
     */
    public BiosynfoniFingerprinter(boolean intraSubOverlapToggle, boolean interSubOverLapToggle, String[] smarts) {
        this.interSubOverLapToggle = interSubOverLapToggle;
        this.intraSubOverlapToggle = intraSubOverlapToggle;
        this.smartsList = smarts;
        this.smartsSize = smarts.length;
    }


    /**
     * {@inheritDoc}
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        BitSet BIOSYNFingerprintBIT = new BitSet();
        List<List<int[]>> filteredMatches = getFilteredMatches(container);
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
        List<List<int[]>> filteredMatches = getFilteredMatches(container);
        final Map<Integer, Integer> map = new TreeMap<>();
        for (int i = 0; i < filteredMatches.size(); i++) {
            map.put(i, filteredMatches.get(i).size());
        }

        final int size = map.size();
        final int[] hash = new int[size];
        final int[] count = new int[size];
        int n = 0;
        for (int h : map.keySet()) {
            hash[n] = h;
            count[n++] = map.get(h);
        }
        return new ICountFingerprint() {
            @Override
            public long size() {
                return count.length;
            }

            @Override
            public int numOfPopulatedbins() {
                return size;
            }

            @Override
            public int getCount(int index) {
                return count[index];
            }

            /**
             * Note the Fingerprint is Key based.
             * The position of the Smart in the SmartsList equals the position of the hash
             * @param index the index of the bin to return the hash for.
             * @return hash from the given feature index
             */
            @Override
            public int getHash(int index) {
                return hash[index];
            }

            @Override
            public void merge(ICountFingerprint fp) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {

            }

            @Override
            public boolean hasHash(int hash) {
                return map.containsKey(hash);
            }

            @Override
            public int getCountForHash(int hash) {
                return map.containsKey(hash) ? map.get(hash) : 0;
            }
        };
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {
        throw new CDKException("Not yet implemented");
    }

    @Override
    public int getSize() {
        return smartsSize;
    }


    /**
     * This method identifies all occurrences of SMARTS patterns within a given molecule and returns them grouped by pattern.
     * it uses by default the SMARTS given by the enum.
     * Purpose
     * - For each SMARTS (default from {@link DefaultBiosynfoniKey} or the custom {@code smartsList}) collect unique substructure matches.
     * - Optionally apply intra- and inter-overlap filters (via {@link #getSubMatches}) based on the configured toggle flags.
     * If custom Smarts are initialized they are used
     * <p>
     * The result is returned as a list of match lists:
     * - The outer list corresponds to the ordered set of SMARTS keys.
     * This method creates only the outer List the inner values are created via {@link #getSubMatches}
     * -Each inner list contains all matches (int[]) found for the respective SMARTS pattern
     * Description
     * - Prepares the input molecule for SMARTS matching using
     * {@code SmartsPattern.prepare(aMolecule)} and applies aromaticity perception
     * via {@link #preprocessMolecule(IAtomContainer)} once before matching.
     * - Initializes a {@code List<List<int[]>>} to store filtered atom-index matches
     * for each SMARTS pattern.
     * - When {@code smartsList == null}, iterates over
     * {@link DefaultBiosynfoniKey} values and creates a
     * {@code SmartsPattern} from each key's SMARTS definition.
     * - Otherwise, creates a {@code SmartsPattern} for every SMARTS string in
     * {@code smartsList}.
     * -
     *
     * @param aMolecule the molecule to analyze
     * @return a list of match lists grouped by SMARTS pattern
     */
    private List<List<int[]>> getFilteredMatches(IAtomContainer aMolecule) {
        SmartsPattern.prepare(aMolecule);
        List<List<int[]>> filteredMatches = new ArrayList<>(smartsSize);
        // prepare aromaticity and hydrogen's once
        IAtomContainer preparedMol = preprocessMolecule(aMolecule);
        if (smartsList == null) {
            for (DefaultBiosynfoniKey key : DefaultBiosynfoniKey.values()) {

                SmartsPattern pattern = SmartsPattern.create(key.smarts);
                List<int[]> subMatches = getSubMatches(pattern, preparedMol, filteredMatches);
                filteredMatches.add(subMatches);
            }

        } else {
            for (String smarts : smartsList) {
                SmartsPattern pattern = SmartsPattern.create(smarts);
                List<int[]> subMatches = getSubMatches(pattern, preparedMol, filteredMatches);
                filteredMatches.add(subMatches);
            }
        }
        return filteredMatches;
    }

    /**
     * Find unique atom-index matches of a SMARTS pattern in a prepared molecule and apply optional overlap filters.
     * This method identifies all unique atom matches for a SMARTS pattern and applies optional
     * overlap filtering based on the configured toggle flags
     * <p>
     * Description
     * - Uses {@code pattern.matchAll(aMolecule).uniqueAtoms().toArray()} to obtain unique atom-index matches (int[][]).
     * - Collects these matches into a {@code List<int[]>}, then optionally applies:
     * - {@link #intraSubOverlap(List)} when {@code intraSubOverlapToggle} is enabled (makes matches atom-disjoint within the same SMARTS),
     * - {@link #interSubOverlap(List, List)} when {@code interSubOverLapToggle} is enabled (prevents reuse of atoms accepted by earlier SMARTS).
     *
     * @param pattern         a Pattern for matching
     * @param aMolecule       one Molecule for Matching (already prepared via getAromaticity)
     * @param filteredMatches list of previously filtered matches for inter-overlap filtering
     * @return List<int[]> containing  all matches for one Patter with atom indices
     */
    private List<int[]> getSubMatches(SmartsPattern pattern, IAtomContainer aMolecule, List<List<int[]>> filteredMatches) {
        List<int[]> subMatches = new ArrayList<>();
        int[][] uniqueMatches = pattern.matchAll(aMolecule)
                .uniqueAtoms()
                .toArray();

        subMatches.addAll(Arrays.asList(uniqueMatches));
        if (intraSubOverlapToggle) {
            subMatches = intraSubOverlap(subMatches);
        }
        if (interSubOverLapToggle) {
            subMatches = interSubOverlap(subMatches, filteredMatches);
        }
        return subMatches;
    }


    /**
     * Filters overlapping Matches from the same Structure
     * Matches are processed in sorted order. A match is accepted only if none of its
     * atom indices have not already been assigned to a previously accepted match.
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
        List<int[]> filteredSubMatches = new ArrayList<>(smartsSize);

        subMatches.sort((a, b) -> {

            if (a[0] != b[0]) {
                return Integer.compare(a[0], b[0]);
            }
            return Integer.compare(a[1], b[1]);
        });

        Set<Integer> blockedAtoms = new HashSet<>();
        for (int[] aMatch : subMatches) {

            if (!hasOverlap(aMatch, blockedAtoms)) {
                filteredSubMatches.add(aMatch);
                addBlockedAtoms(aMatch, blockedAtoms);
            }

        }
        return filteredSubMatches;
    }


    /**
     * Filter matches to prevent overlaping of structures inside the same pattern
     * Description
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
        List<int[]> filteredSubMatches = new ArrayList<>(smartsSize);
        Set<Integer> blockedAtoms = new HashSet<>();

        for (List<int[]> matches : prevMatches) {
            for (int[] aMatch : matches) {
                addBlockedAtoms(aMatch, blockedAtoms);
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
     *  Checks whether a match overlaps with blocked atoms.
     * - Iterates over the atom indices of the provided match. <p>
     * - Checks whether any atom index is contained in the
     *   {@code blockedAtoms} set. <p>
     * - Returns {@code true} as soon as an overlapping atom is found;
     *   otherwise returns {@code false}. <p>
     *
     * @param aMatch the atom-index match ({@code int[]}) to be checked
     * @param blockedAtoms the set of atom indices already assigned to
     *                     accepted matches
     * @return {@code true} if the match overlaps with blocked atoms,
     *         otherwise {@code false}
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
     *   stored only once. <p>
     *
     * @param match the atom-index match ({@code int[]}) whose atoms should be
     *              added to the blocked set
     * @param blockedAtoms the set of blocked atom indices to be updated
     */
    private void addBlockedAtoms(int[] match, Set<Integer> blockedAtoms) {
        for (int atom : match) {
            blockedAtoms.add(atom);
        }
    }
    /**
     *
     * - Iterates over a list of atom-index matches ({@code int[]}) and collects
     *   all atom indices into a {@code Set<Integer>}. <p>
     * - Ensures that each atom index is stored only once, even if it occurs in
     *   multiple matches. <p>
     * - Returns the set of unique atom indices contained in the provided
     *   matches. <p>
     *
     * @param matches the atom-index matches from which atom indices should be
     *                collected
     * @return a set containing all unique atom indices present in the matches
     */
    private Set<Integer> getBlockedAtoms(List<int[]> matches) {
        Set<Integer> blocked = new HashSet<>();
        for (int[] match : matches) {
            for (int atom : match) {
                blocked.add(atom);
            }
        }
        return blocked;
    }


    /**
     * Detects and assigns aromaticity, ring membership information for the given molecule and adds implicit hydrogen for a molecule
     * The method uses the recomended Ring finding method for pattern matching({@code Cycles.sssr()})
     * {@code AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(IAtomContainer)}
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
     * - Returns the same molecule instance with updated atom typing,
     * hydrogen counts, ring membership, and aromaticity information. <p>
     *
     * @param aMolecule the molecule whose aromaticity should be determined
     * @return the same molecule with updated aromaticity
     */
    private IAtomContainer preprocessMolecule(IAtomContainer aMolecule) {
        try {

           

            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
            CDKHydrogenAdder hydrogenAdder = CDKHydrogenAdder.getInstance(aMolecule.getBuilder());

            hydrogenAdder.addImplicitHydrogens(aMolecule);

            for (IAtom atom : aMolecule.atoms()) {
                atom.setIsAromatic(false);
                atom.setIsInRing(false);
            }
            for (IBond bond : aMolecule.bonds()) {
                bond.setIsAromatic(false);
                bond.setIsInRing(false);
            }
            Cycles cycles = Cycles.sssr(aMolecule);
            IRingSet rings = cycles.toRingSet();

            for (IAtomContainer molecule : rings.atomContainers()) {
                for (IAtom atom : molecule.atoms()) {
                    atom.setIsInRing(true);

                }
                for (IBond bond : molecule.bonds()) {
                    bond.setIsInRing(true);
                }
            }
            Aromaticity aromaticity = new Aromaticity(
                    ElectronDonation.cdk(),
                    Cycles.cdkAromaticSet());

            aromaticity.apply(aMolecule);
            return aMolecule;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}

