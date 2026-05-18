package org.openscience.cdk.fingerprint;


import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;

/**
 * Because of the
 */
public class BiosynfoniFingerprinter extends AbstractFingerprinter implements IFingerprinter {
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


        D_PHENYL_C3("d_phenylC3_9_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]~[#6;!$([r6])]"),
        D_PHENYL_C2("d_phenylC2_8_strict", "[#6;R1]~1~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~[#6;R1]~1~[#6]~[#6;!$([r6])]"),
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


    /**
     * {@inheritDoc}
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        BitSet BIOSYNFingerprintBIT = new BitSet();
        List<List<int[]>> filteredMatches = getFilteredMatchesDefault(container);
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
        List<List<int[]>> filteredMatches = getFilteredMatchesDefault(container);
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
                return DefaultBiosynfoniKey.values().length;
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
        return 0;
    }


    /**
     * This method identifies all occurrences of predefined SMARTS patterns within a given molecule and returns them grouped by pattern.
     * <p>
     * The result is returned as a list of match lists:
     * - The outer list corresponds to the ordered set of SMARTS keys.
     * -Each inner list contains all matches (int[]) found for the respective SMARTS pattern
     *
     * @param aMolecule
     * @return
     */
    private List<List<int[]>> getFilteredMatchesDefault(IAtomContainer aMolecule) {

        List<List<int[]>> filteredMatches = new ArrayList<>(DefaultBiosynfoniKey.values().length);
        Set<Integer> intersubBlockedAtoms = new HashSet<>();
        for (DefaultBiosynfoniKey key : DefaultBiosynfoniKey.values()) {
            List<int[]> subMatches = new ArrayList<>();
            SmartsPattern pattern = SmartsPattern.create(key.smarts);

            int[][] uniqueMatches = pattern.matchAll(getAromaticity(aMolecule))
                    .uniqueAtoms()
                    .toArray();

            for (int[] match : uniqueMatches) {
                subMatches.add(match);
            }

            filteredMatches.add(subMatches);
        }
        return filteredMatches;
    }
    /*
     * #To
     */
    private List<int[]> intrasubOverlap(List<int[]> subMatches) {
        List<int[]> filteredMatches = new ArrayList<>(DefaultBiosynfoniKey.values().length);

        subMatches.sort(Collections.reverseOrder());
        Set<Integer> blockedAtoms = new HashSet<>();
        for (int[]  aMatch : subMatches) {

            if(!hasOverlap(aMatch,blockedAtoms)){
                filteredMatches.add(aMatch);
                addAllBlockedAtoms(aMatch, blockedAtoms);
            }

        }
        return filteredMatches;
    }

    private List<int[]> filterIntraSmarter(List<int[]> subMatches){
        List<int[]> filterdMatches = new ArrayList<>(DefaultBiosynfoniKey.values().length);

        //subMatches.sort();
        return filterdMatches;
    }

    private boolean hasOverlap(int[] aMatch, Set<Integer> blockedAtoms) {
        for (int Atom : aMatch) {
            if (blockedAtoms.contains(Atom)) {
                return true;
            }
        }
        return false;
    }

    private void addAllBlockedAtoms(int[] match, Set<Integer> blockedAtoms) {
        for (int Atom : match) {
            blockedAtoms.add(Atom);
        }
    }


    /**
     * Detects and assigins aromaticity information for the given molecule
     * <p>
     * This method adds implicit hydrogens, resets aromaticity flags on bonds and atoms,
     * Applies CDK aromaticity detection
     *
     * @param aMolecule the molcule whose aromaticity should be determined
     * @return the same molecule with updated aromaticity
     */
    public IAtomContainer getAromaticity(IAtomContainer aMolecule) {
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
            IRingSet rings = Cycles.relevant().find(aMolecule).toRingSet();

            for(IAtomContainer molecule : rings.atomContainers()) {
                for(IAtom atom : molecule.atoms()) {
                    atom.setIsInRing(true);

                }
                for(IBond bond : molecule.bonds()) {
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