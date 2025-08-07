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

import org.openscience.cdk.Bond;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IElement;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.Mappings;
import org.openscience.cdk.isomorphism.Transform;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smirks.Smirks;
import org.openscience.cdk.smirks.SmirksTransform;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * Utility class for detecting and extracting sugar moieties from molecular structures.
 *
 * <p>This class extends {@link SugarRemovalUtility} to provide functionality for separating
 * glycosides into their aglycone and sugar components.
 * The main feature is the ability to create copies of both the aglycone (non-sugar backbone)
 * and individual sugar fragments from a given molecule, with proper handling of attachment
 * points and stereochemistry.
 *
 * <p>The extraction process supports:
 * <ul>
 *   <li>Detection and extraction of both circular and linear sugar moieties</li>
 *   <li>Preservation of stereochemistry at connection points</li>
 *   <li>Proper saturation of broken bonds with either R-groups or implicit hydrogens</li>
 *   <li>Post-processing of sugar fragments including bond splitting (O-glycosidic, ether, ester, peroxide)</li>
 *   <li>Handling of connecting heteroatoms (oxygen, nitrogen, sulfur) in glycosidic bonds</li>
 * </ul>
 *
 * <p>All sugar detection and removal operations respect the settings inherited from the
 * parent {@link SugarRemovalUtility} class, including terminal vs. non-terminal sugar
 * removal, preservation mode settings, and various detection thresholds.
 *
 * <p><strong>Usage Example:</strong>
 * <pre>{@code
 * SugarDetectionUtility utility = new SugarDetectionUtility(SilentChemObjectBuilder.getInstance());
 * List<IAtomContainer> fragments = utility.copyAndExtractAglyconeAndSugars(molecule, true, false, false, false);
 * IAtomContainer aglycone = fragments.get(0);  // First element is always the aglycone
 * // Subsequent elements are individual sugar fragments
 * }</pre>
 *
 * @author Jonas Schaub
 */
public class SugarDetectionUtility extends SugarRemovalUtility {

    /**
     * Logger of this class.
     */
    private static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(SugarDetectionUtility.class);

    /**
     * Sole constructor of this class. All settings are set to their default
     * values as declared in the {@link SugarRemovalUtility} class.
     *
     * @param builder IChemObjectBuilder for i.a. parsing SMILES strings of
     *                sugar patterns into atom containers
     */
    public SugarDetectionUtility(IChemObjectBuilder builder) {
        super(builder);
    }

    /**
     * Extracts copies of the aglycone and (specified) sugar parts of the given molecule (if there are any).
     * <p>
     * This method creates a deep copy of the input molecule and removes the specified
     * sugar moieties (circular and/or linear) to produce an aglycone. It then creates
     * a second copy to extract the sugar fragments that were removed. The attachment
     * points between the aglycone and sugars are handled by either adding R-groups
     * (pseudo atoms) or implicit hydrogens to saturate the broken bonds.
     *
     * <p>The method preserves stereochemistry information at connection points and
     * handles glycosidic bonds appropriately. When bonds are broken between sugar
     * moieties and the aglycone, connecting heteroatoms (such as glycosidic oxygen,
     * nitrogen, or sulfur atoms) are copied to both the aglycone and sugar fragments
     * to maintain chemical validity.
     *
     * <p>The extraction process respects all current sugar detection settings as described in
     * {@link SugarRemovalUtility}, including terminal vs. non-terminal sugar removal,
     * preservation mode settings, and various detection thresholds.
     *
     * <p>Note that atom types are not copied, they have to be re-perceived if needed.</p>
     *
     * @param mol The input molecule to separate into aglycone and sugar components.
     *            Must not be null but can be empty; a list containing only the empty given
     *            atom container is returned in the latter case.
     * @return A list of atom containers where the first element is the aglycone
     *         (copy molecule with sugars removed) and subsequent elements are the
     *         individual sugar fragments that were extracted (also copies). If no sugars were
     *         detected or removed, returns a list containing only a copy of the
     *         original molecule. Sugar fragments may be disconnected from each
     *         other if they were not directly linked in the original structure.
     * @throws NullPointerException if the input molecule is null
     * @see SugarRemovalUtility#removeCircularSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeLinearSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeCircularAndLinearSugars(IAtomContainer)
     */
    public List<IAtomContainer> copyAndExtractAglyconeAndSugars(
            IAtomContainer mol) {
        return this.copyAndExtractAglyconeAndSugars(
                mol, true, false, false, false, null, null, null, null);
    }

    /**
     * Extracts copies of the aglycone and (specified) sugar parts of the given molecule (if there are any).
     * <p>
     * This method creates a deep copy of the input molecule and removes the specified
     * sugar moieties (circular and/or linear) to produce an aglycone. It then creates
     * a second copy to extract the sugar fragments that were removed. The attachment
     * points between the aglycone and sugars are handled by either adding R-groups
     * (pseudo atoms) or implicit hydrogens to saturate the broken bonds.
     *
     * <p>The method preserves stereochemistry information at connection points and
     * handles glycosidic bonds appropriately. When bonds are broken between sugar
     * moieties and the aglycone, connecting heteroatoms (such as glycosidic oxygen,
     * nitrogen, or sulfur atoms) are copied to both the aglycone and sugar fragments
     * to maintain chemical validity.
     *
     * <p>The extraction process respects all current sugar detection settings as described in
     * {@link SugarRemovalUtility}, including terminal vs. non-terminal sugar removal,
     * preservation mode settings, and various detection thresholds.
     *
     * <p>Note that atom types are not copied, they have to be re-perceived if needed.</p>
     *
     * @param mol The input molecule to separate into aglycone and sugar components.
     *            Must not be null but can be empty; a list containing only the empty given
     *            atom container is returned in the latter case.
     * @param extractCircularSugars If true, circular sugar moieties will be detected
     *                             and extracted according to current settings.
     * @param extractLinearSugars If true, linear sugar moieties will be detected
     *                           and extracted according to current settings.
     * @return A list of atom containers where the first element is the aglycone
     *         (copy molecule with sugars removed) and subsequent elements are the
     *         individual sugar fragments that were extracted (also copies). If no sugars were
     *         detected or removed, returns a list containing only a copy of the
     *         original molecule. Sugar fragments may be disconnected from each
     *         other if they were not directly linked in the original structure.
     * @throws NullPointerException if the input molecule is null
     * @see SugarRemovalUtility#removeCircularSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeLinearSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeCircularAndLinearSugars(IAtomContainer)
     */
    public List<IAtomContainer> copyAndExtractAglyconeAndSugars(
            IAtomContainer mol,
            boolean extractCircularSugars,
            boolean extractLinearSugars) {
        return this.copyAndExtractAglyconeAndSugars(
                mol, extractCircularSugars, extractLinearSugars, false, false, null, null, null, null);
    }

    /**
     * Extracts copies of the aglycone and (specified) sugar parts of the given molecule (if there are any).
     * <p>
     * This method creates a deep copy of the input molecule and removes the specified
     * sugar moieties (circular and/or linear) to produce an aglycone. It then creates
     * a second copy to extract the sugar fragments that were removed. The attachment
     * points between the aglycone and sugars are handled by either adding R-groups
     * (pseudo atoms) or implicit hydrogens to saturate the broken bonds.
     *
     * <p>The method preserves stereochemistry information at connection points and
     * handles glycosidic bonds appropriately. When bonds are broken between sugar
     * moieties and the aglycone, connecting heteroatoms (such as glycosidic oxygen,
     * nitrogen, or sulfur atoms) are copied to both the aglycone and sugar fragments
     * to maintain chemical validity.
     *
     * <p>The extraction process respects all current sugar detection settings as described in
     * {@link SugarRemovalUtility}, including terminal vs. non-terminal sugar removal,
     * preservation mode settings, and various detection thresholds.
     *
     * <p>Note that atom types are not copied, they have to be re-perceived if needed.</p>
     *
     * @param mol The input molecule to separate into aglycone and sugar components.
     *            Must not be null but can be empty; a list containing only the empty given
     *            atom container is returned in the latter case.
     * @param markAttachPointsByR If true, attachment points where sugars and the aglycone were connected
     *                           are marked with R-groups (pseudo atoms). If false,
     *                           implicit hydrogens are added to saturate the connections.
     * @param extractCircularSugars If true, circular sugar moieties will be detected
     *                             and extracted according to current settings.
     * @param extractLinearSugars If true, linear sugar moieties will be detected
     *                           and extracted according to current settings.
     * @return A list of atom containers where the first element is the aglycone
     *         (copy molecule with sugars removed) and subsequent elements are the
     *         individual sugar fragments that were extracted (also copies). If no sugars were
     *         detected or removed, returns a list containing only a copy of the
     *         original molecule. Sugar fragments may be disconnected from each
     *         other if they were not directly linked in the original structure.
     * @throws NullPointerException if the input molecule is null
     * @see SugarRemovalUtility#removeCircularSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeLinearSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeCircularAndLinearSugars(IAtomContainer)
     */
    public List<IAtomContainer> copyAndExtractAglyconeAndSugars(
            IAtomContainer mol,
            boolean extractCircularSugars,
            boolean extractLinearSugars,
            boolean markAttachPointsByR) {
        return this.copyAndExtractAglyconeAndSugars(
                mol, extractCircularSugars, extractLinearSugars, markAttachPointsByR, false, null, null, null, null);
    }

    /**
     * Extracts copies of the aglycone and (specified) sugar parts of the given molecule (if there are any).
     * <p>
     * This method creates a deep copy of the input molecule and removes the specified
     * sugar moieties (circular and/or linear) to produce an aglycone. It then creates
     * a second copy to extract the sugar fragments that were removed. The attachment
     * points between the aglycone and sugars are handled by either adding R-groups
     * (pseudo atoms) or implicit hydrogens to saturate the broken bonds.
     *
     * <p>The method preserves stereochemistry information at connection points and
     * handles glycosidic bonds appropriately. When bonds are broken between sugar
     * moieties and the aglycone, connecting heteroatoms (such as glycosidic oxygen,
     * nitrogen, or sulfur atoms) are copied to both the aglycone and sugar fragments
     * to maintain chemical validity.
     *
     * <p>The extraction process respects all current sugar detection settings as described in
     * {@link SugarRemovalUtility}, including terminal vs. non-terminal sugar removal,
     * preservation mode settings, and various detection thresholds.
     *
     * <p>Note that atom types are not copied, they have to be re-perceived if needed.</p>
     *
     * @param mol The input molecule to separate into aglycone and sugar components.
     *            Must not be null but can be empty; a list containing only the empty given
     *            atom container is returned in the latter case.
     * @param extractCircularSugars If true, circular sugar moieties will be detected
     *                             and extracted according to current settings.
     * @param extractLinearSugars If true, linear sugar moieties will be detected
     *                           and extracted according to current settings.
     * @param markAttachPointsByR If true, attachment points where sugars and the aglycone were connected
     *                           are marked with R-groups (pseudo atoms). If false,
     *                           implicit hydrogens are added to saturate the connections.
     * @param postProcessSugars If true, postprocessing of sugar fragments is performed, i.e. splitting O-glycosidic
     *                          bonds in circular and splitting ether, ester, and peroxide bonds in linear sugar moieties
     * @return A list of atom containers where the first element is the aglycone
     *         (copy molecule with sugars removed) and subsequent elements are the
     *         individual sugar fragments that were extracted (also copies). If no sugars were
     *         detected or removed, returns a list containing only a copy of the
     *         original molecule. Sugar fragments may be disconnected from each
     *         other if they were not directly linked in the original structure or were disconnected in postprocessing.
     * @throws NullPointerException if the input molecule is null
     * @see SugarRemovalUtility#removeCircularSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeLinearSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeCircularAndLinearSugars(IAtomContainer)
     */
    public List<IAtomContainer> copyAndExtractAglyconeAndSugars(
            IAtomContainer mol,
            boolean extractCircularSugars,
            boolean extractLinearSugars,
            boolean markAttachPointsByR,
            boolean postProcessSugars) {
        return this.copyAndExtractAglyconeAndSugars(
                mol, extractCircularSugars, extractLinearSugars, markAttachPointsByR, postProcessSugars,
                null, null, null, null);
    }

    //TODO: the postprocessing of sugars destroys the mapping, do the SMIRKSTransformations generate completely new atoms?
    //remove this method after all and incorporate the new behaviour into the existing methods? -> better to keep the original behaviour of the original methods
    //do not copy the aglycone? -> too much of a hassle because for postprocessing, we repeatedly need the original structure
    //TODO: simplify this method by encapsulating more code
    //TODO: add special treatment for esters (on the sugar side and on the aglycone side, respectively)?
    //TODO: look at other special cases in the test class that might require additional postprocessing
    //TODO: check doc of all overloaded methods and ensure that they are consistent
    /**
     * Extracts copies of the aglycone and (specified) sugar parts of the given molecule (if there are any).
     * <p>
     * This method creates a deep copy of the input molecule and removes the specified
     * sugar moieties (circular and/or linear) to produce an aglycone. It then creates
     * a second copy to extract the sugar fragments that were removed. The attachment
     * points between the aglycone and sugars are handled by either adding R-groups
     * (pseudo atoms) or implicit hydrogens to saturate the broken bonds.
     *
     * <p>The method preserves stereochemistry information at connection points and
     * handles glycosidic bonds appropriately. When bonds are broken between sugar
     * moieties and the aglycone, connecting heteroatoms (such as glycosidic oxygen,
     * nitrogen, or sulfur atoms) are copied to both the aglycone and sugar fragments
     * to maintain chemical validity.
     *
     * <p>The extraction process respects all current sugar detection settings as described in
     * {@link SugarRemovalUtility}, including terminal vs. non-terminal sugar removal,
     * preservation mode settings, and various detection thresholds.
     *
     * <p>Note that atom types are not copied, they have to be re-perceived if needed.</p>
     *
     * <p>This method additionally gives you the option to supply four maps as parameters that will be filled with a
     * mapping of atoms and bonds in the original molecule to the atoms and bonds in the aglycone and sugar copies. They
     * should be of sufficient size and empty when given.</p>
     *
     * @param mol The input molecule to separate into aglycone and sugar components.
     *            Must not be null but can be empty; a list containing only the empty given
     *            atom container is returned in the latter case.
     * @param extractCircularSugars If true, circular sugar moieties will be detected
     *                             and extracted according to current settings.
     * @param extractLinearSugars If true, linear sugar moieties will be detected
     *                           and extracted according to current settings.
     * @param markAttachPointsByR If true, attachment points where sugars and the aglycone were connected
     *                           are marked with R-groups (pseudo atoms). If false,
     *                           implicit hydrogens are added to saturate the connections.
     * @param postProcessSugars If true, postprocessing of sugar fragments is performed, i.e. splitting O-glycosidic
     *                          bonds in circular and splitting ether, ester, and peroxide bonds in linear sugar moieties
     * @return A list of atom containers where the first element is the aglycone
     *         (copy molecule with sugars removed) and subsequent elements are the
     *         individual sugar fragments that were extracted (also copies). If no sugars were
     *         detected or removed, returns a list containing only a copy of the
     *         original molecule. Sugar fragments may be disconnected from each
     *         other if they were not directly linked in the original structure or were disconnected in postprocessing.
     * @throws NullPointerException if the input molecule is null
     * @see SugarRemovalUtility#removeCircularSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeLinearSugars(IAtomContainer)
     * @see SugarRemovalUtility#removeCircularAndLinearSugars(IAtomContainer)
     */
    public List<IAtomContainer> copyAndExtractAglyconeAndSugars(
            IAtomContainer mol,
            boolean extractCircularSugars,
            boolean extractLinearSugars,
            boolean markAttachPointsByR,
            boolean postProcessSugars,
            Map<IAtom, IAtom> inputAtomToAtomCopyInAglyconeMap,
            Map<IBond, IBond> inputBondToBondCopyInAglyconeMap,
            Map<IAtom, IAtom> inputAtomToAtomCopyInSugarsMap,
            Map<IBond, IBond> inputBondToBondCopyInSugarsMap) {
        //checks
        if (mol == null) {
            throw new NullPointerException("Given molecule is null.");
        }
        if (mol.isEmpty()) {
            List<IAtomContainer> results = new ArrayList<>(1);
            results.add(mol);
            return results;
        }
        //setup and copying for aglycone
        float loadFactor = 0.75f; //default load factor for HashMaps
        //ensuring sufficient initial capacity
        int atomMapInitCapacity = (int)((mol.getAtomCount() / loadFactor) + 3.0f);
        int bondMapInitCapacity = (int)((mol.getBondCount() / loadFactor) + 3.0f);
        if (inputAtomToAtomCopyInAglyconeMap == null) {
            inputAtomToAtomCopyInAglyconeMap = new HashMap<>(atomMapInitCapacity);
        }
        if (inputBondToBondCopyInAglyconeMap == null) {
            inputBondToBondCopyInAglyconeMap = new HashMap<>(bondMapInitCapacity);
        }
        IAtomContainer copyForAglycone = this.deeperCopy(mol, inputAtomToAtomCopyInAglyconeMap, inputBondToBondCopyInAglyconeMap);
        boolean wasSugarRemoved = false;
        if (extractCircularSugars && extractLinearSugars) {
            wasSugarRemoved = this.removeCircularAndLinearSugars(copyForAglycone);
        } else if (extractCircularSugars) {
            wasSugarRemoved = this.removeCircularSugars(copyForAglycone);
        } else if (extractLinearSugars) {
            wasSugarRemoved = this.removeLinearSugars(copyForAglycone);
        } //else: wasSugarRemoved remains false, and input structure is returned, same as when no sugars were detected, see below
        if (!wasSugarRemoved) {
            List<IAtomContainer> results = new ArrayList<>(1);
            results.add(copyForAglycone);
            return results;
        }
        //copying for sugars
        if (inputAtomToAtomCopyInSugarsMap == null) {
            inputAtomToAtomCopyInSugarsMap = new HashMap<>(atomMapInitCapacity);
        }
        if (inputBondToBondCopyInSugarsMap == null) {
            inputBondToBondCopyInSugarsMap = new HashMap<>(bondMapInitCapacity);
        }
        IAtomContainer copyForSugars = this.deeperCopy(mol, inputAtomToAtomCopyInSugarsMap, inputBondToBondCopyInSugarsMap);
        //remove aglycone atoms from sugar container
        //note: instead of copying the whole structure and removing the aglycone atoms, one could only copy those atoms
        // and bonds that are not part of the aglycone to form the sugars to save some memory but the code would be much
        // more complicated, so we don't do it that way for now
        for (IAtom atom : mol.atoms()) {
            if (copyForAglycone.contains(inputAtomToAtomCopyInAglyconeMap.get(atom))) {
                copyForSugars.removeAtom(inputAtomToAtomCopyInSugarsMap.get(atom));
            }
        }
        //note that the four atom and bond maps still hold references to the atoms and bonds that were removed from the
        // two copies to get the aglycone and sugars; important for later queries; only cleared at the end of this method
        boolean hasIdentifiedBrokenBond = false;
        //identify bonds that were broken between sugar moieties and aglycone
        // -> copy connecting hetero atoms (glycosidic O/N/S etc.) from one part (sugar or aglycone) to the other,
        // along with its stereo element
        // -> saturate with R or H, depending on the markAttachPointsByR parameter
        for (IBond bond : mol.bonds()) {
            //bond not in aglycone or sugars, so it was broken during sugar removal
            if (!copyForAglycone.contains(inputBondToBondCopyInAglyconeMap.get(bond))
                    && !copyForSugars.contains(inputBondToBondCopyInSugarsMap.get(bond))) {
                hasIdentifiedBrokenBond = true;
                //if one hetero atom connected sugar and aglycone, copy it to the other side;
                // do nothing for C-C bonds or hetero-hetero connections like peroxides
                if ((this.isCarbonAtom(bond.getBegin()) && this.isHeteroAtom(bond.getEnd()))
                        || (this.isHeteroAtom(bond.getBegin()) && this.isCarbonAtom(bond.getEnd()))) {
                    //-> copy hetero atom to the other side and saturate it with H or R
                    //-> saturate "original" hetero atom with H or R
                    IAtom origHeteroAtom;
                    IAtom origCarbonAtom;
                    if (this.isCarbonAtom(bond.getBegin()) && this.isHeteroAtom(bond.getEnd())) {
                        origHeteroAtom = bond.getEnd();
                        origCarbonAtom = bond.getBegin();
                    } else if (this.isCarbonAtom(bond.getEnd()) && this.isHeteroAtom(bond.getBegin())) {
                        origHeteroAtom = bond.getBegin();
                        origCarbonAtom = bond.getEnd();
                    } else {
                        SugarDetectionUtility.LOGGER.error("Broken bond between sugar and aglycone with one carbon " +
                                "and one hetero atom found but they cannot be assigned, this should not happen!");
                        continue;
                    }
                    boolean isHeteroAtomInAglycone = copyForAglycone.contains(inputAtomToAtomCopyInAglyconeMap.get(origHeteroAtom));
                    boolean isHeteroAtomInSugars = copyForSugars.contains(inputAtomToAtomCopyInSugarsMap.get(origHeteroAtom));
                    if (!(isHeteroAtomInAglycone || isHeteroAtomInSugars)) {
                        SugarDetectionUtility.LOGGER.error("Hetero atom not found in aglycone or sugars, this should not happen!");
                        continue;
                    }
                    //copy hetero atom to the other part
                    IAtom cpyHeteroAtom = this.deeperCopy(origHeteroAtom,
                            isHeteroAtomInSugars? copyForAglycone : copyForSugars);
                    IBond copyBondToHeteroAtom;
                    IAtom carbonAtomToBindTo = isHeteroAtomInSugars?
                            inputAtomToAtomCopyInAglyconeMap.get(origCarbonAtom) : inputAtomToAtomCopyInSugarsMap.get(origCarbonAtom);
                    if (bond.getBegin().equals(origCarbonAtom)) {
                        copyBondToHeteroAtom = carbonAtomToBindTo.getBuilder().newInstance(
                                IBond.class, carbonAtomToBindTo, cpyHeteroAtom, bond.getOrder());
                    } else {
                        copyBondToHeteroAtom = carbonAtomToBindTo.getBuilder().newInstance(
                                IBond.class, cpyHeteroAtom, carbonAtomToBindTo, bond.getOrder());
                    }
                    if (isHeteroAtomInSugars) {
                        copyForAglycone.addBond(copyBondToHeteroAtom);
                        inputAtomToAtomCopyInAglyconeMap.put(origHeteroAtom, cpyHeteroAtom);
                        inputBondToBondCopyInAglyconeMap.put(bond, copyBondToHeteroAtom);
                    } else {
                        copyForSugars.addBond(copyBondToHeteroAtom);
                        inputAtomToAtomCopyInSugarsMap.put(origHeteroAtom, cpyHeteroAtom);
                        inputBondToBondCopyInSugarsMap.put(bond, copyBondToHeteroAtom);
                    }
                    //saturate copied hetero atom with H or R
                    if (markAttachPointsByR) {
                        IPseudoAtom tmpRAtom = cpyHeteroAtom.getBuilder().newInstance(IPseudoAtom.class, "R");
                        tmpRAtom.setAttachPointNum(1);
                        tmpRAtom.setImplicitHydrogenCount(0);
                        IBond bondToR = cpyHeteroAtom.getBuilder().newInstance(
                                IBond.class, cpyHeteroAtom, tmpRAtom, bond.getOrder());
                        cpyHeteroAtom.setImplicitHydrogenCount((int) (cpyHeteroAtom.getImplicitHydrogenCount()
                                + mol.getBondOrderSum(origHeteroAtom) - (1 + bond.getOrder().numeric())));
                        if (isHeteroAtomInSugars) {
                            copyForAglycone.addAtom(tmpRAtom);
                            copyForAglycone.addBond(bondToR);
                        } else {
                            copyForSugars.addAtom(tmpRAtom);
                            copyForSugars.addBond(bondToR);
                        }
                    } else {
                        cpyHeteroAtom.setImplicitHydrogenCount((int) (cpyHeteroAtom.getImplicitHydrogenCount()
                                + mol.getBondOrderSum(origHeteroAtom) - bond.getOrder().numeric()));
                    }
                    //copy stereo elements for the broken bond to preserve the configuration
                    IAtomContainer receivingPart = isHeteroAtomInSugars? copyForAglycone : copyForSugars;
                    Map<IAtom, IAtom> receivingPartOrigAtomToCpy = isHeteroAtomInSugars? inputAtomToAtomCopyInAglyconeMap : inputAtomToAtomCopyInSugarsMap;
                    Map<IBond, IBond> receivingPartOrigBondToCpy = isHeteroAtomInSugars? inputBondToBondCopyInAglyconeMap : inputBondToBondCopyInSugarsMap;
                    for (IStereoElement elem : mol.stereoElements()) {
                        if (elem.contains(bond.getBegin()) && elem.contains(bond.getEnd())
                                && receivingPart.contains(receivingPartOrigAtomToCpy.get(elem.getFocus()))) {
                            boolean carriersAllPresent = true;
                            for (Object object : elem.getCarriers()) {
                                if (object instanceof IAtom) {
                                    if (!receivingPart.contains(receivingPartOrigAtomToCpy.get(object))) {
                                        carriersAllPresent = false;
                                        break;
                                    }
                                } else if (object instanceof IBond) {
                                    if (!receivingPart.contains(receivingPartOrigBondToCpy.get(object))) {
                                        carriersAllPresent = false;
                                        break;
                                    }
                                } else {
                                    carriersAllPresent = false;
                                    break;
                                }
                            }
                            if (carriersAllPresent) {
                                receivingPart.addStereoElement(elem.map(receivingPartOrigAtomToCpy, receivingPartOrigBondToCpy));
                            }
                        }
                    }
                    //saturate the hetero atom in the respective part with H or R
                    if (markAttachPointsByR) {
                        IPseudoAtom tmpRAtom = origHeteroAtom.getBuilder().newInstance(IPseudoAtom.class, "R");
                        tmpRAtom.setAttachPointNum(1);
                        tmpRAtom.setImplicitHydrogenCount(0);
                        IBond tmpNewBond;
                        IAtom partHeteroAtom = isHeteroAtomInAglycone?
                                inputAtomToAtomCopyInAglyconeMap.get(origHeteroAtom) : inputAtomToAtomCopyInSugarsMap.get(origHeteroAtom);
                        if (bond.getBegin().equals(origHeteroAtom)) {
                            tmpNewBond = origHeteroAtom.getBuilder().newInstance(IBond.class, partHeteroAtom, tmpRAtom, bond.getOrder());
                        } else {
                            tmpNewBond = origHeteroAtom.getBuilder().newInstance(IBond.class, tmpRAtom, partHeteroAtom, bond.getOrder());
                        }
                        if (isHeteroAtomInAglycone) {
                            copyForAglycone.addAtom(tmpRAtom);
                            copyForAglycone.addBond(tmpNewBond);
                        } else {
                            copyForSugars.addAtom(tmpRAtom);
                            copyForSugars.addBond(tmpNewBond);
                        }
                    } else {
                        if (isHeteroAtomInAglycone) {
                            IAtom bondAtomInAglycone = inputAtomToAtomCopyInAglyconeMap.get(origHeteroAtom);
                            int implHCount = bondAtomInAglycone.getImplicitHydrogenCount();
                            bondAtomInAglycone.setImplicitHydrogenCount(implHCount + bond.getOrder().numeric());
                        } else {
                            IAtom bondAtomInSugars = inputAtomToAtomCopyInSugarsMap.get(origHeteroAtom);
                            int implHCount = bondAtomInSugars.getImplicitHydrogenCount();
                            bondAtomInSugars.setImplicitHydrogenCount(implHCount + bond.getOrder().numeric());
                        }
                    }
                } else if (markAttachPointsByR) {
                    //broken bond was a C-C or hetero-hetero bond, just saturate both former bond atoms with R if required
                    for (IAtom atom : bond.atoms()) {
                        IPseudoAtom tmpRAtom = atom.getBuilder().newInstance(IPseudoAtom.class, "R");
                        tmpRAtom.setAttachPointNum(1);
                        tmpRAtom.setImplicitHydrogenCount(0);
                        if (copyForAglycone.contains(inputAtomToAtomCopyInAglyconeMap.get(atom))) {
                            copyForAglycone.addAtom(tmpRAtom);
                            IBond bondToR = atom.getBuilder().newInstance(
                                    IBond.class, inputAtomToAtomCopyInAglyconeMap.get(atom), tmpRAtom, bond.getOrder());
                            copyForAglycone.addBond(bondToR);
                            inputAtomToAtomCopyInAglyconeMap.get(atom).setImplicitHydrogenCount(
                                    inputAtomToAtomCopyInAglyconeMap.get(atom).getImplicitHydrogenCount()
                                            + bond.getOrder().numeric() - 1);
                        } else if (copyForSugars.contains(inputAtomToAtomCopyInSugarsMap.get(atom))) {
                            copyForSugars.addAtom(tmpRAtom);
                            IBond bondToR = atom.getBuilder().newInstance(
                                    IBond.class, inputAtomToAtomCopyInSugarsMap.get(atom), tmpRAtom, bond.getOrder());
                            copyForSugars.addBond(bondToR);
                            inputAtomToAtomCopyInSugarsMap.get(atom).setImplicitHydrogenCount(
                                    inputAtomToAtomCopyInSugarsMap.get(atom).getImplicitHydrogenCount()
                                            + bond.getOrder().numeric() - 1);
                        } else {
                            SugarDetectionUtility.LOGGER.error("Bond atom neither found in aglycone nor in sugars, this should not happen!");
                        }
                    }
                } else {
                    //saturate both former bond atoms with implicit hydrogens
                    for (IAtom atom : bond.atoms()) {
                        if (copyForAglycone.contains(inputAtomToAtomCopyInAglyconeMap.get(atom))) {
                            IAtom bondAtomInAglycone = inputAtomToAtomCopyInAglyconeMap.get(atom);
                            int implHCount = bondAtomInAglycone.getImplicitHydrogenCount();
                            bondAtomInAglycone.setImplicitHydrogenCount(implHCount + bond.getOrder().numeric());
                        } else if (copyForSugars.contains(inputAtomToAtomCopyInSugarsMap.get(atom))) {
                            IAtom bondAtomInSugars = inputAtomToAtomCopyInSugarsMap.get(atom);
                            int implHCount = bondAtomInSugars.getImplicitHydrogenCount();
                            bondAtomInSugars.setImplicitHydrogenCount(implHCount + bond.getOrder().numeric());
                        } else {
                            SugarDetectionUtility.LOGGER.error("Bond atom neither found in aglycone nor in sugars, this should not happen!");
                        }
                    }
                }
            } //end of if condition looking for bonds broken during sugar extraction
        } // end of for loop over all bonds in the input molecule
        if (!hasIdentifiedBrokenBond && !copyForAglycone.isEmpty() && ConnectivityChecker.isConnected(mol)) {
            //note for disconnected glycosides, one could process each component separately, but this seems like
            // unnecessary overhead just for the sake of this check
            SugarDetectionUtility.LOGGER.error("No broken bonds found between aglycone and sugars, no saturation performed, this should not happen!");
        }
        if (postProcessSugars) {
            if (extractLinearSugars) {
                copyForSugars = this.splitEtherEsterAndPeroxideBondsPostProcessing(copyForSugars, markAttachPointsByR);
            }
            if (extractCircularSugars) {
                copyForSugars = this.splitOGlycosidicBonds(copyForSugars, markAttachPointsByR);
            }
        }
        //clean up the maps
        for (IAtom atom : mol.atoms()) {
            if (!copyForAglycone.contains(inputAtomToAtomCopyInAglyconeMap.get(atom))) {
                inputAtomToAtomCopyInAglyconeMap.remove(atom);
            }
        }
        for (IBond bond : mol.bonds()) {
            if (!copyForAglycone.contains(inputBondToBondCopyInAglyconeMap.get(bond))) {
                inputBondToBondCopyInAglyconeMap.remove(bond);
            }
        }
        for (IAtom atom : mol.atoms()) {
            if (!copyForSugars.contains(inputAtomToAtomCopyInSugarsMap.get(atom))) {
                inputAtomToAtomCopyInSugarsMap.remove(atom);
            }
        }
        for (IBond bond : mol.bonds()) {
            if (!copyForSugars.contains(inputBondToBondCopyInSugarsMap.get(bond))) {
                inputBondToBondCopyInSugarsMap.remove(bond);
            }
        }
        //return value preparations, partition disconnected sugars
        List<IAtomContainer> resultsList = new ArrayList<>(5);
        resultsList.add(0, copyForAglycone);
        if (ConnectivityChecker.isConnected(copyForSugars)) {
            resultsList.add(copyForSugars);
        } else {
            for (IAtomContainer part : ConnectivityChecker.partitionIntoMolecules(copyForSugars)) {
                if (!part.isEmpty()) {
                    resultsList.add(part);
                }
            }
        }
        return resultsList;
    }

    /**
     * Returns the indices of atoms in the input molecule that correspond to atoms in the given group.
     * <p>
     * This method iterates through all atoms in the input molecule and checks if the corresponding
     * atom (via the provided mapping) exists in the group container. The indices of matching atoms
     * are collected and returned as an array.
     * <p>
     * Note that the group may contain atoms that are not present in the input molecule mapping
     * (e.g., R-groups added during processing), which will be ignored.
     *
     * @param mol The input molecule containing the original atoms
     * @param group The group container to check for atom membership
     * @param inputAtomToAtomCopyMap Map from original atoms to their copies in the group
     * @return Array of atom indices in the input molecule that have corresponding atoms in the group.
     *         Returns empty array if no matching atoms are found or if the group is empty.
     * @throws NullPointerException if any of the parameters is null
     */
    public int[] getGroupAtomIndices(IAtomContainer mol, IAtomContainer group, Map<IAtom, IAtom> inputAtomToAtomCopyMap) {
        if (mol == null || group == null || inputAtomToAtomCopyMap == null) {
            throw new NullPointerException("Given molecule, group, or input atom to atom copy map is null.");
        }
        if (group.isEmpty()) {
            return new int[0];
        }
        //cannot immediately use array because the group may contain atoms that are not in the input molecule, e.g. R atoms
        ArrayList<Integer> groupAtomIndices = new ArrayList<>(group.getAtomCount());
        for (IAtom atom : mol.atoms()) {
            if (group.contains(inputAtomToAtomCopyMap.get(atom))) {
                groupAtomIndices.add(atom.getIndex());
            }
        }
        int[] indices = new int[groupAtomIndices.size()];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = groupAtomIndices.get(i);
        }
        return indices;
    }

    /**
     * Returns the indices of bonds in the input molecule that correspond to bonds in the given group.
     * <p>
     * This method iterates through all bonds in the input molecule and checks if the corresponding
     * bond (via the provided mapping) exists in the group container. The indices of matching bonds
     * are collected and returned as an array.
     * <p>
     * Note that the group may contain bonds that are not present in the input molecule mapping
     * (e.g., bonds to R-groups added during processing), which will be ignored.
     *
     * @param mol The input molecule containing the original bonds
     * @param group The group container to check for bond membership
     * @param inputBondToBondCopyMap Map from original bonds to their copies in the group
     * @return Array of bond indices in the input molecule that have corresponding bonds in the group.
     *         Returns empty array if no matching bonds are found or if the group is empty.
     * @throws NullPointerException if any of the parameters is null
     */
    public int[] getGroupBondIndices(IAtomContainer mol, IAtomContainer group, Map<IBond, IBond> inputBondToBondCopyMap) {
        if (mol == null || group == null || inputBondToBondCopyMap == null) {
            throw new NullPointerException("Given molecule, group, or input bond to bond copy map is null.");
        }
        if (group.isEmpty()) {
            return new int[0];
        }
        //cannot immediately use array because the group may contain bonds that are not in the input molecule, e.g. bonds to R atoms
        ArrayList<Integer> groupBondIndices = new ArrayList<>(group.getBondCount());
        for (IBond bond : mol.bonds()) {
            if (group.contains(inputBondToBondCopyMap.get(bond))) {
                groupBondIndices.add(bond.getIndex());
            }
        }
        int[] indices = new int[groupBondIndices.size()];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = groupBondIndices.get(i);
        }
        return indices;
    }

    /**
     * Creates a relatively deep ("deeper" than cloning) copy of the given atom container mol and fills the given maps
     * with a mappings of the original atoms and bonds to the atoms an d bonds in the copy.
     * Copies:
     * <br>- Atoms (atomic number, implicit hydrogen count, aromaticity flag, valency, atom type name, formal charge, some primitive-based properties)
     * <br>- Bonds (begin and end atom, order, aromaticity flag, stereo, display, in ring flag, some primitive-based properties)
     * <br>- Single electrons
     * <br>- Lone pairs
     * <br>- Stereo elements (mapped to the copied atoms and bonds)
     * <br>- Some primitive-based properties (String, Integer, Boolean)
     * <br>Note: atom types of the original atoms are not copied and hence, some properties will be unset in the copies.
     * If you need atom types and their defining properties, you need to re-perceive them after copying.
     *
     * @param mol the molecule to copy
     * @param origToCopyAtomMap empty map to fill with a mapping of the original atoms to the copied atoms
     * @param origToCopyBondMap empty map to fill with a mapping of the original bonds to the copied bonds
     * @return a relatively deep copy of the given atom container
     */
    protected IAtomContainer deeperCopy(
            IAtomContainer mol,
            Map<IAtom, IAtom> origToCopyAtomMap,
            Map<IBond, IBond> origToCopyBondMap) {
        IAtomContainer copy = mol.getBuilder().newAtomContainer();
        // atoms
        for (IAtom atom : mol.atoms()) {
            IAtom cpyAtom = this.deeperCopy(atom, copy);
            origToCopyAtomMap.put(atom, cpyAtom);
        }
        // bonds
        for (IBond bond : mol.bonds()) {
            IAtom beg = origToCopyAtomMap.get(bond.getBegin());
            IAtom end = origToCopyAtomMap.get(bond.getEnd());
            if (beg == null || end == null || beg.getContainer() != end.getContainer()) {
                continue;
            }
            IBond newBond = this.deeperCopy(bond, beg, end);
            copy.addBond(newBond);
            origToCopyBondMap.put(bond, newBond);
        }
        // single electrons
        for (ISingleElectron se : mol.singleElectrons()) {
            IAtom atom = origToCopyAtomMap.get(se.getAtom());
            if (!Objects.isNull(atom)) {
                atom.getContainer().addSingleElectron(atom.getIndex());
            }
        }
        // lone pairs
        for (ILonePair lp : mol.lonePairs()) {
            IAtom atom = origToCopyAtomMap.get(lp.getAtom());
            if (!Objects.isNull(atom)) {
                atom.getContainer().addLonePair(atom.getIndex());
            }
        }
        // stereo elements
        for (IStereoElement elem : mol.stereoElements()) {
            copy.addStereoElement(elem.map(origToCopyAtomMap, origToCopyBondMap));
        }
        // properties
        for (Map.Entry<Object, Object> entry : mol.getProperties().entrySet()) {
            if ((entry.getKey() instanceof String || entry.getKey() instanceof Integer || entry.getKey() instanceof Boolean)
                    && (entry.getValue() instanceof String || entry.getValue() instanceof Integer || entry.getValue() instanceof Boolean || entry.getValue() == null)) {
                copy.setProperty(entry.getKey(), entry.getValue());
            }
        }
        return copy;
    }

    /**
     *  Creates a relatively deep ("deeper" than cloning) copy of the given atom and adds it to the given container.
     *  Copies:
     *  <br>- atomic number
     *  <br>- implicit hydrogen count
     *  <br>- aromaticity flag
     *  <br>- valency
     *  <br>- atom type name
     *  <br>- formal charge
     *  <br>- some primitive-based properties (String, Integer, Boolean)
     * <br>Note: atom types of the original atoms are not copied and hence, some properties will be unset in the copies.
     * If you need atom types and their defining properties, you need to re-perceive them after copying.
     *
     * @param atom the atom to copy
     * @param container the container to add the copied atom to
     * @return the copied atom
     */
    protected IAtom deeperCopy(IAtom atom, IAtomContainer container) {
        IAtom cpyAtom = container.newAtom(atom.getAtomicNumber(),
                atom.getImplicitHydrogenCount());
        cpyAtom.setIsAromatic(atom.isAromatic());
        cpyAtom.setValency(atom.getValency());
        cpyAtom.setAtomTypeName(atom.getAtomTypeName());
        cpyAtom.setFormalCharge(atom.getFormalCharge());
        // properties
        for (Map.Entry<Object, Object> entry : atom.getProperties().entrySet()) {
            if ((entry.getKey() instanceof String || entry.getKey() instanceof Integer || entry.getKey() instanceof Boolean)
                    && (entry.getValue() instanceof String || entry.getValue() instanceof Integer || entry.getValue() instanceof Boolean || entry.getValue() == null)) {
                cpyAtom.setProperty(entry.getKey(), entry.getValue());
            }
        }
        return cpyAtom;
    }

    /**
     * Creates a relatively deep ("deeper" than cloning) copy of the given bond between the given begin and end atoms.
     * Copies:
     * <br>- order
     * <br>- aromaticity flag
     * <br>- stereo
     * <br>- display
     * <br>- in ring flag
     * <br>- some primitive-based properties (String, Integer, Boolean)
     * <br>Note: The begin and end atoms are not copied, but the given ones are used in the copy.
     * <br>Note also: the created bond must be added to the copy atom container by the calling code!
     *
     * @param bond the bond to copy
     * @param begin the begin atom of the bond in the copy(!)
     * @param end the end atom of the bond in the copy(!)
     * @return the copied bond
     */
    protected IBond deeperCopy(IBond bond, IAtom begin, IAtom end) {
        //using begin.getContainer().newBond() here caused weird issues sometimes
        IBond newBond = new Bond(begin, end, bond.getOrder());
        newBond.setIsAromatic(bond.isAromatic());
        newBond.setStereo(bond.getStereo());
        newBond.setDisplay(bond.getDisplay());
        newBond.setIsInRing(bond.isInRing());
        // properties
        for (Map.Entry<Object, Object> entry : bond.getProperties().entrySet()) {
            if ((entry.getKey() instanceof String || entry.getKey() instanceof Integer || entry.getKey() instanceof Boolean)
                    && (entry.getValue() instanceof String || entry.getValue() instanceof Integer || entry.getValue() instanceof Boolean || entry.getValue() == null)) {
                newBond.setProperty(entry.getKey(), entry.getValue());
            }
        }
        return newBond;
    }

    /**
     * Checks whether the given atom is a hetero-atom (i.e. non-carbon and
     * non-hydrogen). Pseudo (R) atoms will also return false.
     *
     * @param atom the atom to test
     * @return true if the given atom is neither a carbon nor a hydrogen or
     *         pseudo atom
     */
    protected boolean isHeteroAtom(IAtom atom) {
        Integer tmpAtomicNr = atom.getAtomicNumber();
        if (Objects.isNull(tmpAtomicNr)) {
            return false;
        }
        int tmpAtomicNumberInt = tmpAtomicNr;
        return tmpAtomicNumberInt != IElement.H && tmpAtomicNumberInt != IElement.C
                && !this.isPseudoAtom(atom);
    }

    /**
     * Checks whether the given atom is a pseudo atom. Very strict, any atom
     * whose atomic number is null or 0, whose symbol equals "R" or "*", or that
     * is an instance of an IPseudoAtom implementing class will be classified as
     * a pseudo atom.
     *
     * @param atom the atom to test
     * @return true if the given atom is identified as a pseudo (R) atom
     */
    protected boolean isPseudoAtom(IAtom atom) {
        Integer tmpAtomicNr = atom.getAtomicNumber();
        if (Objects.isNull(tmpAtomicNr)) {
            return true;
        }
        String tmpSymbol = atom.getSymbol();
        return tmpAtomicNr == IElement.Wildcard ||
                tmpSymbol.equals("R") ||
                tmpSymbol.equals("*") ||
                atom instanceof IPseudoAtom;
    }

    /**
     * Checks whether the given atom is a carbon atom.
     *
     * @param atom the atom to test
     * @return true if the given atom is a carbon atom
     */
    protected boolean isCarbonAtom(IAtom atom) {
        Integer tmpAtomicNr = atom.getAtomicNumber();
        if (Objects.isNull(tmpAtomicNr)) {
            return false;
        }
        int tmpAtomicNumberInt = tmpAtomicNr;
        return tmpAtomicNumberInt == IElement.C;
    }

    /**
     * Splits O-glycosidic bonds in the given molecule (circular sugar moieties) and optionally marks the attachment points with R-groups.
     * This method identifies O-glycosidic bonds in the molecule using a SMARTS pattern and applies a SMIRKS transformation
     * to break these bonds. The transformation can either mark the attachment points with R-groups or saturate the
     * resulting open valences with implicit H atoms, depending on the `markAttachPointsByR` parameter.
     * If bonds are split, an unconnected atom container is returned. If no O-glycosidic bonds are found, the original
     * molecule is returned unchanged.
     *
     * @param molecule The molecule in which O-glycosidic bonds are to be split.
     * @param markAttachPointsByR If true, the attachment points are marked with R-groups; otherwise, they are saturated with implicit H.
     * @return An unconnected structure with the O-glycosidic bonds split, or the original molecule if no transformation was applied.
     */
    protected IAtomContainer splitOGlycosidicBonds(IAtomContainer molecule, boolean markAttachPointsByR) {
        //an aliphatic C in a ring with degree 3 and no charge, connected to an aliphatic O not in a ring with degree 2 and no charge,
        // connected to an aliphatic C with no charge (this side is left more promiscuous for corner cases)
        String eductPattern = "[C;R;D3;+0:1]-!@[O;!R;D2;+0:2]-!@[C;+0:3]";
        SmirksTransform transformation;
        if (markAttachPointsByR) {
            //transformed into two hydroxy groups with R atoms
            transformation = Smirks.compile(eductPattern + ">>([C:1]-O-*).(*-O-[C:3])");
        } else {
            //transformed into two saturated hydroxy groups
            transformation = Smirks.compile(eductPattern + ">>([C:1]-[OH1]).([OH1]-[C:3])");
        }
        //to detect rings
        transformation.setPrepare(true);
        for (IAtomContainer result : transformation.apply(molecule, Transform.Mode.Exclusive)) {
            //the iterable is empty if the educt pattern does not match, so this is not reached;
            // if it matches, only one result is in the iterable because of Mode.Exclusive, so this one result should be returned
            molecule = result;
        }
        //return the original molecule if no transformation was applied
        return molecule;
    }

    /**
     * Splits ether, ester, and peroxide bonds in the given molecule (linear sugar moieties) and optionally marks the
     * attachment points with R-groups.
     * This method identifies specific bond types (ether, ester, and peroxide) in the molecule using SMARTS patterns and applies
     * SMIRKS transformations to break these bonds. The transformation can either mark the attachment points with R-groups or
     * saturate the resulting open valences with implicit H atoms, depending on the `markAttachPointsByR` parameter.
     * If bonds are split, an unconnected atom container is returned. If no matching bonds are found, the original molecule
     * is returned unchanged.
     *
     * @param molecule The molecule in which ether, ester, and peroxide bonds are to be split.
     * @param markAttachPointsByR If true, the attachment points are marked with R-groups; otherwise, they are saturated with implicit H.
     * @return An unconnected structure with the specified bonds split, or the original molecule if no transformation was applied.
     */
    protected IAtomContainer splitEtherEsterAndPeroxideBondsPostProcessing(IAtomContainer molecule, boolean markAttachPointsByR) {
        //TODO put patterns into constants
        String esterEductPattern = "[C;!R;+0:1](=!@[O;!R;+0:2])-!@[O;!R;D2;+0:3]-!@[C;!R;+0:4]";
        String etherEductCrosslinkPattern = "[C;!R;+0:1]-!@[O;!R;D2;+0:2]-!@[C;!R;+0:3]-!@[OH1;!R;+0:4]";
        String etherEductPattern = "[C;!R;+0:1]-!@[O;!R;D2;+0:2]-!@[C;!R;+0:3]";
        String peroxideEductPattern = "[C;!R;+0:1]-!@[O;!R;D2;+0:2]-!@[O;!R;D2;+0:3]-!@[C;!R;+0:4]";
        SmirksTransform esterTransformation;
        SmirksTransform etherCrosslinkTransformation;
        SmirksTransform etherTransformation;
        SmirksTransform peroxideTransformation;
        if (markAttachPointsByR) {
            //transformed into two hydroxy groups with R atoms
            esterTransformation = Smirks.compile(esterEductPattern + ">>([C:1](=[O:2])-O-*).(*-O-[C:4])");
            etherCrosslinkTransformation = Smirks.compile(etherEductCrosslinkPattern + ">>([C:1]-[O:2]-*).(*-[C:3]-[OH1:4])");
            etherTransformation = Smirks.compile(etherEductPattern + ">>([C:1]-O-*).(*-O-[C:3])");
            peroxideTransformation = Smirks.compile(peroxideEductPattern + ">>([C:1]-[O:2]-*).(*-[O:3]-[C:4])");
        } else {
            //transformed into two saturated hydroxy groups
            esterTransformation = Smirks.compile(esterEductPattern + ">>([C:1](=[O:2])-[OH1]).([OH1]-[C:4])");
            etherCrosslinkTransformation = Smirks.compile(etherEductCrosslinkPattern + ">>([C:1]-[OH1:2]).([H][C:3]-[OH1:4])");
            etherTransformation = Smirks.compile(etherEductPattern + ">>([C:1]-[OH1]).([OH1]-[C:3])");
            peroxideTransformation = Smirks.compile(peroxideEductPattern + ">>([C:1]-[OH1:2]).([OH1:3]-[C:4])");
        }
        SmirksTransform[] transformations = new SmirksTransform[] {
                esterTransformation, etherCrosslinkTransformation, etherTransformation, peroxideTransformation
        };
        for (SmirksTransform transformation : transformations) {
            //to detect rings
            transformation.setPrepare(true);
            for (IAtomContainer result : transformation.apply(molecule, Transform.Mode.Exclusive)) {
                //the iterable is empty if the educt pattern does not match, so this is not reached;
                // if it matches, only one result is in the iterable because of Mode.Exclusive, so this one result should be returned
                molecule = result;
            }
        }
        return molecule;
    }

    /**
     *
     */
    protected void splitEsters(IAtomContainer molecule, boolean markAttachPointsByR) {
        //SmartsPattern.prepare should be done before calling this method
        String esterEductPattern = "[C;!R;+0;$(C=!@[O;!R;+0]):1]-!@[O;!R;D2;+0:2]-!@[C;!R;+0:3]";
        Mappings esterMappings = SmartsPattern.create(esterEductPattern).matchAll(molecule).uniqueAtoms();
        if (esterMappings.atLeast(1)) {
            for (IAtomContainer esterGroup : esterMappings.toSubstructures()) {
                IAtom carbonOne = null;
                IAtom connectingOxygen = null;
                for (IAtom atom : esterGroup.atoms()) {
                    if (atom.getAtomicNumber() == IElement.O) {
                        connectingOxygen = atom;
                    } else if (carbonOne == null ) {
                        carbonOne = atom;
                    }
                }
                molecule.removeBond(carbonOne, connectingOxygen);
                IAtom newOxygen = molecule.newAtom(IElement.O);
                molecule.newBond(carbonOne, newOxygen, IBond.Order.SINGLE);
                IAtom[] oxygens = new IAtom[] {connectingOxygen, newOxygen};
                for (IAtom oxygen : oxygens) {
                    if (markAttachPointsByR) {
                        IPseudoAtom tmpRAtom = oxygen.getBuilder().newInstance(IPseudoAtom.class, "R");
                        tmpRAtom.setAttachPointNum(1);
                        tmpRAtom.setImplicitHydrogenCount(0);
                        IBond bondToR = oxygen.getBuilder().newInstance(
                                IBond.class, oxygen, tmpRAtom, IBond.Order.SINGLE);
                        molecule.addAtom(tmpRAtom);
                        molecule.addBond(bondToR);
                        oxygen.setImplicitHydrogenCount(0);
                    } else {
                        oxygen.setImplicitHydrogenCount(1);
                    }
                }
            }
        }
    }

    /**
     *
     */
    protected void splitEtherCrosslinking(IAtomContainer molecule, boolean markAttachPointsByR) {
        String etherEductCrosslinkPattern = "[C;!R;+0:1]-!@[O;!R;D2;+0:2]-!@[C;!R;+0:3]-!@[OH1;!R;+0:4]";
        //>>([C:1]-[O:2]-*).(*-[C:3]-[OH1:4])
        //>>([C:1]-[OH1:2]).([H][C:3]-[OH1:4])
    }
}
