/* Copyright (C) 2006-2007  Miguel Rojas <miguel.rojas@uni-koeln.de>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.tools;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.reaction.type.PiBondingMovementReaction;
import org.openscience.cdk.reaction.type.RearrangementAnionReaction;
import org.openscience.cdk.reaction.type.RearrangementCationReaction;
import org.openscience.cdk.reaction.type.RearrangementLonePairReaction;
import org.openscience.cdk.reaction.type.RearrangementRadicalReaction;
import org.openscience.cdk.reaction.type.SharingLonePairReaction;
import org.openscience.cdk.reaction.type.parameters.IParameterReact;
import org.openscience.cdk.reaction.type.parameters.SetReactionCenter;

/**
 * <p>This class try to generate resonance structure for a determinate molecule.</p>
 * <p>Make sure that the molecule has the corresponding lone pair electrons
 * for each atom. You can use the method: <pre> LonePairElectronChecker </pre>
 * <p>It is needed to call the addExplicitHydrogensToSatisfyValency
 *  from the class tools.HydrogenAdder.</p>
 * <p>It is based on rearrangements of electrons and charge</p>
 * <p>The method is based on call by reactions which occur in a resonance.</p>
 *
 * <pre>
 * StructureResonanceGenerator srG = new StructureReseonanceGenerator(true,true,true,true,false);
 * MoleculeSet setOf = srG.getResonances(new Molecule());
 * </pre>
 *
 * <p>We have the possibility to localize the reactive center. Good method if you
 * want to localize the reaction in a fixed point</p>
 * <pre>atoms[0].setFlag(CDKConstants.REACTIVE_CENTER,true);</pre>
 * <p>Moreover you must put the parameter as true</p>
 * <p>If the reactive center is not localized then the reaction process will
 * try to find automatically the possible reactive center.</p>
 *
 * @author       Miguel Rojas
 * @cdk.created  2006-5-05
 *
 * @see org.openscience.cdk.reaction.IReactionProcess
 */
public class StructureResonanceGenerator {

    private final ILoggingTool           logger        = LoggingToolFactory
                                                         .createLoggingTool(StructureResonanceGenerator.class);
    private List<IReactionProcess> reactionsList = new ArrayList<>();
    /**Generate resonance structure without looking at the symmetry*/
    private final boolean                lookingSymmetry;
    /** TODO: REACT: some time takes too much time. At the moment fixed to 50 structures*/
    private int                    maxStructures = 50;

    /**
     * Construct an instance of StructureResonanceGenerator. Default restrictions
     * are initiated.
     *
     * @see #setDefaultReactions()
     */
    public StructureResonanceGenerator() {
        this(false);
    }

    /**
     * Construct an instance of StructureResonanceGenerator. Default restrictions
     * are initiated.
     *
     * @param lookingSymmetry  Specify if the resonance generation is based looking at the symmetry
     * @see #setDefaultReactions()
     */
    public StructureResonanceGenerator(boolean lookingSymmetry) {
        logger.info("Initiate StructureResonanceGenerator");
        this.lookingSymmetry = lookingSymmetry;
        setDefaultReactions();

    }

    /**
     * Set the reactions that must be used in the generation of the resonance.
     *
     * @param newReactionsList  The IReactionsProcess's to use
     *
     * @see #getReactions()
     * @see #setReactions(java.util.List)
     * @see IReactionProcess
     */
    public void setReactions(List<IReactionProcess> newReactionsList) {
        reactionsList = newReactionsList;
    }

    /**
     * Get the reactions that must be presents in the generation of the resonance.
     *
     * @return The reactions to be imposed
     *
     *
     * @see #setDefaultReactions()
     */
    public List<IReactionProcess> getReactions() {
        return this.reactionsList;
    }

    /**
     * Set the number maximal of resonance structures to be found. The
     * algorithm breaks the process when is came to this number.
     *
     * @param maxStruct The maximal number
     */
    public void setMaximalStructures(int maxStruct) {
        maxStructures = maxStruct;
    }

    /**
     * Get the number maximal of resonance structures to be found.
     *
     * @return The maximal number
     */
    public int getMaximalStructures() {
        return maxStructures;
    }

    /**
     * Set the default reactions that must be presents to generate the resonance.
     *
     * @see #getReactions()
     */
    public void setDefaultReactions() {
        callDefaultReactions();

    }

    private void callDefaultReactions() {
        List<IParameterReact> paramList = new ArrayList<>();
        IParameterReact param = new SetReactionCenter();
        param.setParameter(Boolean.FALSE);
        paramList.add(param);

        IReactionProcess type = new SharingLonePairReaction();
        try {
            type.setParameterList(paramList);
        } catch (CDKException e) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e);
        }
        reactionsList.add(type);

        type = new PiBondingMovementReaction();
        List<IParameterReact> paramList2 = new ArrayList<>();
        IParameterReact param2 = new SetReactionCenter();
        param2.setParameter(Boolean.FALSE);
        paramList2.add(param2);
        try {
            type.setParameterList(paramList2);
        } catch (CDKException e) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e);
        }
        reactionsList.add(type);

        type = new RearrangementAnionReaction();
        try {
            type.setParameterList(paramList);
        } catch (CDKException e) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e);
        }
        reactionsList.add(type);

        type = new RearrangementCationReaction();
        try {
            type.setParameterList(paramList);
        } catch (CDKException e) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e);
        }
        reactionsList.add(type);

        type = new RearrangementLonePairReaction();
        try {
            type.setParameterList(paramList);
        } catch (CDKException e) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e);
        }
        reactionsList.add(type);

        type = new RearrangementRadicalReaction();
        try {
            type.setParameterList(paramList);
        } catch (CDKException e) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e);
        }
        reactionsList.add(type);

    }

    /**
     * Get the resonance structures from an {@link IAtomContainer}.
     *
     * @param molecule The IAtomContainer to analyze
     * @return         The different resonance structures
     */
    public IAtomContainerSet getStructures(IAtomContainer molecule) {
        int countStructure = 0;
        IAtomContainerSet setOfMol = molecule.getBuilder().newInstance(IAtomContainerSet.class);
        setOfMol.addAtomContainer(molecule);

        for (int i = 0; i < setOfMol.getAtomContainerCount(); i++) {
            IAtomContainer mol = setOfMol.getAtomContainer(i);
            for (IReactionProcess aReactionsList : reactionsList) {
                IReactionProcess reaction = aReactionsList;
                IAtomContainerSet setOfReactants = molecule.getBuilder().newInstance(IAtomContainerSet.class);
                setOfReactants.addAtomContainer(mol);
                try {
                    IReactionSet setOfReactions = reaction.initiate(setOfReactants, null);
                    if (setOfReactions.getReactionCount() != 0)
                        for (int k = 0; k < setOfReactions.getReactionCount(); k++)
                            for (int j = 0; j < setOfReactions.getReaction(k).getProducts().getAtomContainerCount(); j++) {
                                IAtomContainer product = setOfReactions.getReaction(k).getProducts()
                                        .getAtomContainer(j);
                                if (!existAC(setOfMol, product)) {
                                    setOfMol.addAtomContainer(product);
                                    countStructure++;
                                    if (countStructure > maxStructures) return setOfMol;
                                }
                            }
                } catch (CDKException e) {
                    LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                                      .warn("Unexpected Error:", e);
                }
            }
        }
        return setOfMol;
    }

    /**
     * Get the container which is found resonance from a {@link IAtomContainer}.
     * It is based on looking if the order of the bond changes.
     *
     * @param molecule The IAtomContainer to analyze
     * @return         The different containers
     */
    public IAtomContainerSet getContainers(IAtomContainer molecule) {
        IAtomContainerSet setOfCont = molecule.getBuilder().newInstance(IAtomContainerSet.class);
        IAtomContainerSet setOfMol = getStructures(molecule);

        if (setOfMol.getAtomContainerCount() == 0) return setOfCont;

        /* extraction of all bonds which has been produced a changes of order */
        List<IBond> bondList = new ArrayList<>();
        for (int i = 1; i < setOfMol.getAtomContainerCount(); i++) {
            IAtomContainer mol = setOfMol.getAtomContainer(i);
            for (int j = 0; j < mol.getBondCount(); j++) {
                IBond bond = molecule.getBond(j);
                if (!mol.getBond(j).getOrder().equals(bond.getOrder())) {
                    if (!bondList.contains(bond)) bondList.add(bond);
                }
            }
        }

        if (bondList.size() == 0) return null;

        int[] flagBelonging = new int[bondList.size()];
        for (int i = 0; i < flagBelonging.length; i++)
            flagBelonging[i] = 0;
        int[] position = new int[bondList.size()];
        int maxGroup = 1;

        /* Analysis if the bond are linked together */
        List<IBond> newBondList = new ArrayList<>();
        newBondList.add(bondList.get(0));

        int pos = 0;
        for (int i = 0; i < newBondList.size(); i++) {

            if (i == 0)
                flagBelonging[i] = maxGroup;
            else {
                if (flagBelonging[position[i]] == 0) {
                    maxGroup++;
                    flagBelonging[position[i]] = maxGroup;
                }
            }

            IBond bondA = newBondList.get(i);
            for (int ato = 0; ato < 2; ato++) {
                IAtom atomA1 = bondA.getAtom(ato);
                List<IBond> bondA1s = molecule.getConnectedBondsList(atomA1);
                for (IBond bondB : bondA1s) {
                    if (!newBondList.contains(bondB)) for (int k = 0; k < bondList.size(); k++)
                        if (bondList.get(k).equals(bondB)) if (flagBelonging[k] == 0) {
                            flagBelonging[k] = maxGroup;
                            pos++;
                            newBondList.add(bondB);
                            position[pos] = k;

                        }
                }
            }
            //if it is final size and not all are added
            if (newBondList.size() - 1 == i) for (int k = 0; k < bondList.size(); k++)
                if (!newBondList.contains(bondList.get(k))) {
                    newBondList.add(bondList.get(k));
                    position[i + 1] = k;
                    break;
                }
        }
        /* creating containers according groups */
        for (int i = 0; i < maxGroup; i++) {
            IAtomContainer container = molecule.getBuilder().newInstance(IAtomContainer.class);
            for (int j = 0; j < bondList.size(); j++) {
                if (flagBelonging[j] != i + 1) continue;
                IBond bond = bondList.get(j);
                IAtom atomA1 = bond.getBegin();
                IAtom atomA2 = bond.getEnd();
                if (!container.contains(atomA1)) container.addAtom(atomA1);
                if (!container.contains(atomA2)) container.addAtom(atomA2);
                container.addBond(bond);
            }
            setOfCont.addAtomContainer(container);
        }
        return setOfCont;
    }

    /**
     * Get the container which the atom is found on resonance from a {@link IAtomContainer}.
     * It is based on looking if the order of the bond changes. Return null
     * is any is found.
     *
     * @param molecule The IAtomContainer to analyze
     * @param atom     The IAtom
     * @return         The container with the atom
     */
    public IAtomContainer getContainer(IAtomContainer molecule, IAtom atom) {
        IAtomContainerSet setOfCont = getContainers(molecule);
        if (setOfCont == null) return null;

        for (IAtomContainer container : setOfCont.atomContainers()) {
            if (container.contains(atom)) return container;
        }

        return null;
    }

    /**
     * Get the container which the bond is found on resonance from a {@link IAtomContainer}.
     * It is based on looking if the order of the bond changes. Return null
     * is any is found.
     *
     * @param molecule The IAtomContainer to analyze
     * @param bond     The IBond
     * @return         The container with the bond
     */
    public IAtomContainer getContainer(IAtomContainer molecule, IBond bond) {
        IAtomContainerSet setOfCont = getContainers(molecule);
        if (setOfCont == null) return null;

        for (IAtomContainer container : setOfCont.atomContainers()) {
            if (container.contains(bond)) return container;
        }

        return null;
    }

    /**
     * Search if the setOfAtomContainer contains the atomContainer
     *
     *
     * @param set            ISetOfAtomContainer object where to search
     * @param atomContainer  IAtomContainer to search
     * @return   			 True, if the atomContainer is contained
     */
    private boolean existAC(IAtomContainerSet set, IAtomContainer atomContainer) {

        IAtomContainer acClone = null;
        try {
            acClone = atomContainer.clone();
            if (!lookingSymmetry) { /* remove all aromatic flags */
                for (IAtom atom : acClone.atoms())
                    atom.setFlag(IChemObject.AROMATIC, false);
                for (IBond bond : acClone.bonds())
                    bond.setFlag(IChemObject.AROMATIC, false);
            }
        } catch (CloneNotSupportedException e1) {
            LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                              .warn("Unexpected Error:", e1);
        }

        for (int i = 0; i < acClone.getAtomCount(); i++)
            //			if(acClone.getAtom(i).getID() == null)
            acClone.getAtom(i).setID("" + acClone.indexOf(acClone.getAtom(i)));

        if (lookingSymmetry) {
            try {
                Aromaticity.cdkLegacy().apply(acClone);
            } catch (CDKException e) {
                LoggingToolFactory.createLoggingTool(StructureResonanceGenerator.class)
                                  .warn("Unexpected Error:", e);
            }
        } else {
            if (!lookingSymmetry) { /* remove all aromatic flags */
                for (IAtom atom : acClone.atoms())
                    atom.setFlag(IChemObject.AROMATIC, false);
                for (IBond bond : acClone.bonds())
                    bond.setFlag(IChemObject.AROMATIC, false);
            }
        }
        for (int i = 0; i < set.getAtomContainerCount(); i++) {
            IAtomContainer ss = set.getAtomContainer(i);
            for (int j = 0; j < ss.getAtomCount(); j++)
                //				if(ss.getAtom(j).getID() == null)
                ss.getAtom(j).setID("" + ss.indexOf(ss.getAtom(j)));

            try {

                if (!lookingSymmetry) {
                    QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolChargeIDQueryContainer(acClone);
                    if (new UniversalIsomorphismTester().isIsomorph(ss, qAC)) {
                        QueryAtomContainer qAC2 = QueryAtomContainerCreator
                                .createSymbolAndBondOrderQueryContainer(acClone);
                        if (new UniversalIsomorphismTester().isIsomorph(ss, qAC2)) return true;
                    }
                } else {
                    QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(acClone);
                    Aromaticity.cdkLegacy().apply(ss);
                    if (new UniversalIsomorphismTester().isIsomorph(ss, qAC)) return true;
                }

            } catch (CDKException e1) {
                System.err.println(e1);
                logger.error(e1.getMessage());
                logger.debug(e1);
            }
        }
        return false;
    }
}
