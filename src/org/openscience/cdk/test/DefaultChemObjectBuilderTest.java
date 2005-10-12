/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 * 
 * Copyright (C) 2002-2005  The Chemistry Development Kit (CDK) project
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. 
 * 
 */

package org.openscience.cdk.test;

import junit.framework.Test;
import junit.framework.TestSuite;

import org.openscience.cdk.interfaces.Atom;
import org.openscience.cdk.interfaces.AtomParity;
import org.openscience.cdk.interfaces.AtomType;
import org.openscience.cdk.interfaces.BioPolymer;
import org.openscience.cdk.interfaces.ChemFile;
import org.openscience.cdk.interfaces.ChemModel;
import org.openscience.cdk.interfaces.ChemObject;
import org.openscience.cdk.interfaces.ChemSequence;
import org.openscience.cdk.interfaces.ElectronContainer;
import org.openscience.cdk.interfaces.Element;
import org.openscience.cdk.interfaces.Isotope;
import org.openscience.cdk.interfaces.LonePair;
import org.openscience.cdk.interfaces.Molecule;
import org.openscience.cdk.interfaces.Monomer;
import org.openscience.cdk.interfaces.Polymer;
import org.openscience.cdk.interfaces.PseudoAtom;
import org.openscience.cdk.interfaces.Reaction;
import org.openscience.cdk.interfaces.Ring;
import org.openscience.cdk.interfaces.RingSet;
import org.openscience.cdk.interfaces.SetOfAtomContainers;
import org.openscience.cdk.interfaces.SetOfMolecules;
import org.openscience.cdk.interfaces.SetOfReactions;
import org.openscience.cdk.interfaces.SingleElectron;
import org.openscience.cdk.interfaces.Strand;

/**
 * Checks the funcitonality of the Crystal.
 *
 * @cdk.module test
 */
public class DefaultChemObjectBuilderTest extends CDKTestCase {

	private static org.openscience.cdk.ChemObject rootObject;
	
    public DefaultChemObjectBuilderTest(String name) {
        super(name);
        rootObject = new org.openscience.cdk.ChemObject();
    }

    public void setUp() {}

    public static Test suite() {
        return new TestSuite(DefaultChemObjectBuilderTest.class);
    }

	public void testNewAtom() {
		Object object = rootObject.getBuilder().newAtom();
		assertNotNull(object);
		assertTrue(object instanceof org.openscience.cdk.ChemObject);

		assertTrue(object instanceof Atom);
	}
	    
}
