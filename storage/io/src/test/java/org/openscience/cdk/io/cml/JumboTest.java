/* Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *  */
package org.openscience.cdk.io.cml;

import java.io.InputStream;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.test.CDKTestCase;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * TestCase for reading CML files using a few test files
 * in data/cmltest as found in the original Jumbo3 release
 * (http://www.xml-cml.org/).
 *
 */
class JumboTest extends CDKTestCase {

    private final ILoggingTool logger = LoggingToolFactory.createLoggingTool(JumboTest.class);

    /**
     * Now come the actual tests...
     */

    /**
     * Special CML characteristics:
     * - <atomArray><atom/><atom/></atomArray>
     * - X2D only
     */
    @Test
    void testCuran() throws Exception {
        String filename = "curan.xml";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getResourceAsStream(filename);
        CMLReader reader = new CMLReader(ins);
        IChemFile chemFile = reader.read(new org.openscience.cdk.ChemFile());
        reader.close();

        // test the resulting ChemFile content
        Assertions.assertNotNull(chemFile);
        Assertions.assertEquals(chemFile.getChemSequenceCount(), 1);
        IChemSequence seq = chemFile.getChemSequence(0);
        Assertions.assertNotNull(seq);
        Assertions.assertEquals(seq.getChemModelCount(), 1);
        IChemModel model = seq.getChemModel(0);
        Assertions.assertNotNull(model);
        Assertions.assertEquals(model.getMoleculeSet().getAtomContainerCount(), 1);

        // test the molecule
        IAtomContainer mol = model.getMoleculeSet().getAtomContainer(0);
        Assertions.assertNotNull(mol);
        Assertions.assertEquals(mol.getAtomCount(), 24);
        Assertions.assertEquals(mol.getBondCount(), 28);
        Assertions.assertFalse(GeometryUtil.has3DCoordinates(mol));
        Assertions.assertTrue(GeometryUtil.has2DCoordinates(mol));
    }

    /**
     * Special CML characteristics:
     * - use of cml: namespace
     * - X2D only
     */
    @Test
    void testCephNS() throws Exception {
        String filename = "ceph-ns.xml";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getResourceAsStream(filename);
        CMLReader reader = new CMLReader(ins);
        IChemFile chemFile = reader.read(new org.openscience.cdk.ChemFile());
        reader.close();

        // test the resulting ChemFile content
        Assertions.assertNotNull(chemFile);
        Assertions.assertEquals(chemFile.getChemSequenceCount(), 1);
        IChemSequence seq = chemFile.getChemSequence(0);
        Assertions.assertNotNull(seq);
        Assertions.assertEquals(seq.getChemModelCount(), 1);
        IChemModel model = seq.getChemModel(0);
        Assertions.assertNotNull(model);
        Assertions.assertEquals(model.getMoleculeSet().getAtomContainerCount(), 1);

        // test the molecule
        IAtomContainer mol = model.getMoleculeSet().getAtomContainer(0);
        Assertions.assertNotNull(mol);
        Assertions.assertEquals(mol.getAtomCount(), 15);
        Assertions.assertEquals(mol.getBondCount(), 16);
        Assertions.assertFalse(GeometryUtil.has3DCoordinates(mol));
        Assertions.assertTrue(GeometryUtil.has2DCoordinates(mol));
    }

    /**
     * Special CML characteristics:
     * - <atomArray><stringArray builtin="atomId"/></atomArray>
     * - <bondArray><stringArray builtin="atomRef"/></atomArray>
     * - no coords
     */
    @Test
    void testNucleustest() throws Exception {
        String filename = "nucleustest.xml";
        logger.info("Testing: " + filename);
        InputStream ins = this.getClass().getResourceAsStream(filename);
        CMLReader reader = new CMLReader(ins);
        IChemFile chemFile = reader.read(new org.openscience.cdk.ChemFile());
        reader.close();

        // test the resulting ChemFile content
        Assertions.assertNotNull(chemFile);
        Assertions.assertEquals(chemFile.getChemSequenceCount(), 1);
        IChemSequence seq = chemFile.getChemSequence(0);
        Assertions.assertNotNull(seq);
        Assertions.assertEquals(seq.getChemModelCount(), 1);
        IChemModel model = seq.getChemModel(0);
        Assertions.assertNotNull(model);
        Assertions.assertEquals(model.getMoleculeSet().getAtomContainerCount(), 1);

        // test the molecule
        IAtomContainer mol = model.getMoleculeSet().getAtomContainer(0);
        Assertions.assertNotNull(mol);
        Assertions.assertEquals(11, mol.getAtomCount(), "Incorrect number of atoms");
        Assertions.assertEquals(12, mol.getBondCount(), "Incorrect number of bonds");
        Assertions.assertFalse(GeometryUtil.has3DCoordinates(mol), "File does not have 3D coordinates");
        Assertions.assertFalse(GeometryUtil.has2DCoordinates(mol), "File does not have 2D coordinates");
    }

}
