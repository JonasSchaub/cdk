/*
 * Copyright (C) 2012 John May <jwmay@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version. All we ask is that proper credit is given for our
 * work, which includes - but is not limited to - adding the above copyright
 * notice to the beginning of your source code files, and to any copyright
 * notice that you may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.ringsearch;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.templates.TestMoleculeFactory;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 * ring search unit tests for a branched aliphatic compounds
 *
 * @author John May
 */
final class RingSearch_NonCyclicTest {

    private final IAtomContainer nonCyclic = TestMoleculeFactory.makeBranchedAliphatic();

    @Test
    void testCyclic() {
        assertThat(new RingSearch(nonCyclic).cyclic().length, is(0));
    }

    @Test
    void testCyclic_Int() {
        int n = nonCyclic.getAtomCount();
        RingSearch ringSearch = new RingSearch(nonCyclic);
        for (int i = 0; i < n; i++) {
            Assertions.assertFalse(ringSearch.cyclic(i));
        }
    }

    @Test
    void testIsolated() {
        assertThat(new RingSearch(nonCyclic).isolated().length, is(0));
    }

    @Test
    void testFused() {
        assertThat(new RingSearch(nonCyclic).fused().length, is(0));
    }

    @Test
    void testRingFragments() {
        Assertions.assertTrue(new RingSearch(nonCyclic).ringFragments().isEmpty());
    }

    @Test
    void testIsolatedRingFragments() {
        Assertions.assertTrue(new RingSearch(nonCyclic).isolatedRingFragments().isEmpty());
    }

    @Test
    void testFusedRingFragments() {
        Assertions.assertTrue(new RingSearch(nonCyclic).fusedRingFragments().isEmpty());
    }

}
