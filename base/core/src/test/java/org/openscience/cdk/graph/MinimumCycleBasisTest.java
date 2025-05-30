/*
 * Copyright (c) 2013 European Bioinformatics Institute (EMBL-EBI)
 *                    John May <jwmay@users.sf.net>
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
package org.openscience.cdk.graph;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.openscience.cdk.graph.InitialCyclesTest.anthracene;
import static org.openscience.cdk.graph.InitialCyclesTest.bicyclo;
import static org.openscience.cdk.graph.InitialCyclesTest.cyclophane_even;
import static org.openscience.cdk.graph.InitialCyclesTest.naphthalene;
import static org.openscience.cdk.graph.InitialCyclesTest.norbornane;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * @author John May
 */
class MinimumCycleBasisTest {

    @Test
    void noGraph() {
        Assertions.assertThrows(NullPointerException.class,
                                () -> {
                                    new MinimumCycleBasis((int[][]) null);
                                });
    }

    @Test
    void noInitialCycles() {
        Assertions.assertThrows(NullPointerException.class,
                                () -> {
                                    new MinimumCycleBasis((InitialCycles) null);
                                });
    }

    @Test
    void paths_norbornane() {
        int[][] norbornane = norbornane();
        MinimumCycleBasis mcb = new MinimumCycleBasis(norbornane);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{5, 6, 2, 1, 0, 5}, {5, 6, 2, 3, 4, 5}};
        assertThat(paths, is(expected));
    }

    @Test
    void paths_bicyclo() {
        int[][] bicyclo = bicyclo();
        MinimumCycleBasis mcb = new MinimumCycleBasis(bicyclo);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{5, 0, 1, 2, 3, 4, 5}, {5, 0, 1, 2, 7, 6, 5}};
        assertThat(paths, is(expected));
    }

    @Test
    void paths_napthalene() {
        int[][] napthalene = naphthalene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(napthalene);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{5, 0, 1, 2, 3, 4, 5}, {5, 4, 7, 8, 9, 6, 5}};
        assertThat(paths, is(expected));
    }

    @Test
    void paths_anthracene() {
        int[][] anthracene = anthracene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(anthracene);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(3));
        int[][] expected = new int[][]{{5, 0, 1, 2, 3, 4, 5}, {9, 6, 5, 4, 7, 8, 9}, {9, 8, 10, 11, 12, 13, 9}};
        assertThat(paths, is(expected));
    }

    @Test
    void paths_cyclophane_even() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis mcb = new MinimumCycleBasis(cyclophane_even);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{3, 2, 1, 0, 5, 4, 3}, {3, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3}};
        assertThat(paths, is(expected));
    }

    @Test
    void paths_cyclophane_odd() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis mcb = new MinimumCycleBasis(cyclophane_even);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
        int[][] expected = new int[][]{{3, 2, 1, 0, 5, 4, 3}, {3, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3}};
        assertThat(paths, is(expected));
    }

    @Test
    void size_norbornane() {
        int[][] norbornane = norbornane();
        MinimumCycleBasis mcb = new MinimumCycleBasis(norbornane);
        int[][] paths = mcb.paths();
        assertThat(paths.length, is(2));
    }

    @Test
    void size_bicyclo() {
        int[][] bicyclo = bicyclo();
        MinimumCycleBasis mcb = new MinimumCycleBasis(bicyclo);
        assertThat(mcb.size(), is(2));
    }

    @Test
    void size_napthalene() {
        int[][] napthalene = naphthalene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(napthalene);
        assertThat(mcb.size(), is(2));
    }

    @Test
    void size_anthracene() {
        int[][] anthracene = anthracene();
        MinimumCycleBasis mcb = new MinimumCycleBasis(anthracene);
        assertThat(mcb.size(), is(3));
    }

    @Test
    void size_cyclophane_even() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis relevant = new MinimumCycleBasis(cyclophane_even);
        assertThat(relevant.size(), is(2));
    }

    @Test
    void size_cyclophane_odd() {
        int[][] cyclophane_even = cyclophane_even();
        MinimumCycleBasis mcb = new MinimumCycleBasis(cyclophane_even);
        assertThat(mcb.size(), is(2));
    }
}
