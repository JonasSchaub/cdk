/* Copyright (C) 2008  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@slists.sourceforge.net
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
 */
package org.openscience.cdk.io.formats;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 */
public class PubChemSubstancesXMLFormatTest extends ChemFormatMatcherTest {

    PubChemSubstancesXMLFormatTest() {
        super.setChemFormatMatcher((IChemFormatMatcher) PubChemSubstancesXMLFormat.getInstance());
    }

    /**
     * @cdk.bug 2832835
     */
    @Test
    @Override
    void testMatches() throws Exception {
        String header = "<?xml version=\"\"?><PC-Substances/>";
        Assertions.assertTrue(matches(header));
    }

}
