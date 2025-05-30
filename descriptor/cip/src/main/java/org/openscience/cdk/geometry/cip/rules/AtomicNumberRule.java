/* Copyright (C) 2010  Egon Willighagen <egonw@users.sf.net>
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
 */
package org.openscience.cdk.geometry.cip.rules;

import org.openscience.cdk.geometry.cip.ILigand;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 * Compares to {@link ILigand}s based on atomic numbers.
 *
 */
class AtomicNumberRule implements ISequenceSubRule<ILigand> {

    /** {@inheritDoc} */
    @Override
    public int compare(ILigand ligand1, ILigand ligand2) {
        return getAtomicNumber(ligand1).compareTo(getAtomicNumber(ligand2));
    }

    private Integer getAtomicNumber(ILigand ligand) {
        Integer atomNumber = ligand.getLigandAtom().getAtomicNumber();
        if (atomNumber != null) return atomNumber;
        return PeriodicTable.getAtomicNumber(ligand.getLigandAtom().getSymbol());
    }
}
