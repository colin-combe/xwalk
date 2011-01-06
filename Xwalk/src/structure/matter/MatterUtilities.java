/*
 * (C) 2010 Abdullah Kahraman
 *
 * This software is part of the open-source project "Xwalk". You can use this
 * software under the terms of the
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 * (http://creativecommons.org/licenses/by-nc-sa/3.0/).
 * This means that you
 * 1.) can copy, modify, distribute the software
 * 2.) must give credit to the author
 * 3.) must not use this work for commercial purposes
 * 4.) must license derivative works under the same or a similar license.
 *
 */

package structure.matter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import structure.constants.Constants;
import structure.constants.Constants.BondTypes;
import structure.constants.Constants.ElementTypes;
import structure.math.Mathematics;
import structure.math.Point3d;
import structure.matter.parameter.Element;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;


/**
 * A generic class that holds various methods to operate on classes defined in
 * the matter package.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public abstract class MatterUtilities {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected MatterUtilities() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------

    /**
     * Assess whether two atoms are of the same type.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object to assess the equivalence with.
     * @return {@code TRUE} only if both atoms have the same atom name and
     *         residue name, {@code FALSE} otherwise.
     */
    public static boolean equalsType(final Atom atom1, final Atom atom2) {
        if (atom1.getName().trim().equals(atom2.getName().trim())
            &&
            atom1.getResidueName().trim().equals(atom2.getResidueName().trim())
           ) {
               return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Assess whether two atoms are in the same residue.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object to assess the equivalence with.
     * @return {@code TRUE} only if both atoms have the same residue name,
     *         residue number and chain Id, {@code FALSE} otherwise.
     */
    public static boolean equalsResidue(final Atom atom1, final Atom atom2) {
        if (atom1.getResidueName().trim().equals(atom2.getResidueName().trim())
            &&
           (atom1.getResidueNumber() == atom2.getResidueNumber())
            &&
           (atom1.getChainId() == atom2.getChainId())) {
            return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Assess whether two atoms occupy the same point in Cartesian space.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object to assess the equivalence with.
     * @return {@code TRUE} only if both atoms have the same Cartesian
     *         coordinates, {@code FALSE} otherwise.
     */
    public static boolean equalsPosition(final Atom atom1, final Atom atom2) {
        return atom1.getPoint3d().equals(atom2.getPoint3d());
    }
    //--------------------------------------------------------------------------

    /**
     * Returns those two atoms that are closest in two AtomList objects.
     * @param list1
     *        - AtomList object holding a first list of atom coordinates.
     * @param list2
     *        - AtomList object holding a second list of atom coordinates.
     * @return AtomList object holding the coordinates of the two closest atoms
     *         in both atom lists.
     */
     public static AtomList getClosestAtomPair(final AtomList list1,
                                               final AtomList list2) {
        double minDist = Integer.MAX_VALUE;
        AtomList minList = new AtomList();
        minList.add(list1.get(0));
        minList.add(list2.get(0));

        for (Atom atom1 : list1) {
            for (Atom atom2 : list2) {
                double dist = Mathematics.distance(atom1.getPoint3d(),
                                                   atom2.getPoint3d()
                                                  );
                if (dist < minDist) {
                    minDist = dist;
                    minList.set(0, atom1);
                    minList.set(1, atom2);
                }
            }
        }
    return minList;
    }

    //--------------------------------------------------------------------------

    /**
     * Calculates the maximum Cartesian coordinates of an atom list.
     * @param coords
     *        - AtomList object holding the Cartesian coordinates of a list of
     *          atoms.
     * @return Point3d object holding the maximum coordinates in each Cartesian
     *         dimension in its XYZ fields.
     */
    public static Point3d getMaximumCooridnate(final AtomList coords) {
        double maxX = Integer.MIN_VALUE;
        double maxY = Integer.MIN_VALUE;
        double maxZ = Integer.MIN_VALUE;
        for (Atom atom : coords) {
             Point3d xyz = atom.getPoint3d();
             double r = atom.getVanDerWaalsRadius();
             maxX = Math.max(maxX, xyz.getX() + r);
             maxY = Math.max(maxY, xyz.getY() + r);
             maxZ = Math.max(maxZ, xyz.getZ() + r);
        }
        return new Point3d(maxX, maxY, maxZ);
    }
    //--------------------------------------------------------------------------

    /**
     * Calculates the minimum Cartesian coordinates of an atom list.
     * @param coords
     *        AtomList object holding the Cartesian coordinates of a list of
     *        atoms.
     * @return Point3d object holding the minimum coordinates in each Cartesian
     *         dimension in its XYZ fields.
     */
    public static Point3d getMinimumCooridnate(final AtomList coords) {
        double minX = Integer.MAX_VALUE;
        double minY = Integer.MAX_VALUE;
        double minZ = Integer.MAX_VALUE;
        for (Atom atom : coords) {
             Point3d xyz = atom.getPoint3d();
             double r = atom.getVanDerWaalsRadius();
             minX = Math.min(minX, xyz.getX() - r);
             minY = Math.min(minY, xyz.getY() - r);
             minZ = Math.min(minZ, xyz.getZ() - r);
        }
        return new Point3d(minX, minY, minZ);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the dimension of a protein complex.
     * @param complex
     *        - PolyPeptideList object to which dimension should be calculated.
     * @return double value representing the dimension of the protein complex.
     */
    public static double getDimension(final PolyPeptideList complex) {
        AtomList allAtoms = complex.getAllAtoms();
        Point3d max = MatterUtilities.getMaximumCooridnate(allAtoms);
        Point3d min = MatterUtilities.getMinimumCooridnate(allAtoms);
        Point3d diff = new Point3d(max.getX() - min.getX(),
                                   max.getY() - min.getY(),
                                   max.getZ() - min.getZ());
        double sum = Math.pow(diff.getX(), 2)
                   + Math.pow(diff.getY(), 2)
                   + Math.pow(diff.getZ(), 2);

        double dim = Math.sqrt(sum);
    return dim;
    }
    //--------------------------------------------------------------------------

    /**
     * Connects a list of atoms depending on the distance and their
     * van der Waals radii + PDB uncertainty factor.
     * @param atoms
     *        - AtomList object holding the coordinates of a PDB molecule.
     * @return Bond array holding all potential covalently bound atoms.
     */
    public static ArrayList < Bond > calculateBonds(final AtomList atoms) {
        // Covalent bonds in organic molecules should have a distance less than
        // 1.54 Angstroem between both atom centers.
        // (see http://en.wikipedia.org/wiki/Bond_length).
        // To account for resolution error distance is taken to be less than
        // 1.8 (0.28 Angstroem estimated standard error for X-ray structures)
        ArrayList < Bond > bonds = new ArrayList < Bond >();
        for (Atom atom1 : atoms) {
            for (Atom atom2 : atoms) {
                if (atom1 != atom2) {
                    if (atom1.getElement().getType() == ElementTypes.METAL
                        ||
                        atom2.getElement().getType() == ElementTypes.METAL) {
                        continue;
                    }
                    double dist = Mathematics.distance(atom1.getPoint3d(),
                                                       atom2.getPoint3d());
                    double maxDist;
                    if (atom1.getElement() == Element.HYDROGEN
                        ||
                        atom2.getElement() == Element.HYDROGEN) {
                        maxDist = Constants.BOND_LENGTH_TO_HYDROGEN;
                    } else {
                        maxDist = (
                                   atom1.getVanDerWaalsRadius()
                                   +
                                   atom2.getVanDerWaalsRadius()
                                  )
                                + (
                                   Constants.COORDINATE_UNCERTAINTY * 2
                                  );
                    }
                    if (dist <= maxDist) {
                       bonds.add(new Bond(atom1, atom2, BondTypes.SINGLE_BOND));
                    }
                }
            }
        }
        return bonds;
    }
    //--------------------------------------------------------------------------
    /**
     * Sort a list of PolyPeptide objects according to their sequence length.
     * @param peptides
     *        - List of PolyPeptide object to be sorted.
     */
    public static void sort(final ArrayList < PolyPeptide > peptides) {
        Collections.sort(peptides, new Comparator<PolyPeptide>() {
            public int compare(final PolyPeptide p1, final PolyPeptide p2) {
                if (p1.size() < p2.size()) {
                    return 1;
                }
                if (p1.size() > p2.size()) {
                    return -1;
                }
                return 0;
            } });
    }
}

