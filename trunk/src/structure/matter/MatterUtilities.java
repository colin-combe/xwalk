package structure.matter;

import java.util.ArrayList;

import structure.constants.Constants;
import structure.constants.Constants.BondTypes;
import structure.constants.Constants.ElementTypes;
import structure.math.Mathematics;
import structure.math.Point3d;
import structure.matter.parameter.Element;
import structure.matter.pdb.ProteinComplex;


/**
 * A generic class that holds various methods to operate on classes defined in
 * the matter package.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
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
     *        - ProteinComplex object to which dimension should be calculated.
     * @return double value representing the dimension of the protein complex.
     */
    public static double getDimension(final ProteinComplex complex) {
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
        // 1.54 Å between both atom centers.
        // (see http://en.wikipedia.org/wiki/Bond_length).
        // To account for resolution error distance is taken to be less than
        // 1.8 (0.28 Å estimated standard error for X-ray structures)
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
                        maxDist = Constants.BOND_TO_HYDROGEN;
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
}

