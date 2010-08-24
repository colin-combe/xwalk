package xwalk.matter.hetgroups;

import xwalk.constants.Constants;
import xwalk.matter.Atom;
import xwalk.matter.AtomList;
import xwalk.matter.Molecule;

/**
 * Class representing any hetgroups/small molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 *
 */
public class SmallMolecule extends Molecule {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor. Checks whether all atoms really belong to a single small
     * molecule. Furthermore assigns the amino acid to one of various amino acid
     * types.
     * @param atoms
     *        - AtomList object holding the atomic coordinates of this amino
     *        acid.
     */
    public SmallMolecule(final AtomList atoms) {
        super(atoms);
        // first check whether molecule really consists of a single amino acid.
        try {
            Atom preAtom = this.getAtom(0);
            for (Atom atom : this.getAllAtoms()) {
                if (!atom.getResidueName().equals(preAtom.getResidueName())
                    &&
                    atom.getResidueNumber() != preAtom.getResidueNumber()) {
                    throw new Exception("Data contains more than one residue "
                                      + "information.");
                }
                if (!atom.getFlag().equals("HETATM")) {
                    throw new Exception("Data has non HETATM entries.");
                }
            }
        } catch (Exception e) {
            System.err.print("ERROR: in reading in SmallMolecule information. "
                           + e + Constants.LINE_SEPERATOR);
        }
    }
}
