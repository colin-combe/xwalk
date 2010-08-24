package xwalk.matter.pdb;

import xwalk.matter.AtomList;
import xwalk.matter.Molecule;

/**
 * Dummy class as a place holder for future implementations on small molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 *
 */
public class SmallMolecule extends Molecule {

    //--------------------------------------------------------------------------
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atoms
     *        - List of Atom objects that this small molecule consists of.
     */
    public SmallMolecule(final AtomList atoms) {
        super(atoms);
    }
    //--------------------------------------------------------------------------
}
