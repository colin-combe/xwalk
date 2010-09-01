package structure.matter;

import java.util.ArrayList;

import structure.matter.parameter.Element;



/**
 * A generic class holding important functions for classes that store a list
 * of Atom objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class AtomList extends ArrayList < Atom > {
    //--------------------------------------------------------------------------
    /**
     * default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Returns all PDB related information of all atoms in PDB format.
     * @return String object holding the text information of all atoms in PDB
     *         format.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (Atom atom : this) {
            output.append(atom.toString());
        }
    return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether this AtomList contains atom2, according to PDB
     * information.
     * @param atom2
     *        - Atom object holding the information of atom2.
     * @return {@code TRUE} if atom is in this AtomList, {@code FALSE}
     *         otherwise.
     */
    public final boolean contains(final Atom atom2) {
        for (Atom atom : this) {
            if (atom.equals(atom2)) {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Retrieves all hydrogen atoms from this list of atoms.
     * @return AtomList object with hydrogen atoms as elements. If none are
     *         found then the AtomList object will be empty.
     * @see #getHeavyAtoms()
     */
    public final AtomList getHydrogenAtoms() {
        AtomList hydrogens = new AtomList();
        for (Atom atom : this) {
             if (atom.getElement() == Element.HYDROGEN) {
                 hydrogens.add(atom);
            }
        }
        return hydrogens;
    }
    //--------------------------------------------------------------------------
    /**
     * Retrieves all heavy atoms, i.e. non-hydrogen atoms from this list of
     * atoms.
     * @return AtomList object with heavy atoms as elements. If none are found
     *         then the AtomList object will be empty.
     * @see #getHydrogenAtoms()
     */
    public final AtomList getHeavyAtoms() {
        AtomList heavyAtoms = new AtomList();
        for (Atom atom : this) {
             if (atom.getElement() == Element.HYDROGEN) {
                 heavyAtoms.add(atom);
            }
        }
        return heavyAtoms;
    }
    //--------------------------------------------------------------------------
}
