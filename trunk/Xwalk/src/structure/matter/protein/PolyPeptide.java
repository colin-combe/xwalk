package structure.matter.protein;

import java.util.ArrayList;

import structure.matter.Atom;
import structure.matter.hetgroups.SmallMolecule;
/**
 * Class representing Proteins and polypeptides in general.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class PolyPeptide extends ArrayList <AminoAcid> {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * List of small molecules that are associated to this polypeptides.
     */
    private ArrayList < SmallMolecule > smallMolecules =
                                             new ArrayList < SmallMolecule >();
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param chain
     *        - List of AminoAcid object that this polypeptide object consists
     *          of.
     */
    public PolyPeptide(final ArrayList < AminoAcid > chain) {
        this.addAll(chain);
    }
    //--------------------------------------------------------------------------

    /**
     * Adds small molecules to this polypeptide object.
     * @param hetgroups
     *        - List of SmallMolecule objects.
     */
    public final void addSmallMolecules(
                                     final ArrayList < SmallMolecule > hetgroups
                                       ) {
        this.smallMolecules = hetgroups;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all SmallMolecule object that are associated to this polypeptide.
     * @return List of SmallMolecule objects.
     */
    public final ArrayList < SmallMolecule > getSmallMolecules() {
        return this.smallMolecules;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all PDB related information of all atoms in PDB format.
     * @return String object holding the text information of all atoms in PDB
     *         format.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (AminoAcid aa : this) {
            for (Atom atom : aa.getAllAtoms()) {
                output.append(atom.toString());
            }
        }
    return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the one letter code sequence of this PolyPeptide object.
     * @return String object holding the one letter code.
     */
    public final String toStringOneLetterCode() {
        StringBuffer buffer = new StringBuffer();
        for (AminoAcid aa : this) {
            buffer.append(aa.getType().getOneLetterCode());
        }
    return buffer.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this PolyPeptide object.
     * @return Copy of this PolyPeptide object.
     */
    public final PolyPeptide copy() {
        ArrayList<AminoAcid> aminoAcidListCopy = new ArrayList<AminoAcid>();
        for (AminoAcid aa : this) {
            aminoAcidListCopy.add(aa.copy());
        }
        PolyPeptide copy = new PolyPeptide(aminoAcidListCopy);
        ArrayList<SmallMolecule> smallMolsCopy = new ArrayList<SmallMolecule>();
        for (SmallMolecule sm : this.getSmallMolecules()) {
            smallMolsCopy.add(sm.copy());
        }

        copy.smallMolecules = smallMolsCopy;
    return copy;
    }
}

