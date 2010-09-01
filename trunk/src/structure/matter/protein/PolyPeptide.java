package structure.matter.protein;

import java.util.ArrayList;

import structure.constants.Constants;
import structure.matter.Atom;
import structure.matter.hetgroups.SmallMolecule;
/**
 * Class representing Proteins and polypeptides in general.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class PolyPeptide {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------

    /**
     * List of all amino acids in this polypeptide.
     */
    private ArrayList < AminoAcid > aminoAcidChain =
                                                 new ArrayList < AminoAcid >();
    //--------------------------------------------------------------------------

    /**
     * List of small molecules that are associated to this polypeptides.
     */
    private ArrayList < SmallMolecule > smallMolecule =
                                             new ArrayList < SmallMolecule >();
    //--------------------------------------------------------------------------

    /**
     * Constructor.
     * @param chain
     *        - List of AminoAcid object that this polypeptide object consists
     *          of.
     */
    public PolyPeptide(final ArrayList < AminoAcid > chain) {
        this.aminoAcidChain = chain;
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
        this.smallMolecule = hetgroups;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns all AminoAcid object that this polypeptide consists of.
     * @return List of AminoAcid objects.
     */
    public final ArrayList < AminoAcid > getAminoAcids() {
        return this.aminoAcidChain;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all SmallMolecule object that are associated to this polypeptide.
     * @return List of SmallMolecule objects.
     */
    public final ArrayList < SmallMolecule > getHetGroups() {
        return this.smallMolecule;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all PDB related information of all atoms in PDB format.
     * @return String object holding the text information of all atoms in PDB
     *         format.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (AminoAcid aa : this.getAminoAcids()) {
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
        for (AminoAcid a : this.getAminoAcids()) {
            buffer.append(a.getType().getOneLetterCode());
        }
        buffer.append(Constants.LINE_SEPERATOR);
    return buffer.toString();
    }

}

