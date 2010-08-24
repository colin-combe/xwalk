package xwalk.matter.pdb;

import java.util.ArrayList;
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
                                                 new ArrayList < AminoAcid > ();
    //--------------------------------------------------------------------------

    /**
     * List of small molecules that are associated to this polypeptides.
     */
    private ArrayList < SmallMolecule > smallMolecule =
                                             new ArrayList < SmallMolecule > ();
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
}
