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

package structure.matter.protein;

import structure.constants.Constants;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.Molecule;
import structure.matter.parameter.AminoAcidType;

/**
 * Class representing amino acid molecules.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class AminoAcid extends Molecule {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Type of this amino acid.
     */
    private AminoAcidType aminoAcidType;
    //--------------------------------------------------------------------------
    /**
     * Number of this amino acid within an amino acid sequence.
     */
    private int number;
    //--------------------------------------------------------------------------
    /**
     * Conservation grade of this amino acid.
     */
    private int conservationGrade;
    //--------------------------------------------------------------------------
    /**
     * Constructor. Checks whether all atoms really belong to a single amino
     * acid. Furthermore assigns the amino acid to one of various amino acid
     * types.
     * @param atoms
     *        - AtomList object holding the atomic coordinates of this amino
     *          acid.
     */
    public AminoAcid(final AtomList atoms) {
        super(atoms);
        // first check whether molecule really consists of a single amino acid.
        try {
            Atom preAtom = this.getAtom(0);
            for (Atom atom : this.getAllAtoms()) {
                if (!atom.getResidueName().equals(preAtom.getResidueName())
                    &&
                    atom.getResidueNumber() != preAtom.getResidueNumber()) {
                    throw new Exception("Data contains more than one residue "
                                      + "information");
                }
                if (!atom.getFlag().equals("ATOM  ")) {
                    throw new Exception("Data has non ATOM entries.");
                }
            }
        } catch (Exception e) {
            System.err.print("Error in reading in Residue information. " + e
                           + Constants.LINE_SEPERATOR);
        }

        // Determine the type/name of this amino acid
        for (AminoAcidType type : AminoAcidType.values()) {
            if (type.getThreeLetterCode().equals(
                                                this.getAtom(0).getResidueName()
                                                )) {
                this.aminoAcidType = type;
                break;
            }
        }
        this.number = atoms.get(0).getResidueNumber();
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this AminoAcid object.
     * @return Copy of this AminoAcid object.
     */
    public final AminoAcid copy() {
        AtomList atomListCopy = new AtomList();
        for (Atom atom : this.getAllAtoms()) {
            atomListCopy.add(atom.copy());
        }
        AminoAcid copy = new AminoAcid(atomListCopy);
        copy.number = this.number;
        copy.aminoAcidType = this.aminoAcidType;
        copy.conservationGrade = this.conservationGrade;

    return copy;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the type of this amino acid.
     * @return AminoAcidType object holding the type of this amino acid.
     */
    public final AminoAcidType getType() {
        return this.aminoAcidType;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the number residue number of this residue as defined in the PDB
     * file.
     * @return integer representing the residue number.
     */
    public final int getNumber() {
        return this.number;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of the amino acid's name, number and
     * chain Id separated each by a dash.
     * @param atom
     *        - Atom object being one of the atoms of an AminoAcid object or any
     *          other matter object.
     * @return String object holding the above mentioned amino acid information.
     */
    public static String getAminoAcidId(final Atom atom) {
        String residueId = atom.getResidueName().trim() + "-"
                         + atom.getResidueNumber();
        if (atom.getChainId() == ' ') {
            residueId += "-_";
        } else {
            residueId += "-" + atom.getChainId();
        }
    return residueId;
    }
    //--------------------------------------------------------------------------
    /**
     * Method to assign all atoms of a residue a rank position reflecting the
     * rank of the amino acid within the PDB sequence.
     * @param rank
     *        - Integer representing the rank position of this amino acid.
     */
    public final void setRank(final int rank) {
        for (Atom atom : this.getAllAtoms()) {
            atom.setRank(rank);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the conservation grade of this amino acid.
     * @param grade
     *        - integer value representing this amino acid's conservation grade.
     */
    public final void setConservationGrade(final int grade) {
        this.conservationGrade = grade;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the conservation grade of this amino acid.
     * @return integer value representing this amino acid's conservation grade.
     * @see #setConservationGrade(int)
     */
    public final int getConservationGrade() {
        return conservationGrade;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether this amino acid is likely to be observed at a
     * protein-protein interface.
     * @return {@code TRUE} if amino acid is ILE, HIS, CYS, PHE, MET, TYR or
     *         TRP
     */
    public final boolean isInterfaceTypic() {
        if (this.getType().equals(AminoAcidType.ISOLEUCINE)
         || this.getType().equals(AminoAcidType.HISTIDINE)
         || this.getType().equals(AminoAcidType.CYSTEINE)
         || this.getType().equals(AminoAcidType.PHENYLALANINE)
         || this.getType().equals(AminoAcidType.METHIONINE)
         || this.getType().equals(AminoAcidType.TYROSINE)
         || this.getType().equals(AminoAcidType.TRYPTOPHANE)) {
            return true;
        }
    return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether this amino acid is unlikely to be observed at a
     * protein-protein interface.
     * @return {@code TRUE} if amino acid is LYS, GLU, ASP, SER, THR, ALA, PRO
     *         or GLN
     */
    public final boolean isNonInterfaceTypic() {
        if (this.getType().equals(AminoAcidType.LYSINE)
         || this.getType().equals(AminoAcidType.GLUTAMIC_ACID)
         || this.getType().equals(AminoAcidType.ASPARTIC_ACID)
         || this.getType().equals(AminoAcidType.SERINE)
         || this.getType().equals(AminoAcidType.THREONINE)
         || this.getType().equals(AminoAcidType.ALANINE)
         || this.getType().equals(AminoAcidType.PROLINE)
         || this.getType().equals(AminoAcidType.GLUTAMINE)
        ) {
        return true;
        }
    return false;
    }
}
