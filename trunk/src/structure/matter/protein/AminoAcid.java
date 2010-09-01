package structure.matter.protein;

import structure.constants.Constants;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.Molecule;
import structure.matter.parameter.AminoAcidType;

/**
 * Class representing amino acid molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
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
}
