package xwalk.constants;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import structure.matter.parameter.AminoAcidType;

/**
 * Class holding generic constant values related to Xwalk execution.
 * @author abdullah
 *
 */
public class Constants {
    /**
     * Constructor.
     */
    protected Constants() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    // CONSTANT NUMERIC VALUES
    //--------------------------------------------------------------------------
    /**
     * Maximum dimension of protein complex after which local grid
     * calculation is used.
     */
    public static final int MAX_PROTEIN_DIMENSION = 150;
    //--------------------------------------------------------------------------
    /**
     * Default cross-linker length.
     */
    public static final double DEFAULT_CROSS_LINKER_LENGTH = 27;
    //--------------------------------------------------------------------------
    /**
     * Minimum length of peptide in order to be detected by xQuest.
     */
    public static final int MIN_PEPTIDE_LENGTH = 5;
    //--------------------------------------------------------------------------
    /**
     * Maximum length of peptide in order to be detected by xQuest.
     */
    public static final int MAX_PEPTIDE_LENGTH = 40;
    //--------------------------------------------------------------------------
    /**
     * Number of miscleavages allowed in a digested peptide.
     */
    public static final int MAXIMUM_NUMBER_OF_MISCLEAVAGES = 1;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // CONSTANT STRING VALUES
    //--------------------------------------------------------------------------
    /**
     * Regular expression for checking peptide sequences for being
     * cross-linkable. Criteria for cross-linkable are:
     * <ol>
     *     <li>Tryptic peptide, i.e. C-terminus must be either arginine or
     *         lysine.
     *     </li>
     *     <li>One central lysine residue which represents the cross-linked
     *         amino acid.
     *     </li>
     *     <li>Up to one mis-cleavage is allowed. Here the miscleavage is
     *         positioned prior to the cross-linked lysine residue.
     *     </li>
     * </ol>
     */
    public static final String CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION1 =
         "[^"
        + AminoAcidType.LYSINE.getOneLetterCode()
        + AminoAcidType.ARGININE.getOneLetterCode()
        + "]*"
        + "["
        + AminoAcidType.LYSINE.getOneLetterCode()
        + AminoAcidType.ARGININE.getOneLetterCode()
        + "]{0,1}"
        + "[^"
        + AminoAcidType.LYSINE.getOneLetterCode()
        + AminoAcidType.ARGININE.getOneLetterCode()
        + "]*"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + "{1}"
       + "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + "["
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]$";
    //--------------------------------------------------------------------------
    /**
     * Regular expression for checking peptide sequences for being
     * cross-linkable. Criteria for cross-linkable are:
     * <ol>
     *     <li>Tryptic peptide, i.e. C-terminus must be either arginine or
     *         lysine.
     *     </li>
     *     <li>One central lysine residue which represents the cross-linked
     *         amino acid.
     *     </li>
     *     <li>Up to one mis-cleavage is allowed. Here the miscleavage is
     *         positioned post to the cross-linked lysine residue.
     *     </li>
     * </ol>
     */
    public static final String CROSS_LINKABLE_PEPTIDE_SEQUENCE_EXPRESSION2 =
         "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + "{1}"
       + "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + "["
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]{0,1}"
       + "[^"
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]*"
       + "["
       + AminoAcidType.LYSINE.getOneLetterCode()
       + AminoAcidType.ARGININE.getOneLetterCode()
       + "]$";
    //--------------------------------------------------------------------------
    // CONSTANT ENUM SETS
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // OTHER CONSTANT VALUES
    //--------------------------------------------------------------------------
    /**
     * Distances are given up to a single digit after comma.
     */
    public static final NumberFormat DISTANCE_DEC_FORMAT =
                                                       new DecimalFormat("0.0");

}
