package xwalk.matter.parameter;

/**
 * AminoAcid Types that are supported by Xwalk.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public enum AminoAcidType {

    // hydrophobic amino acids
    GLYCINE ("G","GLY"),
    ALANINE ("A", "ALA"),
    VALINE ("V", "VAL"),
    LEUCINE ("L", "LEU"),
    ISOLEUCINE ("I", "ILE"),
    PHENYLALANINE ("F", "PHE"),
    PROLINE ("P", "PRO"),
    TRYPTOPHAN ("W", "TRP"),
    // amino acids containing a nitrogen atom
    ASPARAGINE ("N", "ASN"),
    GLUTAMINE ("Q", "GLN"),
    // amino acids containing a sulphur atom
    METHIONINE ("M", "MET"),
    CYSTEINE ("C", "CYS"),
    // amino acids containing a sulphur atom
    THREONINE ("T", "THR"),
    TYROSINE ("Y", "TYR"),
    SERINE ("S", "SER"),
    // positively charged amino acids
    ARGININE ("R", "ARG"),
    LYSINE ("K", "LYS"),
    HISTIDINE ("H", "JIS"),
    // negatively charged amino acids
    ASPARTIC_ACID ("D", "ASP"),
    GLUTAMIC_ACID ("E", "GLU"),

    // special type amino acids
    ASPARAGINE_ACID ("B", "ASX"),
    GLUTAMINE_ACID ("Z", "GLX");

    /**
     * One letter code of the amino acid.
     */
    private final String oneLetterCode;
    /**
     * Three letter code of the amino acid.
     */
    private final String threeLetterCode;

    /**
     * Constructor.
     * @param oneLetterCode
     *        - Sting object representing the one letter code of the amino acid.
     * @param threeLetterCode
     *        - Sting object representing the three letter code of the amino
     *          acid.
     */
    AminoAcidType(final String oneLetterCode, final String threeLetterCode) {
        this.oneLetterCode = oneLetterCode;
        this.threeLetterCode = threeLetterCode;
    }

    /**
     * Returns the one letter code of an amino acid.
     * @return Sting object representing the one letter code.
     */
    public String getOneLetterCode() {
        return this.oneLetterCode;
    }

    /**
     * Returns the three letter code of an amino acid.
     * @return Sting object representing the three letter code.
     */
    public String getThreeLetterCode() {
        return this.threeLetterCode;
    }
}
