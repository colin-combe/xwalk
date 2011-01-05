package mm.constants;

import structure.matter.parameter.AminoAcidType;
/**
 * Various constants that are used by various classes.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Constants {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Constants() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    // CONSTANT NUMERIC VALUES
    //--------------------------------------------------------------------------
    /**
     * XlogP value of glycine residue, which serves as a reference for
     * calculating the XlogP value of all amino acids.
     */
    public static final double GLYCINE_LOGP = -1.128564;
    /**
     * XlogP value added to organic ions like COO- in order to correct for their
     * charge and thus higher solvation tendency.
     */
    public static final double ORGANIC_ION_CORRECTION_FACTOR = -1.083;
    /**
     * XlogP value given to all metal ions. Value is more or less arbitrary, but
     * should be lower than any charged amino acid.
     */
    public static final double METAL_XLOGP = -3;
    /**
     * Radius above which any physicochemical property calculations will be
     * omitted.
     */
    public static final double PHYSICOCHEMICAL_INFLUENCE_RADIUS = 9.0;

    //--------------------------------------------------------------------------
    // CONSTANT STRING VALUES
    //--------------------------------------------------------------------------
    /**
     * Text string identifying XlogP energy values in a hash.
     */
    public static final String XLOGP_ENERGY_PROPERTY_TAG = "XlogPenergy";
    /**
     * Text string identifying XlogP potential energy values in a hash.
     */
    public static final String XLOGP_POTENTIAL_ENERGY_PROPERTY_TAG =
                                                               "XlogPpotEnergy";

    //--------------------------------------------------------------------------
    // ENUM TYPES
    //--------------------------------------------------------------------------

    /**
     * Values taken from <a href="http://www.bioinformatics.cm-uj.krakow.pl/
     * reveal/">model</a>.
     */
    public enum PhysicoChemicalProperty {

        // hydrophobic amino acids
        GLYCINE(AminoAcidType.GLYCINE, 0.550, 0.00,  0.000),
        ALANINE(AminoAcidType.ALANINE, 0.572, 0.31, 0.570),
        VALINE(AminoAcidType.VALINE, 0.811, 1.22, 1.226),
        LEUCINE(AminoAcidType.LEUCINE ,0.783, 1.70, 2.191),
        ISOLEUCINE(AminoAcidType.ISOLEUCINE, 0.883, 1.80, 1.903),
        PHENYLALANINE(AminoAcidType.PHENYLALANINE, 0.906, 1.79, 1.970),
        PROLINE(AminoAcidType.PROLINE, 0.300, 0.72, 0.639),
        TRYPTOPHANE(AminoAcidType.TRYPTOPHANE, 0.856, 2.25, 2.892),
        // amino acids containing a nitrogen atom
        ASPARAGINE(AminoAcidType.ASPARAGINE, 0.278, -0.6, -1.961),
        GLUTAMINE(AminoAcidType.GLUTAMINE, 0.250, -0.22, -1.560),
        // amino acids containing a sulphur atom
        METHIONINE(AminoAcidType.METHIONINE, 0.828, 1.23, 1.558),
        CYSTEINE(AminoAcidType.CYSTEINE, 1.000, 1.54, 1.783),
        // amino acids containing a oxygen atom
        THREONINE(AminoAcidType.THREONINE, 0.478, 0.26, -0.583),
        TYROSINE(AminoAcidType.TYROSINE, 0.700, 0.96, 1.353),
        SERINE(AminoAcidType.SERINE, 0.422, -0.04, -1.099),
        // positively charged amino acids
        ARGININE(AminoAcidType.ARGININE, 0.272, -1.01, -0.873),
        LYSINE (AminoAcidType.LYSINE, 0.000, -0.99, -0.022),
        HISTIDINE(AminoAcidType.HISTIDINE, 0.628, 0.13, -0.009),
        // negatively charged amino acids
        GLUTAMIC_ACID(AminoAcidType.GLUTAMIC_ACID, 0.083, -0.64, -2.701),
        ASPARTIC_ACID(AminoAcidType.ASPARTIC_ACID, 0.167, -0.77, -3.102);

        /**
         * Fuzzy Oil Drop hydrophobicity value.
         */
        private double fuzzyOilDropHydrophobicity;
        /**
         * Fauchere and Pliska hydrophobicity value.
         */
        private double faucherePliskaHydrophobicity;
        /**
         * XlogP hydrophobicity value.
         */
        private double xlogPhydrophobicity;
        /**
         * Type of amino acid.
         */
        private AminoAcidType aminoAcidType;

        /**
         * Constructor.
         * @param type
         *        AminoAcidType object indicating which type of amino acid this
         *        aminoacid is.
         * @param fuzzyOilDrop
         *        double value representing the amino acids hydrophobicity
         *        value as defined by the
         *        <a href="http://www.bioinformatics.cm-uj.krakow.pl/reveal/">
         *        Fuzzy Oil Drop method</a>.
         * @param faucherePliska
         *        double value representing the amino acids hydrophobicity
         *        value as defined by Fauchere and Pliska, European Journal
         *        of Medicinal Chemistry (1983) vol. 18 (4) pp. 369-375.
         * @param xlogP
         *        double value representing the amino acids hydrophobicity
         *        value as defined by the XlogP method give in Kahraman et al.,
         *        Proteins, (2010).
         */
        PhysicoChemicalProperty(final AminoAcidType type,
                                final double fuzzyOilDrop,
                                final double faucherePliska,
                                final double xlogP) {
            this.fuzzyOilDropHydrophobicity = fuzzyOilDrop;
            this.faucherePliskaHydrophobicity = faucherePliska;
            this.xlogPhydrophobicity = xlogP;
            this.aminoAcidType = type;
        }
    }

}
