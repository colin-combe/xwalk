package xwalk.crosslink;

/**
 * List of CrossLinker types available to Xwalk.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 *
 */
public enum CrossLinkerType {

    /**
     * Available cross linker types.
     */
    DSS(25, 6);

    //--------------------------------------------------------------------------
    /**
     * Length of the cross-linker.
     */
    private double length;

    //--------------------------------------------------------------------------

    /**
     * Diameter of the cross-linker.
     */
    private double diameter;

    //--------------------------------------------------------------------------

    /**
     * Constructor.
     * @param xlLength
     *        - Length of the cross-linker
     * @param xlDiameter
     *        - Diameter of the cross-linker
     */
    CrossLinkerType(final double xlLength, final double xlDiameter) {
        this.length = xlLength;
        this.diameter = xlDiameter;
    }

    //--------------------------------------------------------------------------

    /**
     * Returns the length of the cross-linker.
     * @return double number representing the length of the cross-linker
     */
    public double getLength() {
        return this.length;
    }

    //--------------------------------------------------------------------------

    /**
     * Returns the diameter of the cross-linker.
     * @return double number representing the diameter of the cross-linker
     */
    public double getdiameter() {
        return this.diameter;
    }
}
