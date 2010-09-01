package structure.exceptions;

/**
 * Class that is thrown by the Mathematics class if inconsistencies occur during
 * transformation from Cartesian to spherical coordinates and vice versa.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class ConversionException extends Exception {

    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Error message.
     */
    private String message;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param errorMessage
     *        - String object holding the error message, ideally some text that
     *          indicates the origin of the error.
     */
    public ConversionException(final String errorMessage) {
        super(errorMessage);
    }

    //--------------------------------------------------------------------------
    /**
     * @return String object that holds the error message.
     */
    public final String toString() {
        return this.message;
    }
}
