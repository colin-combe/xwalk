package structure.exceptions;

/**
 * This class is thrown by Reader classes if format of an input file is
 * unknown or erroneous.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class FileFormatException extends Exception {

    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param errorMessage
     *        - String object holding the error message, ideally some text that
     *          indicates the origin of the error.
     */
    public FileFormatException(final String errorMessage) {
        super(errorMessage);
    }
}
