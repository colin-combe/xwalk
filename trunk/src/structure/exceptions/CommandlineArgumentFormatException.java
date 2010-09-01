package structure.exceptions;

/**
 * Class that is thrown by the CommandlineArguments class if commandline
 * parameters are of wrong format.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class CommandlineArgumentFormatException extends Exception {

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
    public CommandlineArgumentFormatException(final String errorMessage) {
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
