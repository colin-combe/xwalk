package xwalk.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Class for handling commandline arguments.
 * @author Abdullah Kahraman
 * @version 3.0
 */
public class Commandline {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Commandline() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Check and get a commandline argument with its value. The argument must
     * not be fully spelled, as long as it matches unambiguously the beginning
     * of an argument on the commandline.
     * @param  args
     *         String array of commandline arguments. Normally args[] array from
     *         the main-method.
     * @param  param
     *         String object of argument to look for, within ARGS[]
     * @param  hasValue
     *         if {@code TRUE} than argument has a value which will be returned.
     * @return String with value of parameter
     *         "EXISTS" if argument exists and hasValue is {@code FALSE}, or
     *         "AMBIGUOUS" if argument can not be determined unambiguously
     *                     from the commandline, or
     *         "ERROR" if argument does not exists.
     * @see    #get()
     */
    public static String get(final String[] args, final String param,
                             final boolean hasValue) {
        int c = 0;
        for (String argument : args) {
             if (param.toLowerCase().startsWith(argument.toLowerCase())) {
                c++;
            }
        }
        if (c > 1) {
            return "AMBIGUOUS";
        }

        // loop over all arguments in the commandline
        for (int i = 0; i < args.length; i++) {
             if (param.toLowerCase().startsWith(args[i].toLowerCase()) 
                 &&
                 !hasValue) {
                 return "EXISTS";
             }
             if (param.toLowerCase().startsWith(args[i].toLowerCase())
                 &&
                 hasValue
                 &&
                 args.length > i + 1) {
                 return args[i + 1];
             }
        }
    return "ERROR";
    }

    //--------------------------------------------------------------------------
    /**
     * Returns a user input from a dialog on the commandline.
     * @return String object that holds the answer of the user to the dialog
     * question. If problems occur during the read in of the user input an empty
     * string object "" is returned.
     */
    public static String get() {
        BufferedReader br = new BufferedReader(new InputStreamReader(
                                                                     System.in
                                                                     )
                                              );
        try {
            return br.readLine();
        } catch (IOException ioe) {
            return "";
        }
    }
}
