import java.io.FileNotFoundException;
import java.util.Locale;

import xwalk.constants.Constants;
import xwalk.crosslink.CrossLinkParameter;
import xwalk.crosslink.CrossLinkSet;
import xwalk.crosslink.CrossLinkUtilities;
import xwalk.exceptions.CommandlineArgumentFormatException;
import xwalk.exceptions.CommandlineArgumentNotFoundException;
import xwalk.io.CommandlineArguments;
import xwalk.io.WriteFile;


/**
 * Main class Xwalk that implements the main method to execute the virtual
 * cross-link calculation.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Xwalk {
    /**
     * Constructor.
     */
    protected Xwalk() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //-------------------------------------------------------------------------
    /**
     * Main program of Xwalk.
     * @param args -
     *        String array of commandline arguments
     */
    public static void main(final String[] args) {
        Locale.setDefault(Locale.US);

        String nl = Constants.LINE_SEPERATOR;

        if (args.length == 0) {
            CommandlineArguments.outputBasicHelpText();
        } else if (CommandlineArguments.isHelpSet(args)) {
            CommandlineArguments.outputVerboseHelpText();
        }

        CommandlineArguments arguments = null;
        try {
            arguments = new CommandlineArguments(args);
        } catch (FileNotFoundException e) {
            System.err.print(nl + "FileNotFoundException: " + e.toString());
        } catch (CommandlineArgumentNotFoundException e) {
            System.err.print(nl + "CommandlineArgumentNotFoundException: "
                             + e.toString());
        } catch (CommandlineArgumentFormatException e) {
            System.err.print(nl + "CommandlineArgumentFormatException: "
                             + e.toString());
        }

        CrossLinkParameter parameter = new CrossLinkParameter(arguments);

        CrossLinkSet set = null;
        try {
            set = CrossLinkUtilities.getVirtualCrossLinks(parameter);
        } catch (FileNotFoundException e) {
            System.err.println("ERROR: PDB file could not be found\n");
        }

        if (arguments.getOutfileArgument().equals("")) {
            if (arguments.isPymolOutputSet()) {
                System.out.print(CrossLinkUtilities.outputPymolScript(set,
                                                                      parameter)
                                );
            } else {
                System.out.print(CrossLinkUtilities.toString(set, parameter));
            }
        } else {
            WriteFile write = new WriteFile();
            write.setFile(arguments.getOutfileArgument());
            if (arguments.isPymolOutputSet()) {
                write.write(CrossLinkUtilities.outputPymolScript(set,
                                                                 parameter)
                           );
            } else {
                write.write(CrossLinkUtilities.toString(set, parameter));
            }
        }
    }
}
