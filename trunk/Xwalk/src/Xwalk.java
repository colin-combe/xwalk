import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;
import java.util.zip.DataFormatException;

import structure.constants.Constants;
import structure.exceptions.CommandlineArgumentFormatException;
import structure.exceptions.CommandlineArgumentNotFoundException;
import structure.exceptions.FileFormatException;
import structure.matter.protein.PolyPeptideList;

import xwalk.crosslink.CrossLinkParameter;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.CrossLinkUtilities;
import xwalk.crosslink.CrossLinkParameter.Parameter;
import xwalk.io.CommandlineArguments;
import xwalk.io.DistanceWriter;


/**
 * Main class Xwalk that implements the main method to execute the virtual
 * cross-link calculation.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Xwalk {
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     */
    protected Xwalk() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Outputs on the STDERR channel that no cross-linker could be found in the
     * structure.
     * @param parameter
     *        - CrossLinkParameter object, holding all parameter that are
     *          necessary for the virtual cross-link calculation.
     */
    private static void outputNoXLfound(final CrossLinkParameter parameter) {
        String infile = parameter.getParameter(Parameter.INFILE_PATH).trim();
        String fileName = new File(infile).getName().trim();
        String residueType1 = parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME1
                                                    ).trim();
        String residueType2 = parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME2
                                                    ).trim();
        double maxDist = Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                  ).trim());

        System.err.println("WARNING: " + fileName + "\tThere is no pair of "
                           + residueType1.replaceAll("#", "") + "\t"
                           + residueType2.replaceAll("#", "")
                           + " residues that has a distance smaller than "
                           + maxDist + ".");

    }
    //--------------------------------------------------------------------------
    /**
     * Reads all user arguments from the commandline. Note System.exit commands
     * are executed if Exception occur during the read process of the
     * commandline.
     * @param args
     *        - Array of String object representing each one word on the
     *          commandline.
     * @return CommandlineArguments object holding all user set commandline
     *         parameter.
     */
    public static CommandlineArguments readCommandline(final String[] args) {

        String nl = Constants.LINE_SEPERATOR;

        CommandlineArguments arguments = null;
        try {
            arguments = new CommandlineArguments(args);
        } catch (FileNotFoundException e) {
            System.err.println(nl + "FileNotFoundException: " + e.toString());
            System.exit(-1);

        } catch (CommandlineArgumentNotFoundException e) {
            System.err.println(nl + "CommandlineArgumentNotFoundException: "
                             + e.toString());
            System.exit(-2);
        } catch (CommandlineArgumentFormatException e) {
            System.err.println(nl + "CommandlineArgumentFormatException: "
                             + e.toString());
            System.exit(-3);
        }
    return arguments;
    }
    //--------------------------------------------------------------------------
    /**
     * Creates virtual cross-links on the protein complexes given by the infile
     * commandline parameter.
     * @param parameter
     *        - CommandlineArguments object holding all user set commandline
     *          parameter.
     * @return List of CrossLink objects found on the infile protein complexes.
     */
    public static CrossLinkList createVirtualCrossLinks(
                                              final CrossLinkParameter parameter
                                                       ) {
        String nl = Constants.LINE_SEPERATOR;

        CrossLinkList list = null;
        try {
            // get all protein complex atom coordinates of the user given
            // inputFile.
            ArrayList < PolyPeptideList > complexes =
                                  CrossLinkUtilities.getComplexesCoordinates(
                                                                       parameter
                                                                            );

            list = CrossLinkUtilities.getVirtualCrossLinks(complexes,
                                                           parameter);
        } catch (FileNotFoundException e) {
            System.err.println(nl
                               + "ERROR: Infile could not be found" + nl + e
                               + nl);
            System.exit(-4);
        } catch (IOException e) {
            System.err.println(nl
                               + "ERROR: Could not read infile" + nl + e
                               + nl);
            System.exit(-5);
        } catch (FileFormatException e) {
            System.err.println(nl
                               + "ERROR: Format exception in input file" + nl
                               + e + nl);
            System.exit(-6);
        } catch (DataFormatException e) {
            System.err.println(nl
                               + "ERROR: GnuZip format exception in" + nl + e
                               + nl);
            System.exit(-7);
    }
    return list;
    }
    //--------------------------------------------------------------------------
    /**
     * Outputs all determined cross-links either on the terminal or into a file
     * depending on user's choice.
     * @param arguments
     *        - CommandlineArguments object holding all user set commandline
     *          parameter.
     * @param parameter
     *        - CommandlineArguments object holding all user set commandline
     *          parameter.
     * @param crossLinks
     *        - List of CrossLink objects found on the infile protein complexes.

     */
    public static void outputVirtualCrossLinks(
                                           final CommandlineArguments arguments,
                                           final CrossLinkParameter parameter,
                                           final CrossLinkList crossLinks
                                              ) {
        if (crossLinks.size() == 0) {
            Xwalk.outputNoXLfound(parameter);
        }

        if (arguments.getOutfileArgument().equals("")) {
            System.out.print(DistanceWriter.toString(crossLinks, parameter));
        } else {
            DistanceWriter write = new DistanceWriter();
            write.setFile(arguments.getOutfileArgument());
            if (arguments.isPymolOutputSet()) {
                write.writePymolScript(crossLinks, parameter);
            } else {
                write.writeFile(crossLinks, parameter);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Main program of Xwalk.
     * @param args -
     *        String array of commandline arguments
     */
    public static void main(final String[] args) {
        Locale.setDefault(Locale.US);

        if (args.length == 0) {
            CommandlineArguments.outputBasicHelpText();
        } else if (CommandlineArguments.isHelpSet(args)) {
            CommandlineArguments.outputVerboseHelpText();
        }

        CommandlineArguments arguments = Xwalk.readCommandline(args);
        CrossLinkParameter parameter = new CrossLinkParameter(arguments);
        // stop calculation if output is declined.
        if (!parameter.getParameter(Parameter.OUTFILE_PATH).equals("")
            &&
            !arguments.isOutputFileToBeCreated()) {
            System.exit(0);
        }
        CrossLinkList list = Xwalk.createVirtualCrossLinks(parameter);
        Xwalk.outputVirtualCrossLinks(arguments, parameter, list);
    }
}
