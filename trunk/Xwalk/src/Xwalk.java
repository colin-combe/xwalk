/*
 * (C) 2010 Abdullah Kahraman
 *
 * This software is part of the open-source project "Xwalk". You can use this
 * software under the terms of the
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 * (http://creativecommons.org/licenses/by-nc-sa/3.0/).
 * This means that you
 * 1.) can copy, modify, distribute the software
 * 2.) must give credit to the author
 * 3.) must not use this work for commercial purposes
 * 4.) must license derivative works under the same or a similar license.
 */

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

import xwalk.crosslink.CrossLink;
import xwalk.crosslink.CrossLinkParameter;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.CrossLinkUtilities;
import xwalk.crosslink.MonoLink;
import xwalk.crosslink.MonoLinkList;
import xwalk.crosslink.CrossLinkParameter.Parameter;
import xwalk.io.CommandlineArguments;
import xwalk.io.DistanceWriter;


/**
 * Main class Xwalk that implements the main method to execute the virtual
 * cross-link calculation.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
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
// Residue information is not displayed any more. No need for this.
/*        String residueType1 = parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME1
                                                    ).trim();
        String residueType2 = parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME2
                                                    ).trim();
        double maxDist = Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                  ).trim());
*/
        System.err.println("WARNING: " + fileName + "\tNo virtual cross-links "
                           + "found.");

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
     * Returns a list of monoLinks for the protein complexes given by the infile
     * commandline parameter.
     * @param parameter
     *        - CommandlineArguments object holding all user set commandline
     *          parameter.
     * @return List of MonoLink objects found on the infile protein complexes.
     */
    public static MonoLinkList getMonoLinks(
                                              final CrossLinkParameter parameter
                                                       ) {
        String nl = Constants.LINE_SEPERATOR;

        MonoLinkList list = null;
        try {
            // get all protein complex atom coordinates of the user given
            // inputFile.
            ArrayList < PolyPeptideList > complexes =
                                  CrossLinkUtilities.getComplexesCoordinates(
                                                                       parameter
                                                                            );

            list = CrossLinkUtilities.getVirtualMonoLinks(complexes,
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
     * @param parameter
     *        - CommandlineArguments object holding all user set commandline
     *          parameter.
     * @param crossLinks
     *        - List of CrossLink objects found on the infile protein complexes.
     * @param monoLinks
     *        - List of MonoLink objects found on the infile protein complexes.
     */
    public static void outputVirtualCrossLinks(
                                           final CrossLinkParameter parameter,
                                           final CrossLinkList crossLinks,
                                           final MonoLinkList monoLinks
                                              ) {
        if (crossLinks.size() == 0) {
            Xwalk.outputNoXLfound(parameter);
        }

        if (parameter.getParameter(Parameter.OUTFILE_PATH).equals("")) {
            System.out.print(DistanceWriter.toString(crossLinks,
                                                     monoLinks,
                                                     parameter));
        } else {
            DistanceWriter write = new DistanceWriter();
            write.setFile(parameter.getParameter(Parameter.OUTFILE_PATH));
            if (Boolean.parseBoolean(parameter.getParameter(
                                                      Parameter.DO_PYMOL_OUTPUT)
                                                           )) {
                System.out.print(DistanceWriter.toString(crossLinks,
                                                         monoLinks,
                                                         parameter)
                                                        );
                write.writePymolScript(crossLinks,
                                       monoLinks,
                                       parameter);
            } else {
                write.writeFile(crossLinks,
                                monoLinks,
                                parameter);
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
        CrossLinkList xlList = Xwalk.createVirtualCrossLinks(parameter);
        MonoLinkList monoList = new MonoLinkList();
        if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_MONO_CROSSLINK
                                                       ))) {
            monoList = Xwalk.getMonoLinks(parameter);
            // show only mono-links for residues that have not been found in
            // cross-links.
            for (CrossLink xl : xlList) {
                for (MonoLink ml : monoList) {
                    if (xl.isInBond(ml)) {
                        monoList.remove(ml);
                        break;
                    }
                }
            }
            // in order to set index of monoList, first find out max index of
            // cross-links.
            int maxIndex = 0;
            for (CrossLink xl : xlList) {
                maxIndex = Math.max(maxIndex, xl.getIndex());
            }
            if (parameter.getParameter(
                                     Parameter.DISTANCE_FILE_PATH).equals("")) {
                for (MonoLink ml : monoList) {
                    ml.setIndex(++maxIndex);
                }
            }
        }
        Xwalk.outputVirtualCrossLinks(parameter, xlList, monoList);
    }
}
