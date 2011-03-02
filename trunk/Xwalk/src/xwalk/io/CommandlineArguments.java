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
 *
 */

package xwalk.io;
import java.io.FileNotFoundException;

import structure.constants.Constants;
import structure.exceptions.CommandlineArgumentFormatException;
import structure.exceptions.CommandlineArgumentNotFoundException;
import structure.io.Commandline;
import structure.io.ReadFile;


/**
 * Class that holds attributes and provides methods for the handling of
 * commandline arguments for the Xwalk class.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 *
 */
public class CommandlineArguments {
    /**
     * Operating system independent new line character.
     */
    private static final String NL = Constants.LINE_SEPERATOR;
    /**
     * String array holing all commandline arguments.
     */
    private String[] arguments;
    /**
     * Path to the input file.
     * Default {@code infile = ""}.
     */
    private String infile = "";
    /**
     * Path to a distance input file.
     * Default {@code distanceInfile = ""}.
     */
    private String distanceInfile = "";
    /**
     * To read in only atoms in the backbone of the atom plus beta-carbon.
     */
    private boolean doBackboneReadOnly = false;
    /**
     * Output only when all distances in a distance file are found.
     * Default {@code findAll = FALSE};
     */
    private boolean findAll = false;
    /**
     * Path to the output file.
     * Default {@code outfile = ""}.
     */
    private String outfile = "";
    /**
     * To determine whether output file is to be created.
     * Default {@code doOutputFile = false}.
     */
    private boolean doOutputFile = false;
    /**
     * Force output into a file without checking whether file already exists.
     * Default {@code force = FALSE}.
     */
    private boolean force = false;
    /**
     * Three letter PDB residue code of 1st type of cross-linked amino acids.
     * Default {@code residueType1 = ""}.
     */
    private String residueType1 = "";
    /**
     * Three letter PDB residue code of 2nd type of cross-linked amino acids.
     * Default {@code residueType2 = ""}.
     */
    private String residueType2 = "";
    /**
     * PDB residue number of 1st type of cross-linked amino acids.
     * Default {@code aa1resNo = -999}.
     */
    private String aa1resNo = "-999";
    /**
     * PDB residue number of 2nd type of cross-linked amino acids.
     * Default {@code aa2resNo = -999}.
     */
    private String aa2resNo = "-999";
    /**
     * PDB chain Id of 1st type of cross-linked amino acids.
     * Default {@code Constants.ALPHANUMERIC}.
     */
    private String chainIds1 = Constants.ALPHANUMERIC;
    /**
     * PDB chain Id of 2nd type of cross-linked amino acids.
     * Default {@code Constants.ALPHANUMERIC}.
     */
    private String chainIds2 = Constants.ALPHANUMERIC;
    /**
     * PDB alternative location Id of 1st type of cross-linked amino acids.
     * Default {@code Constants.ALPHANUMERIC}.
     */
    private String altLocs1 = Constants.ALPHANUMERIC;
    /**
     * PDB alternative location Id of 2nd type of cross-linked amino acids.
     * Default {@code Constants.ALPHANUMERIC}.
     */
    private String altLocs2 = Constants.ALPHANUMERIC;
    /**
     * PDB atom name of 1st cross-linked atom.
     * Default {@code atomType1 = ""}.
     */
    private String atomType1 = "";
    /**
     * PDB atom name of 2nd cross-linked atom.
     * Default {@code atomType2 = ""}.
     */
    private String atomType2 = "";
    /**
     * Maximal distance that is span by the cross-linker.
     * Default {@code maxDist =
     * xwalk.constants.Constants.DEFAULT_CROSS_LINKER_LENGTH}.
     */
    private double maximumDistance =
                          xwalk.constants.Constants.DEFAULT_CROSS_LINKER_LENGTH;
    /**
     * To output a PyMol script that visualizes all cross-links.
     * Default {@code pymolOutput = FALSE}.
     */
    private boolean pymolOutput = false;
    /**
     * To calculate the cross-link associated Solvent-Path-Distance.
     * Default {@code solvDist = TRUE};
     */
    private boolean solventPathDistance = true;
    /**
     * Length of a cubic grid cell edge. Grids are used to calculate
     * Solvent-Path-Distances.
     * Default {@code gridCellLength = Constants.DEFAULT_GRID_CELL_SIZE}.
     */
    private double gridCellLength = Constants.DEFAULT_GRID_CELL_SIZE;
    /**
     * To calculate solvent accessible surface areas and allow only cross-links
     * between residues that are accessible to a solvent molecule with radius
     * solventRadius.
     * Default {@code doSas = TRUE};
     */
    private boolean doSAS = true;
    /**
     * Solvent radius for calculating SAS.
     * Default {@code solventRadius = 1.4}.
     */
    private double solventRadius = Constants.SOLVENT_RADIUS;
    /**
     * To regard the protein complex as a homomeric one, which disregards
     * cross-links that are formed between identical numbered and typed amino
     * acids in different chains.
     * Default {@code homo = FALSE}.
     */
    private boolean homo = false;
    /**
     * To output only intramolecular cross-link distances.
     * Default {@code showIntra = TRUE};
     */
    private boolean showIntra = true;
    /**
     * To output only intermolecular cross-link distances.
     * Default {@code showInter = TRUE};
     */
    private    boolean showInter         = true;
    /**
     * To output some information about the current execution status of the
     * program.
     * Default {@code verbose = FALSE};
     */
    private boolean verbose         = false;
    /**
     * Use local grids for Solvent-Path distance calculations.
     * Default {@code verboseGrid = FALSE};
     */
    private boolean useGlobalGrid = false;
    /**
     * To output all grids that are used to calculate the
     * Solvent-Path distances.
     * Default {@code verboseGrid = FALSE};
     */
    private boolean verboseGrid     = false;
    /**
     * To trypsinate the protein.
     * Default {@code help = FALSE};
     */
    private boolean doTrypsin = false;
    /**
     * To digest according to the exceptions listed in
     * <a href="http://expasy.org/tools/peptidecutter/
     * peptidecutter_enzymes.html">ExPASy</a>.
     */
    private boolean doExpasy = false;
    /**
     * To output probabilities according to experimental data on DSS and BS3.
     */
    private boolean doProbability = false;
    /**
     * To output the help text that includes all commandline arguments of Xwalk.
     * Default {@code help = FALSE};
     */
    private static boolean help            = false;

    //--------------------------------------------------------------------------
    /**
     * Constructor, which is determining all arguments that are set by the user.
     * @param  args
     *         Array of String objects holding all commandline arguments.
     * @throws FileNotFoundException if the input file set by the user could not
     *         be found.
     * @throws CommandlineArgumentNotFoundException <br>
     *         - if input file path is not set by the user <br>
     *         - if residue informations are not set through -r1,r2 or -aa1,aa2
     *           <br>
     *         - if flag -out is set without specifying the output file path.
     * @throws CommandlineArgumentFormatException
     *         if output file and PyMOL script is set but output file suffix is
     *         not equal to {@code .pml},
     */
    public CommandlineArguments(final String[] args) throws
                                           FileNotFoundException,
                                           CommandlineArgumentNotFoundException,
                                           CommandlineArgumentFormatException {
        this.arguments = args;

        this.readAlternativeLocation1Argument();
        this.readAlternativeLocation2Argument();
        this.readAtomType1Argument();
        this.readAtomType2Argument();
        this.readChainIds1Argument();
        this.readChainIds2Argument();
        this.readForceArgument();
        this.readGridCellSizeArgument();
        this.readGridOutputArgument();
        this.readHomomericArgument();
        this.readInfileArgument();
        this.readBackBoneOnlyArgument();
        this.readDistanceInfileArgument();
        this.readFindAllArgument();
        this.readInterMolecularDistanceArgument();
        this.readIntraMolecularDistanceArgument();
        this.readMaximumDistanceArgument();
        this.readSolventAccessibiltyArgument();
        this.doOutputFile = this.readOutfileArgument();
        this.readPymolArgument();
        this.readSolventRadiusArgument();
        this.readAminoAcidNumber1Argument();
        this.readAminoAcidNumber2Argument();
        this.readAminoAcidName1Argument();
        this.readAminoAcidName2Argument();
        this.readSolventPathDistanceArgument();
        this.readGlobalGridArgument();
        this.readTrypsinateArgument();
        this.readProbabilityArgument();
        this.readExpasyArgument();
        this.readVerboseOutputArgument();
    }
    //--------------------------------------------------------------------------

    /**
     * Returns a comprehensive list of in information about Xwalk and its
     * commandline arguments. These information include:
     * <li> an exemplary execution command,
     * <li> version number,
     * <li> about text,
     * <li> output format and
     * <li> a list of commandline arguments.
       * @return String object holding the comprehensive list of information
     */
    private static String getVerboseHelpText() {
        return NL
              + "EXAMPLARY command for program execution:" + NL
              + "Xwalk -in 1brs.pdb -aa1 LYS#ARG -aa2 lys#arg -max 21"
              + NL
              + NL
              + "ABOUT"
              + NL
              + "Version 0.1"
              + NL
              + "Xwalk calculates and outputs distances in Angstroem "
              + "for potential cross-links "
              + NL
              + "between -aa1 type amino acids and -aa2 type amino "
              + "acids in the PDB file -in."
              + NL
              + NL
              + "IMPORTANT"
              + NL
              + "If large protein complexes are processed, the Java "
              + "heap size might need to be"
              + NL
              + "increased from the default 64MB to 256MB, with the "
              + "Java parameter -Xmx256m "
              + NL
              + NL
              + "OUTPUT FORMAT:"
              + NL
              + "IndexNo\tInfileName\tAtom1info\tAtom2info\t"
              + "DistanceInPDBsequence\tEuclideanDistance\t"
              + "SolventPathDistance\t(EucProbability\tSASDprobability)\t"
              + "PeptidePairSequences"
              + NL
              + NL
              + "Virtual cross-links are sorted first by "
              + "decreasing probability, then by increasing SASD and "
              + "finally by increasing Euclidean distance."
              + NL
              + NL
              + "Commandline PARAMETER:"
              + NL
              + "INPUT/OUTPUT:"
              + NL
              + "\t-infile\t<path>\tAny PDB file; .tar, .gz and .tar.gz files "
              + "with PDB file content are also accepted [required]."
              + NL
              + "\t-bb\t[switch]\tReads in only backbone and beta carbon "
              + "atom coordinates from the input file and sets -radius to 2.0. "
              + "This might be of value when virtual cross-links are to be "
              + "created between backbone or beta-carbon atoms [optional]."
              + NL
              + "\t-dist\t<path>\tAny Xwalk distance file, from which all "
              + "residue information will be extracted [optional]."
              + NL
              + "\t-out\t<path>\tWrites output to this file, otherwise "
              + "output is directed to the STDOUT channel. If -pymol is "
              + "set than filename must have .pml filename ending "
              + "[optional]."
              + NL
              + "\t-f\t[switch]\tForces output to be written into a "
              + "file even if file already exists [optional]."
              + NL
              + "\tpymol\t[switch]\tOutputs a PyMOL (http://www.pymol."
              + "org/) script highlighting the calculated distances of "
              + "the potential cross-links [optional]."
              + NL
              + "\t-v\t[switch]\tOutputs various information other "
              + "than distances [optional]."
              + NL
              + "\t-grid\t[switch]\tOutputs on STDOUT channel the grid, which"
              + "is used to calculate the Solvent Accessible Surface Distance. "
              + "Requires the option -global. The grid is in PDB format with "
              + "distances in the B-factor column [optional]."
              + NL
              + NL
              + "RESIDUE/ATOM SELECTION:"
              + NL
              + "\t-aa1\t[String]\tThree letter code of 1st amino "
              + "acid. To specify more than one amino acid use '#' as "
              + "a delimeter [required, if -r1 is not set]."
              + NL
              + "\t-aa2\t[String]\tThree letter code of 2nd amino "
              + "acid. To specify more than one amino acid use '#' as "
              + "a delimeter [required, if -r2 is not set]."
              + NL
              + "\t-r1\t[String]\tAmino acid residue number. To "
              + "specify more than one residue number use '#' as a "
              + "delimeter. [required, if -aa1 is not set]."
              + NL
              + "\t-r2\t[String]\tAmino acid residue number. To "
              + "specify more than one residue number use '#' as a "
              + "delimeter. [required, if -aa2 is not set]."
              + NL
              + "\t-c1\t[String]\tChain ids for -aa1 or -r1. For blank "
              + "chain Id use '_'. To specify more than one chain Id, "
              + "append chain ids to a single string, e.g. ABC "
              + "[optional](default: all chain Ids)."
              + NL
              + "\t-c2\t[String]\tChain ids for -aa2 or -r2. For blank "
              + "chain Id use '_'. To specify more than one chain Id, "
              + "append chain ids to a single string, e.g. ABC "
              + "[optional](default: all chain Ids)."
              + NL
              + "\t-a1\t[String]\tAtom type for -aa1 or -r1. To "
              + "specify more than one atom type use '#' as a "
              + "delimeter. [optional]."
              + NL
              + "\t-a2\t[String]\tAtom type for -aa2 or -r2. To "
              + "specify more than one atom type use '#' as a "
              + "delimeter. [optional]."
              + NL
              + "\t-l1\t[String]\tAlternative location id for -aa1 or "
              + "-r1. To specify more than one alternative location, "
              + "append alternative location ids to a single string, "
              + "e.g. AB [optional]."
              + NL
              + "\t-l2\t[String]\tAlternative location id for -aa2 or "
              + "-r1. To specify more than one alternative location, "
              + "append alternative location ids to a single string, "
              + "e.g. AB [optional]."
              + NL
              + "\t-intra\t[switch]\tOutputs only \"intra-molecular\" "
              + "distances [optional]."
              + NL
              + "\t-inter\t[switch]\tOutputs only \"inter-molecular\" "
              + "distances [optional]."
              + NL
              + "\t-homo\t[double]\tOutputs only shortest distance of "
              + "potential cross-links between equally numbered "
              + "residues. Reduces redundancy if PDB file is a "
              + "homomeric protein complex. [optional]."
              + NL
              + NL
              + "DIGESTION RELATED:"
              + NL
              + "\t-trypsin\t[switch]\tDigests in silico the protein with "
              + "trypsin and excludes peptides that are shorter "
              + "than "
              + xwalk.constants.Constants.MIN_PEPTIDE_LENGTH
              + " AA or larger than "
              + xwalk.constants.Constants.MAX_PEPTIDE_LENGTH
              + " AA [optional]."
              + NL
              + NL
              + "DISTANCE RELATED:"
              + NL
              + "\t-max\t[double]\tCalculates distances in Angstroem "
              + "only up-to this value, where the value must be smaller than "
              + xwalk.constants.Constants.MAX_SASD_DISTANCE
              + " for SASD calculations. (default: "
              + xwalk.constants.Constants.DEFAULT_CROSS_LINKER_LENGTH + ")."
              + NL
              + "\t-euc\t[switch]\tSkips Solvent-Path-Distance "
              + "calculation and outputs only Euclidean distances "
              + "[optional]. "
              + NL
              + "\t-prob\t[switch]\tOutputs probability information for each "
              + "vXL as determined by experimental data on DSS and BS3 cross-"
              + "linking experiments [optional]. "
              + NL
              + "\t-all\t[switch]\tOutputs only distances if all virtual cross-"
              + "links in a distance file (-dist) can be found."
              + NL
              + NL
              + "SOLVENT-PATH-DISTANCE GRID RELATED:"
              + NL
              + "\t-global\t[switch]\tUses a single large grid to calculate "
              + "all Solvent-Path distances, rather than using a smaller local "
              + "grid around each cross-linkable atom. This option will "
              + "reduce computation time for smaller protein complexes. Note "
              + "however that the distances might differ by about an Angstroem "
              + "as compared to a local grid and that you might need to "
              + "increase the Java heap size to for example 512 MB (-Xmx512m)"
              + NL
// No need for -xsas parameter anymore, as SASD makes no sense without checking
// for SAS.
//            + "\t-xsas\t[switch]\tDoes not calculate the solvent accessible "
//            + "surface surface area and thus does not exclude non-accessible "
//            + "amino acids [optional]."
//            + NL
              + "\t-radius\t[double]\tSolvent radius for calculating the "
              + "solvent accessible surface area [optional](default "
              + Constants.SOLVENT_RADIUS
              + ")."
              + NL
              + "\t-space\t[double]\tSpacing in Angstroem between grid "
              + "cells. [optional](default "
              + Constants.DEFAULT_GRID_CELL_SIZE + ")."
              + NL
              + NL;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the information that with -help/-h a list of commandline
     * arguments can be retrieved.
     * @return String object holding the information.
     * @see #outputBasicHelpText()
     */
    private static String getBasicHelpText() {
        return NL
               + "Please type -h or -help for a list of commandline arguments"
               + NL;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns whether -help/-h flag has been set on the commandline.
     * @param args
     *        - String array object holding all arguments typed by the user on
     *          the commandline.
     * @return {@code TRUE} if the user has set the -help/-h flag, {@code FALSE}
     *         otherwise.
     * @see #readHelpArgument(String[])
     * @see #getVerboseHelpText()
     * @see #readVerboseOutputArgument()
     */
    public static boolean isHelpSet(final String[] args) {
        CommandlineArguments.readHelpArgument(args);
        return help;
    }
    //--------------------------------------------------------------------------

    /**
     * Outputs on the STDERR channel a comprehensive list of in information
     * about Xwalk and its commandline arguments.
     * Be aware that after the output a System.exit(0) signal is sent.
     * @see #readHelpArgument(String[])
     * @see #isHelpSet(String[])
     * @see #getVerboseHelpText()
     */
    public static void outputVerboseHelpText() {
        System.err.print(CommandlineArguments.getVerboseHelpText());
        System.exit(0);
    }
    //--------------------------------------------------------------------------

    /**
     * Outputs on the STDERR channel the information that with -help/-h a list
     * of commandline arguments can be retrieved.
     * @see #getBasicHelpText()
     * @see #outputVerboseHelpText()
     */
    public static void outputBasicHelpText() {
        System.err.println(CommandlineArguments.getBasicHelpText());
        System.exit(0);
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the flag -f has been set on the commandline.
     * @see #isForceOutputSet()
     */
    private void readForceArgument() {
        if (Commandline.get(this.arguments, "-f", false).equals("EXISTS")) {
            this.force = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether the output should be forced to be written out,
     * even if file already exists.
     * @return {@code TRUE} if forced output should be performed, {@code FALSE}
     * otherwise.
     * @see #readForceArgument()
     */
    public final boolean isForceOutputSet() {
        return this.force;
    }

    //--------------------------------------------------------------------------
    /**
     * Determines whether the flag -h or -help has been set on the commandline.
     * @param args
     *        - String array object holding all arguments typed by the user on
     *          the commandline.
     * @see #isHelpSet()
     * @see #outputVerboseHelpText()
     * @see #getVerboseHelpText()
     */
    private static void readHelpArgument(final String[] args) {
        if (Commandline.get(args, "-help", false).equals("EXISTS")
            ||
            Commandline.get(args, "-h", false).equals("EXISTS")
           ) {
            CommandlineArguments.help = true;
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Determines whether the argument -infile has been set on the commandline.
     * @throws FileNotFoundException if file could not be found
     * @throws CommandlineArgumentNotFoundException if required commandline
     *         argument could not be found.
     * @see #getInfileArgument()
     */
    private void readInfileArgument() throws
                                          FileNotFoundException,
                                          CommandlineArgumentNotFoundException {
        if (Commandline.get(this.arguments, "-infile", true).equals("ERROR")) {
            throw new CommandlineArgumentNotFoundException(NL
                                                          + "ERROR: Could NOT "
                                                          + "find value for "
                                                          + "parameter \""
                                                          + "-infile\"."
                                                          );
        } else {
            this.infile = Commandline.get(this.arguments,
                                          "-infile",
                                          true
                                         ).trim();
            if (!ReadFile.exists(this.infile)) {
                throw new FileNotFoundException(NL
                                                + "ERROR: File \"" + infile
                                                + "\" NOT found!!!"
                                                + NL
                                                + NL
                                               );
            }
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the path to the input file.
     * @return String object holding the path to the input file.
     * @see #readInfileArgument()
     */
    public final String getInfileArgument() {
        return this.infile;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -dist has been set on the commandline.
     * @throws FileNotFoundException if file could not be found
     * @see #getInfileArgument()
     */
    private void readDistanceInfileArgument() throws
                                          FileNotFoundException {
        if (!Commandline.get(this.arguments, "-dist", true).equals("ERROR")) {
            this.distanceInfile = Commandline.get(this.arguments, "-dist", true
                                                                       ).trim();
            if (!ReadFile.exists(this.distanceInfile)) {
                throw new FileNotFoundException(NL
                                                + "ERROR: File \""
                                                + distanceInfile
                                                + "\" NOT found!!!"
                                                + NL
                                                + NL
                                               );
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to a distance file.
     * @return String object holding the path to the distance file.
     * @see #readInfileArgument()
     */
    public final String getDistanceInfileArgument() {
        return this.distanceInfile;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -all has been set on the commandline.
     * @see #isFindAllSet()
     */
    private void readFindAllArgument() {
        if (Commandline.get(this.arguments, "-all", false).equals("EXISTS")) {
            if (!Commandline.get(
                                  this.arguments, "-dist", true).equals("ERROR")
                                ) {
                this.findAll = true;
            } else {
                System.err.println(NL + "WARNING: Ommiting -all argument. "
                                 + "Please specify -dist argument otherwise."
                                 + NL);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether the all distances in a distance file must be found prior
     * to their output.
     * @return {@code TRUE} if all distances must be found, {@code FALSE}
     * otherwise.
     * @see #readFindAllArgument()
     */
    public final boolean isFindAllSet() {
        return this.findAll;
    }
    //--------------------------------------------------------------------------

    /**
     * Determines whether the argument -pymol has been set on the commandline.
     * Also checks that in case an output file is set, that its suffix is
     * {@code .pml}.
     * @throws CommandlineArgumentFormatException if output file is set but
     *         suffix is not equal to {@code .pml}.
     * @see #isPymolOutputSet()
     */
    private void readPymolArgument() throws CommandlineArgumentFormatException {
        if (Commandline.get(this.arguments, "-pymol", false).equals("EXISTS")) {
            this.pymolOutput = true;
            if (!Commandline.get(this.arguments,
                                 "-out",
                                 true).equals("ERROR")
                                ) {
                if (!Commandline.get(this.arguments,
                                     "-out",
                                     true).endsWith(".pml")) {
                    throw new CommandlineArgumentFormatException(
                                   NL
                                 + "ERROR: Please use the file ending \".pml\" "
                                 + "for the output file \""
                                 + Commandline.get(this.arguments, "-out", true)
                                 + "\""
                                 + NL
                                 + NL);
                }
            } else {
                System.err.print(NL
                               + "WARNING: Please make sure that your "
                               + "filename ending is \".pml\", if STDOUT is "
                               + "redirected into a file. "
                               + NL
                               + NL
                                );
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether a PyMOL script should be output.
     * @return {@code TRUE} if PyMOL script is to be output, {@code FALSE}
     * otherwise.
     * @see #readPymolArgument()
     */
    public final boolean isPymolOutputSet() {
        return this.pymolOutput;
    }
    //--------------------------------------------------------------------------

    /**
     * Determines whether the argument -out has been set on the commandline.
     * If the output file already exists, then the user is asked whether
     * the existing file should be overwritten.
     * @return {@code TRUE} if output file should be created,
     * {@code FALSE} otherwise.
     * @see #getOutfileArgument()
     * @see #isOutputFileToBeCreated()
     */
    private boolean readOutfileArgument() {
        if (!Commandline.get(this.arguments, "-out", true).equals("ERROR")) {
            this.outfile = Commandline.get(this.arguments, "-out", true).trim();
            if (!Commandline.get(this.arguments, "-f", false).equals(
                                                                        "EXISTS"
                                                                    )
               ) {
                if (ReadFile.exists(this.outfile)) {
                    System.err.print(NL
                                  + "File \"" + outfile + "\" already exists. "
                                  + "Overwrite? [y/n]: ");
                    String respond = Commandline.get();
                    System.err.print(NL);

                    if (!respond.equalsIgnoreCase("y")
                        &&
                        !respond.equalsIgnoreCase("yes")) {
                        return false;
                    }
                    return true;
                }
                return true;
            }
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the path to the output file.
     * @return String object holding the path to the input file.
     * @see #readOutfileArgument()
     */
    public final String getOutfileArgument() {
        return this.outfile;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns whether an output file should be created.
     * @return {@code TRUE} if output file is to be created, {@code FALSE}
     * otherwise.
     * @see #readOutfileArgument()
     */
    public final boolean isOutputFileToBeCreated() {
        return this.doOutputFile;
    }
    //--------------------------------------------------------------------------

    /**
     * Determines whether the argument -aa1 has been set on the commandline.
     * # characters are placed at both ends of the parameter string.
     * @throws CommandlineArgumentNotFoundException if neither -aa1 nor -r1 has
     * been set on the commandline.
     * @see #getAminoAcidName1Argument()
     * @see #readAminoAcidName2Argument()
     * @see #readAminoAcidNumber1Argument()
     */
    private void readAminoAcidName1Argument() throws
                                          CommandlineArgumentNotFoundException {
        if (Commandline.get(this.arguments, "-aa1", true).equals("ERROR")) {
          if (Commandline.get(this.arguments, "-r1", true).equals("ERROR")) {
           if (Commandline.get(this.arguments, "-dist", true).equals("ERROR")) {
                throw new CommandlineArgumentNotFoundException(NL
                                                             + "ERROR: Could "
                                                             + "NOT find value "
                                                             + "neither for "
                                                             + "parameter \""
                                                             + "-aa1\" nor for "
                                                             + "\"-r1\"."
                                                             + NL
                                                              );
            }
          }
        } else {
            this.residueType1 = "#"
                              + Commandline.get(this.arguments,
                                                "-aa1",
                                                true).trim().toUpperCase()
                              + "#";
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB three letter codes of the 1st type of cross-linked amino
     * acids.
     * @return String object holding the PDB three letter codes.
     * @see #readAminoAcidName1Argument()
     * @see #getAminoAcidName2Argument()
     */
    public final String getAminoAcidName1Argument() {
        return this.residueType1;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -aa2 has been set on the commandline.
     * # characters are placed at both ends of the parameter string.
     * @throws CommandlineArgumentNotFoundException if neither -aa2 nor -r2 has
     *         been set on the commandline.
     * @see #getAminoAcidName2Argument()
     * @see #readAminoAcidName1Argument()
     * @see #readAminoAcidNumber2Argument()
     */
    private void readAminoAcidName2Argument() throws
                                          CommandlineArgumentNotFoundException {
        if (Commandline.get(this.arguments, "-aa2", true).equals("ERROR")) {
          if (Commandline.get(this.arguments, "-r2", true).equals("ERROR")) {
           if (Commandline.get(this.arguments, "-dist", true).equals("ERROR")) {
                throw new CommandlineArgumentNotFoundException(NL
                                                             + "ERROR: Could "
                                                             + "NOT find value "
                                                             + "neither for "
                                                             + "parameter \""
                                                             + "-aa2\" nor for "
                                                             + "\"-r2\"."
                                                             + NL
                                                              );
            }
          }
        } else {
            this.residueType2 = "#"
                              + Commandline.get(this.arguments,
                                                "-aa2",
                                                true).trim().toUpperCase()
                              + "#";
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB three letter codes of the 2nd type of cross-linked amino
     * acids.
     * @return String object holding the PDB three letter codes.
     * @see #readAminoAcidName2Argument()
     * @see #getAminoAcidName1Argument()
     */
    public final String getAminoAcidName2Argument() {
        return this.residueType2;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -r1 has been set on the commandline.
     * # characters are placed at both ends of the parameter string.
     * Beware, that -aa1 parameters will be discarded if they have been set as
     * well.
     * @throws CommandlineArgumentNotFoundException if neither -r1 nor -aa1 has
     *         been set on the commandline.
     * @see #getAminoAcidNumber1Argument()
     * @see #readAminoAcidNumber2Argument()
     * @see #readAminoAcidName1Argument()
     */
    private void readAminoAcidNumber1Argument() throws
                                          CommandlineArgumentNotFoundException {
        if (Commandline.get(this.arguments, "-r1", true).equals("ERROR")) {
          if (Commandline.get(this.arguments, "-aa1", true).equals("ERROR")) {
           if (Commandline.get(this.arguments, "-dist", true).equals("ERROR")) {
                 throw new CommandlineArgumentNotFoundException(NL
                                                              + "ERROR: Could "
                                                              + "NOT find "
                                                              + "value neither "
                                                              + "for parameter "
                                                              + "\"-r1\" nor "
                                                              + "for \"-aa1\"."
                                                              + NL
                                                               );
            }
          }
        } else {
            this.aa1resNo = "#"
                          + Commandline.get(this.arguments,
                                            "-r1",
                                            true).trim().toUpperCase()
                          + "#";
            if (!this.residueType1.equals("")) {
                this.residueType1 = "";
                System.err.print(NL
                               + "WARNING: -aa1 and -r1 are both set. "
                               + "Disregarding -aa1 and considering only \""
                               + "-r1 " + aa1resNo.replaceAll("#", "") + "\""
                               + NL
                               + NL
                                );
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB residue numbers of the 1st type of cross-linked amino
     * acids.
     * @return String object holding the PDB residue numbers.
     * @see #getAminoAcidNumber2Argument()
     * @see #readAminoAcidNumber1Argument()
     */
    public final String getAminoAcidNumber1Argument() {
        return this.aa1resNo;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -r2 has been set on the commandline.
     * # characters are placed at both ends of the parameter string.
     * Beware, that -aa2 parameter will be discarded if they have been set as
     * well.
     * @throws CommandlineArgumentNotFoundException if neither -r2 nor -aa2 has
     *         been set on the commandline.
     * @see #getAminoAcidNumber2Argument()
     * @see #readAminoAcidNumber1Argument()
     * @see #readAminoAcidName2Argument()
     */
    private void readAminoAcidNumber2Argument() throws
                                          CommandlineArgumentNotFoundException {
        if (Commandline.get(this.arguments, "-r2", true).equals("ERROR")) {
          if (Commandline.get(this.arguments, "-aa2", true).equals("ERROR")) {
           if (Commandline.get(this.arguments, "-dist", true).equals("ERROR")) {
                throw new CommandlineArgumentNotFoundException(NL
                                                             + "ERROR: Could "
                                                             + "NOT find value "
                                                             + "neither for "
                                                             + "parameter \""
                                                             + "-r2\" nor for "
                                                             + "\"-aa2\"."
                                                             + NL
                                                              );
            }
          }
        } else {
            this.aa2resNo = "#"
                          + Commandline.get(this.arguments,
                                            "-r2",
                                            true).trim().toUpperCase()
                          + "#";
            if (!this.residueType2.equals("")) {
                this.residueType2 = "";
                System.err.print(NL
                               + "WARNING: -aa2 and -r2 are both set. "
                               + "Disregarding -aa1 and considering only \""
                               + "-r2 " + aa2resNo.replaceAll("#", "") + "\""
                               + NL
                               + NL
                                );
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB residue numbers of the 2nd type of cross-linked amino
     * acids.
     * @return String object holding the PDB residue numbers.
     * @see #getAminoAcidNumber1Argument()
     * @see #readAminoAcidNumber2Argument()
     */
    public final String getAminoAcidNumber2Argument() {
        return this.aa2resNo;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -a1 has been set on the commandline.
     * # characters are placed at both ends of the parameter string.
     * @see #getAtomType1Argument()
     * @see #readAtomType2Argument()
     */
    private void readAtomType1Argument() {
        if (!Commandline.get(this.arguments, "-a1", true).equals("ERROR")) {
            this.atomType1 = "#"
                           + Commandline.get(this.arguments,
                                             "-a1",
                                             true).trim().toUpperCase()
                           + "#";
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB atom name of the 1st type of cross-linked atoms.
     * @return String object holding the PDB atom names.
     * @see #getAtomType2Argument()
     * @see #readAtomType1Argument()
     */
    public final String getAtomType1Argument() {
        return this.atomType1;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -a2 has been set on the commandline.
     * # characters are placed at both ends of the parameter string.
     * @see #getAtomType2Argument()
     * @see #readAtomType1Argument()
     */
    private void readAtomType2Argument() {
        if (!Commandline.get(this.arguments, "-a2", true).equals("ERROR")) {
            this.atomType2 = "#"
                           + Commandline.get(this.arguments,
                                             "-a2",
                                             true).trim().toUpperCase()
                           + "#";
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB atom name of the 2nd type of cross-linked atoms.
     * @return String object holding the PDB atom names.
     * @see #getAtomType1Argument()
     * @see #readAtomType2Argument()
     */
    public final String getAtomType2Argument() {
        return this.atomType2;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -c1 has been set on the commandline.
     * @see #getChainIds1Argument()
     * @see #readChainIds2Argument()
     */
    private void readChainIds1Argument() {
        if (!Commandline.get(this.arguments, "-c1", true).equals("ERROR")) {
            this.chainIds1 = Commandline.get(this.arguments,
                                             "-c1",
                                             true
                                            ).trim().toUpperCase();
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB chain Id of the 1st type of cross-linked amino acids.
     * @return String object holding the PDB chain Id.
     * @see #getChainIds2Argument()
     * @see #readChainIds1Argument()
     */
    public final String getChainIds1Argument() {
        return this.chainIds1;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -c2 has been set on the commandline.
     * @see #getChainIds2Argument()
     * @see #readChainIds1Argument()
     */
    private void readChainIds2Argument() {
        if (!Commandline.get(this.arguments, "-c2", true).equals("ERROR")) {
            this.chainIds2 = Commandline.get(this.arguments,
                                             "-c2",
                                             true).trim().toUpperCase();
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB chain Id of the 2nd type of cross-linked amino acids.
     * @return String object holding the PDB chain Id.
     * @see #getChainIds1Argument()
     * @see #readChainIds2Argument()
     */
    public final String getChainIds2Argument() {
        return this.chainIds2;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -l1 has been set on the commandline.
     * @see #readAlternativeLocation2Argument()
     * @see #getAlternativeLocation1Argument()
     */
    private void readAlternativeLocation1Argument() {
        if (!Commandline.get(this.arguments, "-l1", true).equals("ERROR")) {
            this.altLocs1 = Commandline.get(this.arguments,
                                            "-l1",
                                            true).trim().toUpperCase();
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB alternative location Id of the 1st type of cross-linked
     * amino acids.
     * @return String object holding the PDB alternative location Id.
     * @see #readAlternativeLocation1Argument()
     * @see #getAlternativeLocation2Argument()
     */
    public final String getAlternativeLocation1Argument() {
        return this.altLocs1;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -l2 has been set on the commandline.
     * @see #readAlternativeLocation1Argument()
     * @see #getAlternativeLocation2Argument()
     */
    private void readAlternativeLocation2Argument() {
        if (!Commandline.get(this.arguments, "-l2", true).equals("ERROR")) {
            this.altLocs2 = Commandline.get(this.arguments,
                                            "-l2",
                                            true).trim().toUpperCase();
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the PDB alternative location Id of the 2nd type of cross-linked
     * amino acids.
     * @return String object holding the PDB alternative location Id.
     * @see #readAlternativeLocation2Argument()
     * @see #getAlternativeLocation1Argument()
     */
    public final String getAlternativeLocation2Argument() {
        return this.altLocs2;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -max has been set on the commandline.
     * @see #getMaximumDistanceArgument()
     */
    private void readMaximumDistanceArgument() {
        double max = xwalk.constants.Constants.MAX_SASD_DISTANCE;
        if (!Commandline.get(this.arguments, "-max", true).equals("ERROR")) {
            this.maximumDistance = Double.parseDouble(
                            Commandline.get(this.arguments, "-max", true).trim()
                                                     );

            // There is an upper bound on the SASD distance. Check that this
            // is met when -max is changed by the user.
            if (this.maximumDistance > max
                &&
             !Commandline.get(this.arguments, "-euc", false).equals("EXISTS")) {
               System.err.print(NL + "WARNING: value for -max exceeds " + max
                             + ". Setting -max to " + max + NL);
               this.maximumDistance = max;
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the maximum distance that the cross-linker is able to span.
     * @return double number representing the maximum distance.
     * @see #readMaximumDistanceArgument()
     */
    public final double getMaximumDistanceArgument() {
        return this.maximumDistance;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -euc has been set on the commandline.
     * @see #isSolventPathDistanceCalculationSet()
     */
    private void readSolventPathDistanceArgument() {
        if (Commandline.get(this.arguments, "-euc", false).equals("EXISTS")) {
            this.solventPathDistance = false;
            this.doSAS = false;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether in addition to the Euclidean distance the Solvent-Path
     * distance should also be calculated.
     * @return {@code TRUE} if the Solvent-Path-Distance is to be calculated,
     * {@code FALSE} otherwise.
     * @see #readMaximumDistanceArgument()
     */
    public final boolean isSolventPathDistanceCalculationSet() {
        return this.solventPathDistance;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -radius has been set on the commandline.
     * @see #getSolventRadiusArgument()
     */
    private void readSolventRadiusArgument() {
        if (!Commandline.get(this.arguments, "-radius", true).equals("ERROR")) {
            String arg = Commandline.get(this.arguments, "-radius", true);
            this.solventRadius = Double.parseDouble(arg);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the radius of the solvent, which should be larger than that of
     * water with 1.4 A, because both reactive ends of a cross-linker are much
     * larger than a tri-atomic water molecule.
     * @return double number representing the radius.
     * @see #readSolventRadiusArgument()
     */
    public final double getSolventRadiusArgument() {
        return this.solventRadius;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -space has been set on the commandline.
     * @see #getGridCellSizeArgument()
     */
    private void readGridCellSizeArgument() {
        if (!Commandline.get(this.arguments, "-space", true).equals("ERROR")) {
            this.gridCellLength = Double.parseDouble(
                          Commandline.get(this.arguments, "-space", true).trim()
                                                    );
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the length of the cubic grid cell edge. Grids are used to
     * calculate Solvent-Path-Distances.
     * @return double number representing the length of the cell edge.
     * @see #readGridCellSizeArgument()
     */
    public final double getGridCellSizeArgument() {
        return this.gridCellLength;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -xsas has been set on the commandline.
     * @see #getSolventAccessibiltyArgument()
     */
    private void readSolventAccessibiltyArgument() {
        if (Commandline.get(this.arguments, "-xsas", false).equals("EXISTS")) {
            this.doSAS = false;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the a boolean expression whether SAS is to be calculated or not.
     * @return {@code TRUE} if SAS is to be calculated, {@code FALSE} otherwise.
     * @see #readSolventAccessibiltyArgument()
     */
    public final boolean getSolventAccessibiltyArgument() {
        return this.doSAS;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -intra has been set on the commandline.
     * @see #isIntraMolecularDistanceSet()
     * @see #readInterMolecularDistanceArgument()
     */
    private void readIntraMolecularDistanceArgument() {
        if (Commandline.get(this.arguments, "-intra", false).equals("EXISTS")) {
            this.showIntra = true;
        } else {
            this.showIntra = false;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether only intramolecular distances should be calculated.
     * @return {@code TRUE} if only intramolecular distances are to be
     *         calculated, {@code FALSE} otherwise.
     * @see #readIntraMolecularDistanceArgument()
     * @see #isInterMolecularDistanceSet()
     */
    public final boolean isIntraMolecularDistanceSet() {
        return this.showIntra;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -inter has been set on the commandline.
     * @see #isInterMolecularDistanceSet()
     * @see #readIntraMolecularDistanceArgument()
     */
    private void readInterMolecularDistanceArgument() {
        if (Commandline.get(this.arguments, "-inter", false).equals("EXISTS")) {
            this.showInter = true;
        } else {
            this.showInter = false;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether only intermolecular distances should be calculated.
     * @return {@code TRUE} if only intermolecular distances are to be
     *         calculated, {@code FALSE} otherwise.
     * @see #readInterMolecularDistanceArgument()
     * @see #isIntraMolecularDistanceSet()
     */
    public final boolean isInterMolecularDistanceSet() {
        return this.showInter;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -homo has been set on the commandline.
     * @see #isHomomericSet()
     */
    private void readHomomericArgument() {
        if (Commandline.get(this.arguments, "-homo", false).equals("EXISTS")) {
            this.homo = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether the input file has been set as a homomeric protein
     * complex.
     * @return {@code TRUE} if input file has been set as a homomeric protein
     * complex, {@code FALSE} otherwise.
     * @see #readHomomericArgument()
     */
    public final boolean isHomomericSet() {
        return this.homo;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -global has been set on the commandline.
     * @see #isGlobalGridSet()
     */
    private void readGlobalGridArgument() {
       if (Commandline.get(this.arguments, "-global", false).equals("EXISTS")) {
           this.useGlobalGrid = true;
       }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether a global grid should be used for Solvent-Path distance
     * calculation.
     * @return {@code TRUE} if global grid should be used,
     *         {@code FALSE} otherwise.
     * @see #readGlobalGridArgument()
     */
    public final boolean isGlobalGridSet() {
        return this.useGlobalGrid;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -v has been set on the commandline.
     * @see #isVerboseOutputSet()
     */
    private void readVerboseOutputArgument() {
        if (Commandline.get(this.arguments, "-v", false).equals("EXISTS")) {
            this.verbose = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether additional information about the calculation process
     * should be output.
     * @return {@code TRUE} if additional information should be output,
     *         {@code FALSE} otherwise.
     * @see #readVerboseOutputArgument()
     */
    public final boolean isVerboseOutputSet() {
        return this.verbose;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -grid has been set on the commandline.
     * @see #isGridOutputSet()
     */
    private void readGridOutputArgument() {
        if (Commandline.get(arguments, "-grid", false).equals("EXISTS")) {
            this.verboseGrid = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns whether all grids in the Solvent-Path-Distance calculation should
     * be output.
     * @return {@code TRUE} if grids should be output, {@code FALSE} otherwise.
     * @see #readGridOutputArgument()
     */
    public final boolean isGridOutputSet() {
        return this.verboseGrid;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -trypsin has been set on the
     * commandline.
     * @see #isTrypsinateArgumentSet()
     */
    private void readTrypsinateArgument() {
        if (Commandline.get(this.arguments,
                            "-trypsin",
                            false).equals("EXISTS")) {
            this.doTrypsin = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the a boolean expression whether protein is to be digested or
     * not.
     * @return {@code TRUE} if digestion is to be performed, {@code FALSE}
     * otherwise.
     * @see #readTrypsinateArgument()
     */
    public final boolean isTrypsinateArgumentSet() {
        return this.doTrypsin;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -expasy has been set on the commandline.
     * @see #isExpasyArgumentSet()
     */
    private void readExpasyArgument() {
        if (Commandline.get(this.arguments,
                            "-expasy",
                            false).equals("EXISTS")) {
            this.doExpasy = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the a boolean expression whether protein is to be digested,
     * according to the exceptions listed in <a href="http://expasy.org/tools/
     * peptidecutter/peptidecutter_enzymes.html">ExPASy</a>.
     * @return {@code TRUE} if ExPASy rules should be used, {@code FALSE}
     * otherwise.
     * @see #readExpasyArgument()
     */
    public final boolean isExpasyArgumentSet() {
        return this.doExpasy;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -prob has been set on the commandline.
     * @see #isProbabilityArgumentSet()
     */
    private void readProbabilityArgument() {
        if (Commandline.get(this.arguments,
                            "-prob",
                            false).equals("EXISTS")) {
            this.doProbability = true;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the a boolean expression whether probabilities according to
     * experimental data on DSS and BS3 cross-links should be added to
     * distance output.
     * @return {@code TRUE} if probabilities should be printed out,
     * {@code FALSE} otherwise.
     * @see #readProbabilityArgument()
     */
    public final boolean isProbabilityArgumentSet() {
        return this.doProbability;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines whether the argument -bb has been set on the commandline.
     * If it has been set, than solvent radius is automatically set to
     * xwalk.constants.Constants.SOLVENT_RADIUS_BACKBONE too.
     * @see #isBackboneOnlyArgumentSet()
     */
    private void readBackBoneOnlyArgument() {
        if (Commandline.get(this.arguments, "-bb", false).equals("EXISTS")) {
            this.doBackboneReadOnly = true;
            this.solventRadius =
                              xwalk.constants.Constants.SOLVENT_RADIUS_BACKBONE;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the a boolean expression whether only backbone and beta-carbon
     * coordinates should be read in from the input file. This might be
     * interesting if virtual cross-links are to be created between backbone
     * or beta carbon atoms.
     * @return {@code TRUE} if only backbone and beta-carbon atoms should be
     * read in, {@code FALSE} otherwise.
     * @see #readBackBoneOnlyArgument()
     */
    public final boolean isBackboneOnlyArgumentSet() {
        return this.doBackboneReadOnly;
    }
    //--------------------------------------------------------------------------
}
