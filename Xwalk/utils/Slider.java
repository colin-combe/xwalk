/*
 * (C) 2011 Abdullah Kahraman
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

import java.io.IOException;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.pdb.PDBreader;
import structure.math.Mathematics;
import structure.math.Point3f;
import structure.math.pdb.Transformation;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.CrossLinkUtilities;
import xwalk.io.DistanceReader;

/**
 * Class holding a main method to slide two protein towards each other given
 * a list of distance constraints between both proteins.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Slider {
    /**
     * Empty Constructor.
     */
    protected Slider() {
    }

    //--------------------------------------------------------------------------

    /**
     * Path to the first PDB file.
     */
    private String refFilePath;
    /**
     * Path to the second PDB file.
     */
    private String mobFilePath;
    /**
     * Path to the distance file that holds the distance constraints.
     */
    private String distFilePath;
    /**
     * verbose information to print out.
     */
    private boolean verbose = false;
    /**
     * maximum temperature in simulated annealing.
     */
    private double maxTemperature = 1000;
    /**
     * minimum temperature in simulated annealing.
     */
    private double minTemperature = 1;
    /**
     * maximum MC iteration cycles.
     */
    private int maxIterationCycles = 10000;
    /**
     * last accepted distance sum.
     */
    private double lastAcceptedDistanceSum = Double.MAX_VALUE;
    /**
     * lowest accepted distance sum.
     */
    private double lowestDistanceSum = Double.MAX_VALUE;
    /**
     * conformation with lowest distance sum.
     */
    private PolyPeptideList proteinMobLowest = null;

    //--------------------------------------------------------------------------
    /**
     * Reads all parameter from the commandline.
     * @param args
     *        String array holding all commandline arguments.
     */
    private void readCommandline(final String[] args) {
        String nL = Constants.LINE_SEPERATOR;

        //-----------user information-------------------------------------------
        if (args.length == 0) {
            System.out.print("\njava " + Slider.class.getName() + " -help");
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.print(nL
                           + "\"Slider\" slides two proteins towards "
                           + "each other while minimizing a list of "
                           + "distance constraints."
                           + nL
                           + nL
                           + "Usage:"
                           + nL
                           + "java " + Slider.class.getName()
                           + " -ref 1brsA.pdb -mob 1brsD.pdb -dst 1brs.dist"
                           + nL
                           + nL
                           + "Parameters:"
                           + nL
                           + "\t-ref\t<path>\tFirst protein, which coordinates "
                           + "will be kept fiexed (required)."
                           + nL
                           + "\t-mob\t<path>\tSecond protein, which "
                           + "coordinates will be moved to fullfill the "
                           + " distance constraints (required)."
                           + nL
                           + "\t-dist\t<path>\tDistance file holding at least "
                           + "the first 4 columns of the Xwalk output format. "
                           + "The file will be used to extract the indices and "
                           + "the residue pairs for the distance calculation "
                           + " (required)"
                           + nL
                           + "\t-temp\t[double]\tMaximum temperature used in "
                           + "simulated annealing (default: " + maxTemperature
                           + ", optional)"
                           + nL
                           + "\t-cycles\t[int]\tNumber of Monte Carlo cycles "
                           + "(default: " + maxIterationCycles + ", optional)"
                           + nL
                           + "\t-v\t[switch]\tPrint out verbose information on "
                           + "the STDERR channel (optional)"
                           + nL
                           + nL
                    );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-ref", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-ref\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.refFilePath = Commandline.get(args, "-ref", true);

            if (!ReadFile.exists(this.refFilePath)) {
                System.err.print(nL
                              + "Couldn't open file \"" + this.refFilePath
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-mob", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-mob\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.mobFilePath = Commandline.get(args, "-mob", true);

            if (!ReadFile.exists(this.mobFilePath)) {
                System.err.print(nL
                              + "Couldn't open file \"" + this.mobFilePath
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-dist", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-dist\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.distFilePath = Commandline.get(args, "-dist", true);

            if (!ReadFile.exists(this.distFilePath)) {
                System.err.print(nL
                              + "Couldn't open file \"" + this.distFilePath
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-v", false).equals("EXISTS")) {
            this.verbose = true;
        }
        //----------------------------------------------------------------------
        if (!Commandline.get(args, "-temp", true).equals("ERROR")) {
            this.maxTemperature = Double.parseDouble(
                                     Commandline.get(args, "-temp", true)
                                                 );
        }
        //----------------------------------------------------------------------
        if (!Commandline.get(args, "-cycles", true).equals("ERROR")) {
            this.maxIterationCycles = Integer.parseInt(
                                     Commandline.get(args, "-cycles", true)
                                                      );
        }
    }
    //--------------------------------------------------------------------------
    private static Point3f getRandomTranslationVector() {
        Point3f translationVector = new Point3f(
                (float) (Math.random() * 4 - 2),
                (float) (Math.random() * 4 - 2),
                (float) (Math.random() * 4 - 2)
                   );
        return translationVector;
    }
    //--------------------------------------------------------------------------
    private static double[][] getRandomRatationMatrix(){
        double[][] rotationMatrix = Mathematics.getEulerRotationMatrix(
                Math.random() * 4 * Math.PI - 2 * Math.PI,
                Math.random() * 4 * Math.PI - 2 * Math.PI,
                Math.random() * 4 * Math.PI - 2 * Math.PI);
        return rotationMatrix;
    }
    
    //--------------------------------------------------------------------------
    private static double getDistSum(final PolyPeptideList proteinRef,
                              final PolyPeptideList proteinMobCopy,
                              final CrossLinkList constraintsList) {
        // calculate distance sum
        PolyPeptideList complexCopy = new PolyPeptideList();
        complexCopy.addAll(proteinRef);
        complexCopy.addAll(proteinMobCopy);
        Hashtable < Atom, AtomList > relevantAtomPairs = null;
        try {
            relevantAtomPairs = CrossLinkUtilities.extractRelevantPairs(
                                                             complexCopy,
                                                             constraintsList
                                                                       );
        } catch (IOException e) {
            System.err.println("ERROR while attempting to extract "
                             + "cross-linked residues from the input PDB "
                             + "files: " + e);
            System.exit(1);
        }

        double distSum = 0;
        int n = 0;
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            for (Atom atom2 : relevantAtomPairs.get(atom1)) {
                distSum += Mathematics.distance(atom1.getXYZ(),
                                                atom2.getXYZ());
                n++;
            }
        }
        distSum /= n;

    return distSum;
    }

    //--------------------------------------------------------------------------

    private static double getBoltzmannProbability(
                                           final double lastAcceptedDistanceSum,
                                           final double distSum,
                                           final double temperature) {
        // check whether distance is smaller, in which case do the move
        // on the original mobile protein
        double boltzmannFactor = (lastAcceptedDistanceSum - distSum)
                                 /
                                 temperature;
        double probability = Math.exp(boltzmannFactor);
    return probability;
    }

    //--------------------------------------------------------------------------
    private static boolean doTransformation(
                                           final double lastAcceptedDistanceSum,
                                           final double distSum,
                                           final double temperature,
                                           final boolean verbose) {

        double probability = Slider.getBoltzmannProbability(
                                                        lastAcceptedDistanceSum,
                                                        distSum,
                                                        temperature);

        boolean doMove = false;
        if (probability >= 1) {
            doMove = true;
            if (verbose) {
                System.err.println("ACCEPT: Better score: "
                          + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
            }
        } else if (probability >= Math.random()) {
            doMove = true;
            if (verbose) {
                System.err.println("ACCEPT: By thermal probability: "
                        + Constants.CARTESIAN_DEC_FORMAT.format(probability)
                        + "\t"
                        + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
            }
        } else {
            doMove = false;
            if (verbose) {
                System.err.println("FAILED: "
                        + Constants.CARTESIAN_DEC_FORMAT.format(probability)
                        + "\t"
                        + Constants.CARTESIAN_DEC_FORMAT.format(distSum));
            }
        }
        return doMove;
    }

    //--------------------------------------------------------------------------
    private void doMC(final PolyPeptideList proteinRef,
                      final PolyPeptideList proteinMob,
                      final CrossLinkList constraintsList,
                      int maxIterationCycles,
                      final double temperature) {
        // continue move attempts for a maximum number of move attempts or
        // until move resulted in a sufficient conformation
        while (maxIterationCycles-- > 0) {
            //----------------------RANDOM TRANSLATION--------------------------
            // generate a random translation vector with -1 to 1 coordinate
            // values.
            Point3f translationVector = Slider.getRandomTranslationVector();

            // create copy of mobile protein to test move
            PolyPeptideList proteinMobCopy = new PolyPeptideList();
            for (PolyPeptide protein : proteinMob) {
                proteinMobCopy.add(protein.copy());
            }

            // do the random move
            Transformation.move(proteinMobCopy.getAllAtoms(),
                                translationVector);


            double distSum = Slider.getDistSum(proteinRef,
                                               proteinMobCopy,
                                               constraintsList);


            boolean doMove = Slider.doTransformation(
                                                   this.lastAcceptedDistanceSum,
                                                   distSum,
                                                   temperature,
                                                   this.verbose);

            if (doMove) {
                Transformation.move(proteinMob.getAllAtoms(),
                                    translationVector);
                this.lastAcceptedDistanceSum = distSum;
            }
            if (distSum < this.lowestDistanceSum) {
                this.proteinMobLowest = proteinMobCopy;
                this.lowestDistanceSum = distSum;
            }

            //----------------------RANDOM ROTATION-----------------------------
            // do the random move
            // do rotation with Euler angles for which we first need to
            // translate the protein to the coordinate center.
            double[][] rotationMatrix = Slider.getRandomRatationMatrix();

            // create copy of mobile protein to test move
            proteinMobCopy = null;
            proteinMobCopy = new PolyPeptideList();
            for (PolyPeptide protein : proteinMob) {
                proteinMobCopy.add(protein.copy());
            }

            Transformation.rotateAtOrigin(proteinMobCopy.getAllAtoms(),
                                          rotationMatrix);

            distSum = Slider.getDistSum(proteinRef,
                                        proteinMobCopy,
                                        constraintsList);

            boolean doRotation = Slider.doTransformation(
                                                   this.lastAcceptedDistanceSum,
                                                   distSum,
                                                   temperature,
                                                   this.verbose);
            if (doRotation) {
                Transformation.rotateAtOrigin(proteinMob.getAllAtoms(),
                                              rotationMatrix);
                this.lastAcceptedDistanceSum = distSum;
            }

            if (distSum < this.lowestDistanceSum) {
                this.proteinMobLowest = proteinMobCopy;
                this.lowestDistanceSum = distSum;
            }

        }
    }
    //--------------------------------------------------------------------------
    /**
     * Main method. Reads in a PDB file, whose path is given as a first argument
     * on the commandline and calculates binding interfaces to all protein
     * chains given in the PDB file.
     * @param args
     *        Array of Strings holding all commandline arguments.
     */
    public static void main(final String[] args) {
        Slider slider = new Slider();
        slider.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        PDBreader readersRef = null;
        PDBreader readersMob = null;
        CrossLinkList constraintsList = null;
        try {
            readersRef = PDBreader.createPDBreaders(slider.refFilePath).get(0);
            readersMob = PDBreader.createPDBreaders(slider.mobFilePath).get(0);
            constraintsList = DistanceReader.getCrossLinks(slider.distFilePath,
                                                           false,
                                                           false,
                                                           true);
        } catch (Exception e) {
            System.err.println("ERROR while reading in input files: " + e);
            System.exit(1);
        }

        PolyPeptideList proteinRef =
                                    readersRef.getEntireProteinComplex().get(0);
        PolyPeptideList proteinMob =
                                    readersMob.getEntireProteinComplex().get(0);

        double temperature = slider.maxTemperature;
        double temperatureDiff = slider.maxIterationCycles
                                 /
                                 slider.maxTemperature;
        int numberOfCycles = Math.round((float)
                                        (slider.maxIterationCycles
                                        /
                                        (slider.maxTemperature
                                         /
                                         temperatureDiff)));

        while (temperature >= 0) {
            System.out.println("TEMPERATURE: " + temperature);
            System.out.println("CYCLES: " + (Math.round((float) numberOfCycles)));

            slider.doMC(proteinRef, proteinMob, constraintsList,
                        numberOfCycles, temperature);
            temperature -= temperatureDiff;
        }

        if (slider.verbose) {
            System.err.println("FINAL: "
                    + Constants.CARTESIAN_DEC_FORMAT.format(
                                                     slider.lowestDistanceSum));
        }
        System.out.print(slider.proteinMobLowest);
    }
}
