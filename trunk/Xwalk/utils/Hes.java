
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
import java.text.NumberFormat;
import java.util.ArrayList;

import mm.hydrophobicity.Hydrophobicity;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.pdb.PDBreader;
import structure.math.Point3f;
import structure.matter.Atom;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptideList;

/**
 * Class holding for calculating the total HES of a PDB file.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Hes {
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file of the protein complex.
     */
    private String pdbFile;
    /**
     * Do XlogP hydrophobicity calculation.
     */
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
            System.out.println(nL
                             + "java " + Interface.class.getName() + " -help"
                             + nL);
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.println(nL + nL
                           + "java " + Interface.class.getName()
                           + " -in 1b14.pdb"
                           + nL
                           + "This program calculates the total HES of a "
                           + "PDB file." + nL
                           + "Parameters:" + nL
                           + "\t-in <path>\tany structure file in PDB format "
                           + "(required)." + nL
                           + nL
                           );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-in", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-in\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.pdbFile = Commandline.get(args, "-in", true);

            if (!ReadFile.exists(this.pdbFile)) {
                System.err.print(nL
                              + "Couldn't open infile \"" + this.pdbFile
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
    }
    /**
     * Returns REMARK lines with statistics on the non-interface part of a
     * protein surface.
     * @param aminoAcids
     *        List of AminoAcids that form the non-interface surface of a
     *        protein complex.
     * @param hesValues
     *        List of float values for the nonInterfaceSurfaceAminoAcids
     *        object, where each float value represents the HES value for the
     *        same amino acid.
     * @return String object holding PDB-like REMARK lines with statistics on
     *         the non-interface part of the protein surface.
     */
    private String getRemarks(
                       final ArrayList<AminoAcid> aminoAcids,
                       final ArrayList<Float> hesValues) {
        String nL = Constants.LINE_SEPERATOR;
        NumberFormat dec = xwalk.constants.Constants.DISTANCE_DEC_FORMAT;
        StringBuffer output = new StringBuffer();

        int aaCount = 0;
        int atomCount = 0;
        float hesSum = 0;

        int i = 0;
        if (aminoAcids != null) {
            for (AminoAcid aa : aminoAcids) {
                aaCount++;
                for (Atom atom : aa.getAllAtoms()) {
                    atomCount++;
                    if (hesValues != null) {
                        hesSum += hesValues.get(i++);
                    }
                }
            }
        }

        output.append("HEADER   " + new File(this.pdbFile).getName()
                    + nL);
        output.append("REMARK   0  AMINO ACID COUNT: "
                    + aaCount + nL);
        output.append("REMARK   0  ATOM COUNT: "
                    + atomCount + nL);
        output.append("REMARK   0  HES SUM: "
                    + dec.format(hesSum) + nL);
        output.append("REMARK   0  AVERAGE: "
                    + dec.format(hesSum / atomCount)
                    + nL);
        output.append("REMARK" + nL);
        output.append("REMARK   B-factor Column: HES value" + nL);
        output.append("REMARK" + nL);
    return output.toString();
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
        Hes hes = new Hes();
        hes.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(hes.pdbFile);

        } catch (Exception e) {
            System.err.println(e.getMessage());
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);

        //----------------------------------------------------------------------
        // Extract Cartesian coordinates of PDB atoms
        //----------------------------------------------------------------------
        ArrayList<Point3f> complexCoords = new ArrayList<Point3f>();
        for (AminoAcid aa : proteinComplex.getAllAminoAcids()) {
            for (Atom atom : aa.getAllAtoms()) {
                 if (!atom.getElement().getSymbol().equals("H")) {
                     complexCoords.add(atom.getXYZ());
                 }
            }
        }
        //------------------------------------------------------------------
        // Calculate HES
        //------------------------------------------------------------------
        Hydrophobicity hydrophobicity = new Hydrophobicity(proteinComplex);
        ArrayList<Float> hesValues = hydrophobicity.mapHydrophobicity(
                                                                 complexCoords);

        //------------------------------------------------------------------
        // Output interface with HES and conservation grades listed in the
        // occupancy and temperature factor columns.
        //------------------------------------------------------------------

        String output = hes.getRemarks(proteinComplex.getAllAminoAcids(),
                                       hesValues);
        System.out.print(output);
    }
}
