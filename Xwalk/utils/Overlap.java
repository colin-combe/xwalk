import java.io.IOException;
import java.util.ArrayList;


import structure.constants.Constants;
import structure.grid.AtomGrid;
import structure.grid.GridCell;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.pdb.PDBreader;
import structure.math.Point3i;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.protein.PolyPeptideList;

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

/**
 * Class holding a main method to assess whether a protein has overlapping
 * volume with a set of other protein. The calculation is based on a grid.
 * a protein complex.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class Overlap {
    /**
     * Empty Constructor.
     */
    protected Overlap() {
    }
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file which will be checked for an overlap.
     */
    private String pdbFile;
    /**
     * Path to a text file holding paths to a set of PDB files, which will
     * be used to generate a grid, which in turn will be used to check whether
     * pdbFile overlaps with the grid.
     */
    private String gridPDBfiles;
    /**
     * Size of grid cell edge length.
     */
    private float gridCellLength = 1.0f;
    /**
     * Threshold for pdbFile to be considered as overlapping.
     */
    private int overlapThreshold = 2;
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
                             + "java " + this.getClass().getName() + " -help"
                             + nL);
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.println(nL + nL
                           + "java " + this.getClass().getName()
                           + " -in 1b14.pdb"
                           + nL
                           + "This program checks whether a PDB file "
                           + "overlaps with a set of other proteins. "
                           + "Parameters:" + nL
                           + "\t-in <path>\tany structure file in PDB format "
                           + "which will be checked for the overlap (required)."
                           + nL
                           + "\t-grid <path>\ta text file holding paths "
                           + "to other PDB files, which will be used to "
                           + "generate a grid, which in turn will be used to "
                           + "check for the overlap (required)."
                           + nL
                           + "\t-cellSize [float]\tsize of grid cells. Default "
                           + "is 1 Angstroem (optional)."
                           + nL
                           + "\t-max [int]\tNumber of times -in PDB file can "
                           + "overlap with -grid PDB files, before it it "
                           + "considered as overlapping. Default is 2 "
                           + "(optional)."
                           + nL + nL
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
                              + "Couldn't open -in file \"" + this.pdbFile
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
        //----------------------------------------------------------------------
        if (Commandline.get(args, "-grid", true).equals("ERROR")) {
            System.err.println(nL + "Error while reading in parameter \"-grid\""
                           + "!!!" + nL);
            System.exit(1);
        } else {

            this.gridPDBfiles = Commandline.get(args, "-grid", true);

            if (!ReadFile.exists(this.gridPDBfiles)) {
                System.err.print(nL
                              + "Couldn't open -grid file \""
                              + this.gridPDBfiles
                              + "\" !!!" + nL + nL);
                System.exit(1);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns an AtomGrid object that holds in the I indices of its grid cell
     * index objects the number of times a protein in the textFileWithPDBpaths
     * file is occupying a grid cell in the I index of the grid cells.
     * @param textFileWithPDBpaths
     *        String object holding the path to the text file that holds a list
     *        of PDB files, which will be used to generate the grid.
     * @return AtomGrid holding the number of times a protein in the
     *         textFileWithPDBpaths file is occupying a grid cell in the I index
     *         of the grid cells.
     */
    private AtomGrid getGrid(final String textFileWithPDBpaths) {
        ArrayList<PolyPeptideList> allProteinComplexes =
                                               new ArrayList<PolyPeptideList>();
        try {
            ReadFile gridPDBfilePath = new ReadFile(textFileWithPDBpaths);
            //------------------------------------------------------------------
            // Read in PDB file
            //------------------------------------------------------------------
            ArrayList<PDBreader> readers = null;
            for (String pdbFilePath : gridPDBfilePath) {
                try {
                    pdbFilePath = pdbFilePath.trim();
                    pdbFilePath.replaceAll("\n", "");
                    readers = PDBreader.createPDBreaders(pdbFilePath);
                } catch (Exception e) {
                    System.err.println("ERROR while reading in PDB file: " + e);
                }

                allProteinComplexes.add(
                               readers.get(0).getEntireProteinComplex().get(0)
                                       );

            }
        } catch (IOException e) {
            System.err.println("ERROR while reading in text file: " + e);
        }

        AtomList allAtoms = new AtomList();
        for (PolyPeptideList proteins : allProteinComplexes) {
            allAtoms.addAll(proteins.getAllAtoms());
        }

        AtomGrid grid = new AtomGrid(allAtoms,
                                     this.gridCellLength,
                                     this.gridCellLength,
                                     false);

        for (PolyPeptideList proteins : allProteinComplexes) {
            grid.reset();
            for (Atom atom : proteins.getAllAtoms()) {
                ArrayList<GridCell> cells = grid.getAllGridCells(atom);
                for (GridCell cell : cells) {
                    if (!cell.isOccupied()) {
                        cell.setIndices(
                               new Point3i(cell.getIndices().getI() + 1, 0, 0));
                        cell.setOccupation();
                    }
                }
            }
        }
        return grid;
    }
    //--------------------------------------------------------------------------
    /**
     * Main method. Reads in a text file, holding a set of PDB files,
     * generates based on these PDB files a grid, and checks whether another
     * PDB file overlaps with the set of PDB files.
     * @param args
     *        Array of Strings holding all commandline arguments.
     */
    public static void main(final String[] args) {
        Overlap overlap = new Overlap();
        overlap.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB file
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(overlap.pdbFile);
        } catch (Exception e) {
            System.err.println(e);
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);


        AtomGrid grid = overlap.getGrid(overlap.gridPDBfiles);
        for (Atom atom : proteinComplex.getAllAtoms()) {
            GridCell closestCell = grid.get(atom);
            if (closestCell != null) {
               if (closestCell.getIndices().getI() > overlap.overlapThreshold) {
                    System.out.println("1: " + overlap.pdbFile
                                     + " is overlapping with PDBs in "
                                     + overlap.gridPDBfiles + ".");
                    System.exit(0);
                }
            }
        }
        System.out.println("0: " + overlap.pdbFile
                         + " is NOT overlapping with PDBs in "
                         + overlap.gridPDBfiles + ".");
    }
}
