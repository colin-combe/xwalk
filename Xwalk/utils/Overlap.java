import java.io.IOException;
import java.util.ArrayList;


import structure.constants.Constants;
import structure.constants.Constants.ParameterSets;
import structure.grid.AtomGrid;
import structure.grid.GridCell;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.pdb.PDBreader;
import structure.math.Point3i;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
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
    /**
     * Include only these residues for the overlap calculation.
     */
    private String gridPDBzone = "";
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
                           + nL
                           + "\t-zone [string]\tInclude only these "
                           + "residues in the -grid PDB files for the overlap "
                           + "calculation (optional)."
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
        //----------------------------
        if (!Commandline.get(args, "-cellSize", true).equals("ERROR")) {
            this.gridCellLength = Float.parseFloat(
                                        Commandline.get(args, "-cellSize", true)
                                                  );
        }
        //----------------------------
        if (!Commandline.get(args, "-max", true).equals("ERROR")) {
            this.overlapThreshold = Integer.parseInt(
                                             Commandline.get(args, "-max", true)
                                                    );
        }
        //----------------------------
        if (!Commandline.get(args, "-zone", true).equals("ERROR")) {
            this.gridPDBzone =  Commandline.get(args, "-zone", true);
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Reads in all proteins and if selected only a certain region within the
     * proteins and returns them as a list of PolyPeptideList objects.
     * @return A list of PolyPeptideList objects holding all atom coordinates
     * for building the grid.
     */
    private ArrayList<PolyPeptideList> getSelectedAtomCoordinates() {
        ArrayList<PolyPeptideList> allProteinComplexes =
                new ArrayList<PolyPeptideList>();
        try {
            ReadFile gridPDBfilePath = new ReadFile(this.gridPDBfiles);
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

                if (this.gridPDBzone.equals("")) {
                    allProteinComplexes.add(
                            readers.get(0).getEntireProteinComplex().get(0)
                                           );
                } else {
                    String[] zone = this.gridPDBzone.split("-");
                    PolyPeptideList proteinsAll =
                            readers.get(0).getEntireProteinComplex().get(0);
                    PolyPeptideList proteinsSelected = new PolyPeptideList();
                    for (PolyPeptide proteinAll : proteinsAll) {
                        ArrayList<AminoAcid> proteinSelected =
                                                     new ArrayList<AminoAcid>();
                        for (AminoAcid aa : proteinAll) {
                            if (aa.getAtom(0).getResidueNumber()
                                    >=
                                Integer.parseInt(zone[0])
                                &&
                                aa.getAtom(0).getResidueNumber()
                                    <=
                                Integer.parseInt(zone[zone.length - 1])) {
                                proteinSelected.add(aa);
                            }
                        }
                        proteinsSelected.add(new PolyPeptide(proteinSelected));
                    }
                    allProteinComplexes.add(proteinsSelected);
                }

            }
        } catch (IOException e) {
            System.err.println("ERROR while reading in text file: " + e);
        }
        return allProteinComplexes;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns an AtomGrid object that holds in the I indices of its grid cell
     * index objects the number of times a protein in the textFileWithPDBpaths
     * file is occupying a grid cell in the I index of the grid cells.
     * @return AtomGrid holding the number of times a protein in the
     *         textFileWithPDBpaths file is occupying a grid cell in the I index
     *         of the grid cells.
     */
    private AtomGrid calculateGrid() {

        ArrayList<PolyPeptideList> allProteinComplexes =
                                              this.getSelectedAtomCoordinates();


        AtomList allAtoms = new AtomList();
        for (PolyPeptideList proteins : allProteinComplexes) {
            allAtoms.addAll(proteins.getAllAtoms());
        }

        AtomGrid grid = new AtomGrid(allAtoms,
                                     this.gridCellLength,
                                     this.gridCellLength,
                                     false);

        // reset grid's distance values in order to use them for our own
        // purposes.
        for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
            for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                 for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {
                      GridCell cell = grid.get(i, j, k);
                      cell.setDistance(0.0f);
                 }
            }
        }

        for (PolyPeptideList proteins : allProteinComplexes) {
            try {
                proteins.setAtomRadii(ParameterSets.CHARMM);
            } catch (IOException e) {
                System.err.println("ERROR while assigning CHARMM atom radii: "
                                 + e);
            }
            // reset grid indices in order to use them for our own purposes.
            for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
                for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                     for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {
                          GridCell cell = grid.get(i, j, k);
                          cell.unsetOccupation();
                     }
                }
            }
            for (Atom atom : proteins.getAllAtoms()) {
                ArrayList<GridCell> cells = grid.getAllGridCells(atom);
                for (GridCell cell : cells) {
                    if (!cell.isOccupied()) {
                        cell.setDistance(cell.getDistance() + 1);
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


        AtomGrid grid = overlap.calculateGrid();

        for (Atom atom : proteinComplex.getAllAtoms()) {
            GridCell closestCell = grid.get(atom);
            if (closestCell != null) {
               if (Math.round(closestCell.getDistance())
                   >
                   overlap.overlapThreshold) {
                    System.out.println("1: " + overlap.pdbFile
                                     + " is overlapping at least "
                                     + Math.round(closestCell.getDistance()) 
                                     + " times with PDBs in "
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
