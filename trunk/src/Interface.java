import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;

import external.ExternCommand;

import mm.evolution.Consurf;
import mm.hydrophobicity.Hydrophobicity;

import structure.constants.Constants;
import structure.io.Commandline;
import structure.io.ReadFile;
import structure.io.pdb.PDBreader;
import structure.math.Point3d;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.parameter.AminoAcidType;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;
import structure.sas.BindingInterface;

/**
 * Class holding a main method to calculate the binding interfaces for
 * a protein complex.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Interface {
    /**
     * Empty Constructor.
     */
    protected Interface() {
    }
    //--------------------------------------------------------------------------
    /**
     * Path to the PDB formatted file of the protein complex.
     */
    private String pdbFile;
    /**
     * Do XlogP hydrophobicity calculation.
     */
    private boolean doHes = false;
    /**
     * Do evolutionary conservation estimation with ConSurf.
     */
    private boolean doConsurf = false;
    /**
     * Path to the text file with list of ConSurf files to be used in the
     * conservation calculation.
     */
    private String consurfFile;
    /**
     * PDB ID of the protein complex.
     */
    private String pdbId;
    /**
     * Path to the naccess program.
     */
    private String naccessPath;
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
            System.out.print("\njava " + Interface.class.getName() + " -help");
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-help", false).equals("EXISTS")) {
            System.err.print(nL + nL
                           + "java " + Interface.class.getName()
                           + " -in 1b14.pdb"
                           + nL
                           + "This program calculates all interfaces of a "
                           + "protein complex where the interface is defined "
                           + "as all amino acids that are within 9 Angstroem"
                           + "to an atom of another protein chain. In addition "
                           + "a hydrophobic environment score and the "
                           + "evolutionary conservation of the interface amino "
                           + "acids are calculated." + nL
                           + "Parameters:" + nL
                           + "\t-in <path>\tany structure file in PDB format "
                           + "(required)." + nL
                           + "\t-hes [switch]\tAssings the hydrophobic "
                           + "enviroment score for each interface atom "
                           + "to the occupancy column (optional)." + nL
                           + "\t-consurf [switch]\tAssings the evolutionary "
                           + "conservation of each interface amino acid "
                           + "to the temperature factor column (optional)." + nL
                           + "grade files to all protein components of the "
                           + "structure file (optional)." + nL
                           + "\t-consurf <path>\ttext file listing ConSurf "
                           + "grade files to all protein components of the "
                           + "structure file (optional)." + nL
                           + "\t-id [string]\tPDB Id of structure file, if "
                           + "name of file is not its PDB id "
                           + nL
                           + "\t-naccess <string>\tPath to the naccess "
                           + "executable, which will be used to calculate "
                           + "the buried surface area at the interfaces and "
                           + "the average HES and conservation grades for "
                           + "non-interface surface amino acids (optional)."
                           + nL
                           );
            System.exit(0);
        }

        //----------------------------------------------------------------------
        if (Commandline.get(args, "-in", true).equals("ERROR")) {
            System.err.print("\nError while reading in parameter \"-in\""
                           + "!!!\n\n");
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
        //----------------------------
        if (Commandline.get(args, "-hes", false).equals("EXISTS")) {
            doHes = true;
        }
        //----------------------------
        if (Commandline.get(args, "-consurf", false).equals("EXISTS")) {
            doConsurf = true;

            // check whether path has been given by the user
            this.consurfFile = Commandline.get(args, "-consurf", true);
            if (this.consurfFile.startsWith("-")) {
                this.consurfFile = null;
            }
        }
        //----------------------------
        if (!Commandline.get(args, "-id", true).equals("ERROR")) {
            this.pdbId = Commandline.get(args, "-id", true);
        }
        //----------------------------
        if (!Commandline.get(args, "-naccess", true).equals("ERROR")) {
            this.naccessPath = Commandline.get(args, "-naccess", true);
            if (!ExternCommand.exists(this.naccessPath)) {
                System.err.print("ERROR: " + this.naccessPath + " is not "
                               + "executable" + nL);
                System.exit(1);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Determines and sets the conservation grade for each protein in a protein
     * complex using either the ConSurf-DB or local copies of ConSurf
     * conservation grade files.
     * @param proteinComplex
     *        Protein complex to which conservation will be determined.
     * @throws IOException when an error occurs while reading in ConSurf files
     *         from the local drive or from the ConSurf web database.
     */
    private void setConservation(final PolyPeptideList proteinComplex)
                                                            throws IOException {

        String fileName = new File(this.pdbFile).getName().replaceAll(
                                                                     "\\..*", ""
                                                                     );
        String pdbIdent = null;
        if (fileName.matches("[1-9][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z]")
            &&
            fileName.length() == 4) {
            pdbIdent = fileName;
        } else if (this.pdbId != null) {
            pdbIdent = this.pdbId;
        }
        if (pdbIdent != null) {
            for (PolyPeptide protein : proteinComplex) {
                try {
                    Consurf.setConservation(protein, pdbIdent);
                } catch (IOException e) {
                    throw new IOException("ERROR: Could not download ConSurf "
                                        + " grade file for PDB id " + fileName
                                        + ". " + e);
                }
            }
        } else if (this.consurfFile != null) {
            ReadFile read = new ReadFile(this.consurfFile);
            for (PolyPeptide protein : proteinComplex) {
                char proteinChainId = protein.get(0).getAtom(0).getChainId();
                boolean found = false;

                Consurf consurf = null;
                for (String line : read) {
                     consurf = new Consurf(line);
                     if (consurf.getChainId() == proteinChainId) {
                         found = true;
                         break;
                     }
                }
                if (!found) {
                    System.err.println("ERROR: ConSurf file for chain "
                                      + proteinChainId + " is missing.");
                    System.exit(0);
                }
                consurf.assignConservation(protein);
            }
        } else {
            System.err.println("Please specify a ConSurf grade file");
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all BindingInterfaces found between two proteins in the protein
     * complex.
     * @param proteinComplex
     *        Protein complex from which binding sites will be extracted.
     * @return List of BindingInterface object found in the protein complex.
     */
    private ArrayList<BindingInterface> getInterfaces(
                                            final PolyPeptideList proteinComplex
                                                     ) {
        ArrayList<BindingInterface> complexInterfaces =
                                              new ArrayList<BindingInterface>();
        for (int j = 0; j < proteinComplex.size(); j++) {
            for (int k = j + 1; k < proteinComplex.size(); k++) {
                BindingInterface bi = new BindingInterface(
                                                          proteinComplex.get(j),
                                                          proteinComplex.get(k)
                                                          );
                complexInterfaces.add(bi);
            }
        }
        return complexInterfaces;
    }
    //--------------------------------------------------------------------------
    /**
     * Determines all amino acids on the surface of a protein complex using
     * the third party external application NACCESS and a relative ASA of 5%.
     * See Janin, Miller, Chothia (1988), JMB.
     * @param proteinComplex
     *        ProteinComplex object
     * @return List of Aminoacid where each amino acid is found at the surface
     *         as given by its relative ASA.
     */
    private ArrayList<AminoAcid> getSurfaceAminoAcids(
                                            final PolyPeptideList proteinComplex
                                                     ) {
        ArrayList<AminoAcid> surfaceAminoAcids = new ArrayList<AminoAcid>();

        if (this.naccessPath != null) {
            String fileName = new File(this.pdbFile).getName();
            ExternCommand.execute(this.naccessPath + " " + this.pdbFile, false);
            String fileNameBase = fileName.replaceAll("\\..*$", "");

            try {
                ReadFile read = new ReadFile(fileNameBase + ".rsa");
                int i = -1;
                for (String line : read) {
                    if (line.startsWith("RES")) {
                        i++;
                        Atom dummy = new Atom();
                        dummy.setResidueName(line.substring(4,7).trim());
                        dummy.setResidueNumber(Integer.parseInt(
                                                    line.substring(9, 13).trim()
                                                               ));
                        dummy.setChainId(line.substring(8, 9).charAt(0));
                        double relAsa = Double.parseDouble(
                                line.substring(22, 28).trim()
                                       );
                        if (relAsa < 5) {
                            continue;
                        }
                        AminoAcid protAminoAcid =
                                       proteinComplex.getAllAminoAcids().get(i);

                        if (MatterUtilities.equalsResidue(
                                                       dummy,
                                                       protAminoAcid.getAtom(0))
                                                         ) {
                            protAminoAcid.setSas(relAsa);
                        } else {
                            boolean found = false;
                            for (PolyPeptide protein : proteinComplex) {
                                for (AminoAcid aa : protein) {
                                    if (MatterUtilities.equalsResidue(
                                                                  dummy,
                                                                  aa.getAtom(0))
                                                                     ) {
                                        aa.setSas(relAsa);
                                        protAminoAcid = aa;
                                        found = true;
                                        break;
                                    }
                                }
                                if (found) {
                                    break;
                                }
                            }
                            if (!found) {
                                System.err.println("WARNING: Could not find "
                                                 + "ASA amino acid " + line
                                                 + " in the PDB file");
                            }
                        }
                        surfaceAminoAcids.add(protAminoAcid);
                    }
                }
            } catch (IOException e) {
                System.err.println("ERROR: Failed to open " + fileNameBase
                                 + ".rsa : " + e);
                System.exit(1);
            }
        }
        return surfaceAminoAcids;
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
        Interface interfase = new Interface();
        interfase.readCommandline(args);

        //----------------------------------------------------------------------
        // Read in PDB files
        //----------------------------------------------------------------------
        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(interfase.pdbFile);

        } catch (Exception e) {
            System.err.println(e);
        }

        PolyPeptideList proteinComplex =
                                readers.get(0).getEntireProteinComplex().get(0);
        //----------------------------------------------------------------------
        // Calculate interface
        //----------------------------------------------------------------------
        ArrayList<BindingInterface> complexInterfaces =
                                        interfase.getInterfaces(proteinComplex);
        ArrayList<AminoAcid> interfaceAminoAcids =
                                                     new ArrayList<AminoAcid>();
        for (BindingInterface bi : complexInterfaces) {
            interfaceAminoAcids.addAll(bi.getInterface().get(0));
            interfaceAminoAcids.addAll(bi.getInterface().get(1));
        }
        //----------------------------------------------------------------------
        // Calculate non-interface surface amino acids with NACCESS
        //----------------------------------------------------------------------

        ArrayList<AminoAcid> nonInterfaceSurfaceAminoAcids =
                                                     new ArrayList<AminoAcid>();

        if (interfase.naccessPath != null) {
            ArrayList<AminoAcid> surfaceAminoAcids =
                                 interfase.getSurfaceAminoAcids(proteinComplex);

            for (AminoAcid surfaceAA : surfaceAminoAcids) {
                boolean found = false;
                for (AminoAcid interfaceAA : interfaceAminoAcids) {
                    if (interfaceAA == surfaceAA) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    nonInterfaceSurfaceAminoAcids.add(surfaceAA);
                }
            }
        }
        //----------------------------------------------------------------------
        // Extract Cartesian coordinates of interface atoms
        //----------------------------------------------------------------------
        ArrayList<Point3d> interfaceCoords = new ArrayList<Point3d>();
        for (BindingInterface bi : complexInterfaces) {
            for (ArrayList<AminoAcid> halfInterface : bi.getInterface()) {
                for (AminoAcid aa : halfInterface) {
                    for (Atom atom : aa.getAllAtoms()) {
                        interfaceCoords.add(atom.getPoint3d());
                    }
                }
            }
        }
        ArrayList<Double> hesInterfaceCoords = null;
        ArrayList<Double> hesSurfaceCoords = null;
        if (interfase.doHes) {
            //------------------------------------------------------------------
            // Calculate HES for each interface atom
            //------------------------------------------------------------------
            Hydrophobicity hes = new Hydrophobicity(proteinComplex);
            hesInterfaceCoords = hes.mapHydrophobicity(interfaceCoords);
            if (interfase.naccessPath != null) {
                ArrayList<Point3d> nonInterfaceSurfaceCoords =
                                                       new ArrayList<Point3d>();
                for (AminoAcid nonInterfaceSurfaceAA
                                              : nonInterfaceSurfaceAminoAcids) {
                    for (Atom atom : nonInterfaceSurfaceAA.getAllAtoms()) {
                        nonInterfaceSurfaceCoords.add(atom.getPoint3d());
                    }
                }
                hesSurfaceCoords = hes.mapHydrophobicity(
                                                       nonInterfaceSurfaceCoords
                                                        );
            }
        }
        if (interfase.doConsurf) {
            //------------------------------------------------------------------
            // Determine conservation grade for all amino acids in the complex.
            //------------------------------------------------------------------
            try {
                interfase.setConservation(proteinComplex);
            } catch (IOException e) {
                System.err.println("ERROR: Problems while assigning "
                                 + "conservation grades: " + e);
            }
        }

        //------------------------------------------------------------------
        // Output interface with HES and conservation grades listed in the
        // occupancy and temperature factor columns.
        //------------------------------------------------------------------
        int nonInterfaceAaCount = 0;
        int nonInterfaceAtomCount = 0;
        double nonInterfaceSumHes = 0;
        int nonInterfaceSumConservation = 0;
        int nonInterfaceSumTypicInterfaceAA = 0;
        int nonInterfaceSumTypicNonInterfaceAA = 0;

        int i = 0;
        for (AminoAcid aa : nonInterfaceSurfaceAminoAcids) {
            nonInterfaceAaCount++;
            nonInterfaceSumConservation += aa.getConservationGrade();
            if (aa.isInterfaceTypic()) {
                nonInterfaceSumTypicInterfaceAA++;
            }
            if (aa.isNonInterfaceTypic()) {
                nonInterfaceSumTypicNonInterfaceAA++;
            }
            for (Atom atom : aa.getAllAtoms()) {
                nonInterfaceAtomCount++;
                nonInterfaceSumHes += hesSurfaceCoords.get(i++);
            }
        }

        String nL = Constants.LINE_SEPERATOR;
        NumberFormat dec = xwalk.constants.Constants.DISTANCE_DEC_FORMAT;
        StringBuffer output = new StringBuffer();
        output.append("HEADER   " + new File(interfase.pdbFile).getName()
                    + nL);
        output.append("REMARK   0  AMINO ACID COUNT AT NON-INTERFACE: "
                    + nonInterfaceAaCount + nL);
        output.append("REMARK   0  ATOM COUNT AT NON-INTERFACE: "
                    + nonInterfaceAtomCount + nL);
        output.append("REMARK   0  INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: " + nonInterfaceSumTypicInterfaceAA + nL);
        output.append("REMARK   0  NON-INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: "
                    + nonInterfaceSumTypicNonInterfaceAA + nL);
        output.append("REMARK   0  AVERATE INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: "
                    + dec.format((double) nonInterfaceSumTypicInterfaceAA
                                 /
                                 nonInterfaceAaCount) + nL);
        output.append("REMARK   0  AVERAGE NON-INTERFACE TYPIC RESIDUES AT "
                    + "NON-INTERFACE: "
                    + dec.format((double) nonInterfaceSumTypicNonInterfaceAA
                                 /
                                 nonInterfaceAaCount) + nL);
        output.append("REMARK   0  HES SUM AT NON-INTERFACE: "
                    + dec.format(nonInterfaceSumHes) + nL);
        output.append("REMARK   0  CONSERVATION SUM AT NON-INTERFACE: "
                    + nonInterfaceSumConservation + nL);
        output.append("REMARK   0  AVERAGE HES AT NON-INTERFACE: "
                    + dec.format(nonInterfaceSumHes / nonInterfaceAtomCount)
                    + nL);
        output.append("REMARK   0  AVERAGE CONSERVATION AT NON-INTERFACE: "
                    + dec.format((double) nonInterfaceSumConservation
                                 /
                                 nonInterfaceAaCount) + nL);
        output.append("REMARK" + nL);
        output.append("REMARK   Occupancy Column: HES value" + nL);
        output.append("REMARK   B-factor Column: Conservation grade" + nL);
        output.append("REMARK" + nL);

        int j = 0;
        int k = 1;
        for (BindingInterface bi : complexInterfaces) {
            int interfaceAaCount = 0;
            int interfaceAtomCount = 0;
            double interfaceSumHes = 0;
            int interfaceSumConservation = 0;
            int interfaceTypicInterfaceAA = 0;
            int interfaceTypicNonInterfaceAA = 0;
            AtomList atoms = new AtomList();

            for (ArrayList<AminoAcid> halfInterface : bi.getInterface()) {
                for (AminoAcid aa : halfInterface) {
                    interfaceAaCount++;
                    interfaceSumConservation += aa.getConservationGrade();
                    if (aa.isInterfaceTypic()) {
                        interfaceTypicInterfaceAA++;
                    }
                    if (aa.isNonInterfaceTypic()) {
                        interfaceTypicNonInterfaceAA++;
                    }
                    for (Atom atom : aa.getAllAtoms()) {
                        interfaceAtomCount++;
                        if (interfase.doHes) {
                            atom.setOccupancy(hesInterfaceCoords.get(j++));
                            interfaceSumHes += atom.getOccupancy();
                        }
                        if (interfase.doConsurf) {
                            atom.setTemperatureFactor(
                                                       aa.getConservationGrade()
                                                     );
                        }
                        atoms.add(atom);
                    }
                }
            }
            double bsa = -1;
            if (interfase.naccessPath != null) {
                try {
                    bsa = bi.calculateBSA(interfase.naccessPath);
                } catch (IOException e) {
                    System.err.println("WARNING: Failed to successfully run "
                                   + "NACCESS: " + e);
                }
            }
            output.append("REMARK   " + k + "  TITLE " + " "
                       + bi.getInterface().get(0).get(0).getAtom(0).getChainId()
                       + "-"
                       + bi.getInterface().get(1).get(0).getAtom(0).getChainId()
                       + nL);
            output.append("REMARK   " + k + "  AMINO ACID COUNT AT INTERFACE: "
                        + interfaceAaCount + nL);
            output.append("REMARK   " + k + "  ATOM COUNT AT INTERFACE: "
                        + interfaceAtomCount + nL);
            output.append("REMARK   " + k + "  BURIED SURFACE AREA AT "
                        + "INTERFACE: "
                        + dec.format(bsa) + nL);
            output.append("REMARK   " + k + "  INTERFACE TYPIC RESIDUES AT "
                        + "INTERFACE: " + interfaceTypicInterfaceAA + nL);
            output.append("REMARK   " + k + "  NON-INTERFACE TYPIC RESIDUES AT "
                        + "INTERFACE: " + interfaceTypicNonInterfaceAA + nL);
            output.append("REMARK   " + k + "  AVERAGE INTERFACE TYPIC "
                        + "RESIDUES AT INTERFACE: "
                        + dec.format((double) interfaceTypicInterfaceAA
                                     /
                                     interfaceAaCount) + nL);
        output.append("REMARK   " + k + "  AVERAGE NON-INTERFACE TYPIC "
                        + "RESIDUES AT INTERFACE: "
                        + dec.format((double) interfaceTypicNonInterfaceAA
                                     /
                                     interfaceAaCount) + nL);
            output.append("REMARK   " + k + "  HES SUM AT INTERFACE: "
                        + dec.format(interfaceSumHes) + nL);
            output.append("REMARK   " + k + "  CONSERVATION SUM AT INTERFACE: "
                        + interfaceSumConservation + nL);
            output.append("REMARK   " + k + "  AVERAGE HES AT INTERFACE: "
                        + dec.format(interfaceSumHes / interfaceAtomCount)
                        + nL);
            output.append("REMARK   " + k + "  AVERAGE CONSERVATION AT "
                        + "INTERFACE: "
                        + dec.format((double) interfaceSumConservation
                                     /
                                     interfaceAaCount) + nL);
            output.append("REMARK" + nL);
            output.append(atoms.toString() + "END" + nL);

            k++;
        }
        System.out.print(output.toString());
    }
}
