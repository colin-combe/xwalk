package structure.sas;

import java.io.IOException;
import java.util.ArrayList;

import external.Naccess;

import structure.constants.Constants;
import structure.io.pdb.PDBwriter;
import structure.math.Mathematics;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;


/**
 * Class for handling Interfaces between protein molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class BindingInterface {

    //--------------------------------------------------------------------------
    /**
     * First protein of two to which interface is to be calculated.
     */
    private PolyPeptide proteinA;
    //--------------------------------------------------------------------------
    /**
     * Second protein of two to which interface is to be calculated.
     */
    private PolyPeptide proteinB;
    //--------------------------------------------------------------------------
    /**
     * The binding interface consists of two list of amino acids, where each
     * list comes from one protein partner.
     */
    private ArrayList<ArrayList<AminoAcid>> bindingInterface =
                                         new ArrayList<ArrayList<AminoAcid>>(2);
    //--------------------------------------------------------------------------
    /**
     * Constructor to calculate the amino acids at the binding interface
     * between two protein structures.
     * @param protein1
     *        First protein of two to which interface is to be calculated.
     * @param protein2
     *        Second protein of two to which interface is to be calculated.
     */
    public BindingInterface(final PolyPeptide protein1,
                            final PolyPeptide protein2) {
        this.proteinA = protein1;
        this.proteinB = protein2;
        this.setInterface();
    }

    //--------------------------------------------------------------------------
    /**
     * Determines all amino acids at the binding interface of two proteins.
     */
    private void setInterface() {

        this.bindingInterface = new ArrayList<ArrayList<AminoAcid>>();

        ArrayList<AminoAcid> aa1 = new ArrayList<AminoAcid>();
        ArrayList<AminoAcid> aa2 = new ArrayList<AminoAcid>();

        for (AminoAcid aminoAcid1 : this.proteinA) {
            for (AminoAcid aminoAcid2 : this.proteinB) {
                AtomList atomPair = MatterUtilities.getClosestAtomPair(
                                                       aminoAcid1.getAllAtoms(),
                                                       aminoAcid2.getAllAtoms()
                                                                      );
                if (Mathematics.distance(atomPair.get(0).getPoint3d(),
                                         atomPair.get(1).getPoint3d())
                    <
                    Constants.BINDING_INTERFACE_RADIUS) {
                    if (!aa1.contains(aminoAcid1)) {
                        aa1.add(aminoAcid1);
                    }
                    if (!aa2.contains(aminoAcid2)) {
                        aa2.add(aminoAcid2);
                    }
                }
            }
        }

        this.bindingInterface.add(aa1);
        this.bindingInterface.add(aa2);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the binding interface between protein1 and protein2.
     * @return List of AminoAcid objects that are found at the interface of the
     *         dimer.
     */
    public final ArrayList<ArrayList<AminoAcid>> getInterface() {
        return this.bindingInterface;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the buried surface area at this interface using the
     * Simon Hubbard's NACCESS application.
     * @param naccessPath
     *        String holding the path to the NACCESS application.
     * @return double value representing the buried surface area.
     * @throws IOException
     *         if error occurs while executing the NACCESS application.
     */
    public final double calculateBSA(final String naccessPath)
                                                            throws IOException {
        String proteinAname = "proteinA.pdb";
        String proteinBname = "proteinB.pdb";
        String proteinABname = "proteinAB.pdb";

        PDBwriter writer = new PDBwriter();
        writer.setFile(proteinAname);
        writer.write(this.proteinA);

        Naccess naccess = new Naccess(naccessPath);
        double bsaA = naccess.getTotalSurfaceArea(proteinAname);

        writer.setFile(proteinBname);
        writer.write(this.proteinB);
        double bsaB = naccess.getTotalSurfaceArea(proteinBname);

        writer.setFile(proteinABname);
        writer.write(this.proteinA);
        writer.setFile(proteinABname, true);
        writer.write(this.proteinB);
        double bsaAB = naccess.getTotalSurfaceArea(proteinABname);

        return bsaA + bsaB - bsaAB;
    }
    //--------------------------------------------------------------------------
}
