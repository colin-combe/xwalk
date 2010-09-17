package structure.sas;

import java.util.ArrayList;

import structure.constants.Constants;
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

        this.setInterface(protein1, protein2);
    }

    //--------------------------------------------------------------------------
    /**
     * Determines all amino acids at the binding interface of a protein dimer.
     * @param protein1
     *        First protein of two to which interface is to be calculated.
     * @param protein2
     *        Second protein of two to which interface is to be calculated.
     */
    private void setInterface(final PolyPeptide protein1,
                              final PolyPeptide protein2) {

        this.bindingInterface = new ArrayList<ArrayList<AminoAcid>>();

        ArrayList<AminoAcid> aa1 = new ArrayList<AminoAcid>();
        ArrayList<AminoAcid> aa2 = new ArrayList<AminoAcid>();

        for (AminoAcid aminoAcid1 : protein1) {
            for (AminoAcid aminoAcid2 : protein2) {
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
}
