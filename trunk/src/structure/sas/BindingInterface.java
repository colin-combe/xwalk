package structure.sas;

import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.pdb.PDBreader;
import structure.math.Mathematics;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;


/**
 * Class for handling Interfaces between protein molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class BindingInterface {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected BindingInterface() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all amino acids at the binding interface of a protein dimer.
     * @param protein1
     *        First protein of two to which interface is to be calculated.
     * @param protein2
     *        Second protein of two to which interface is to be calculated.
     * @return List of AminoAcid objects that are found at the interface of the
     *         dimer.
     */
    public static ArrayList<ArrayList<AminoAcid>> getInterface(
                                                     final PolyPeptide protein1,
                                                     final PolyPeptide protein2
                                                              ) {

        ArrayList<ArrayList<AminoAcid>> bindingInterface =
                                        new ArrayList<ArrayList<AminoAcid>>();

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

        bindingInterface.add(aa1);
        bindingInterface.add(aa2);

    return bindingInterface;
    }
    //--------------------------------------------------------------------------
    /**
     * Main method. Reads in a PDB file, whose path is give as a first argument
     * on the commandline and calculates binding interfaces to all protein
     * chains given in the PDB file.
     * @param args
     *        Array of Strings holding all commandline arguments.
     */
    public static void main(final String[] args) {
        String nL = Constants.LINE_SEPERATOR;
        if (args.length == 0) {
            System.err.print(nL + "Usage: " + nL + nL + "java "
                           + BindingInterface.class.getName()
                           + " protein.pdb " + nL + "or " + nL + "java "
                           + BindingInterface.class.getName()
                           + " decoys.tar.gz map.pdb " + nL + nL);
            System.exit(1);
        }

        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(args[0]);

        } catch (Exception e) {
            System.err.print(e + nL);
        }

        AtomList allInterfaces = new AtomList();
        String allInterfacesIds = "";
        Hashtable<String, Integer> idsCount = new Hashtable<String, Integer>();
        int complexCount = 0;
        for (PDBreader files : readers) {
            ArrayList<PolyPeptideList> fileMolecules =
                                                files.getEntireProteinComplex();
            for (int i = 0; i < fileMolecules.size(); i++) {
                PolyPeptideList proteinComplex = fileMolecules.get(i);

                AtomList complexInterfaceAtoms = new AtomList();
                String complexInterfaceAtomsIds = "";
                complexCount++;

                for (int j = 0; j < proteinComplex.size(); j++) {
                    for (int k = j + 1; k < proteinComplex.size(); k++) {
                        ArrayList<ArrayList<AminoAcid>> complexInterface =
                                               BindingInterface.getInterface(
                                                          proteinComplex.get(j),
                                                          proteinComplex.get(k)
                                                                            );
                        for (ArrayList<AminoAcid> interfaceHalf
                                                           : complexInterface) {
                            for (AminoAcid aa : interfaceHalf) {
                                for (Atom atom : aa.getAllAtoms()) {
                                    String id = "#"
                                              + atom.getName()
                                              + atom.getResidueNumber()
                                              + atom.getResidueName()
                                              + "-"
                                              + atom.getChainId()
                                              + "#";
                                    if (allInterfacesIds.indexOf(id) == -1
                                        &&
                                        complexInterfaceAtomsIds.indexOf(id)
                                                                        == -1) {
                                        idsCount.put(id, 1);
                                        atom.setTemperatureFactor(
                                                                idsCount.get(id)
                                                                 );

                                        complexInterfaceAtoms.add(atom);
                                        complexInterfaceAtomsIds += id;
                                    } else if (complexInterfaceAtomsIds.indexOf(
                                                                              id
                                                                               )
                                                                        == -1) {
                                        idsCount.put(id, idsCount.get(id) + 1);
                                        atom.setTemperatureFactor(
                                                                idsCount.get(id)
                                                                 );

                                        complexInterfaceAtoms.add(atom);
                                        complexInterfaceAtomsIds += id;
                                    }
                                }
                            }
                        }
                    }
                }
                allInterfaces.addAll(complexInterfaceAtoms);
                allInterfacesIds += complexInterfaceAtomsIds;
            }
        }
        if (args.length == 1) {
            System.out.print(allInterfaces);
        }
        if (args.length == 2) {
            PDBreader map = null;
            try {
                map = new PDBreader(args[1]);
            } catch (Exception e) {
                System.err.print(e + nL);
            }
            // init occupancy and temperature factor values to 0.0.
            ArrayList<PolyPeptideList> complexes =
                                                  map.getEntireProteinComplex();
            for (PolyPeptideList complex : complexes) {
                for (PolyPeptide protein : complex) {
                    for (AminoAcid aa : protein) {
                        for (Atom atom : aa.getAllAtoms()) {
                            atom.setOccupancy(0.0);
                            atom.setTemperatureFactor(0.0);
                        }
                    }
                }
            }

            for (PolyPeptideList complex : complexes) {
                for (PolyPeptide protein : complex) {
                    for (AminoAcid aa : protein) {
                        for (Atom atom1 : aa.getAllAtoms()) {

                            for (Atom atom2 : allInterfaces) {
                                if (MatterUtilities.equalsResidue(atom1, atom2)
                                   &&
                                   atom1.getName().equals(atom2.getName())) {
                                    atom1.setOccupancy(
                                                    atom2.getTemperatureFactor()
                                                      );
                                    atom1.setTemperatureFactor(
                                                    atom2.getTemperatureFactor()
                                                               /
                                                    complexCount
                                                              );
                                }
                            }
                        }
                    }
                }
                System.out.print(complex + "END" + nL);
            }
        }
    }
}
