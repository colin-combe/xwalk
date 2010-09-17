import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.pdb.PDBreader;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;
import structure.sas.BindingInterface;

/**
 * Class holding a main method to calculate the average binding interface for
 * a protein complex have different conformational topologies.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class AvgInterface {

    /**
     * List of Atom objects forming all binding site atoms.
     */
    private AtomList allInterfacesAtoms = new AtomList();

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param readers
     *        List of PDBreader objects holding each the content of a single PDB
     *        file. The PDB files might be compressed into a tar or gzip file.
     */
    public AvgInterface(final ArrayList<PDBreader> readers) {
        this.setUnionOfBindingInterfaces(readers);
    }
    //--------------------------------------------------------------------------
    /**
     * Determines all binding interfaces found in a list of PDB files.
     * @param readers
     *        List of PDBreader objects holding each the content of a single PDB
     *        file. The PDB files might be compressed into a tar or gzip file.
     */
    private void setUnionOfBindingInterfaces(
                                             final ArrayList<PDBreader> readers
                                            ) {

        String allInterfacesIds = "";
        Hashtable<String, Integer> idsCount = new Hashtable<String, Integer>();


        for (PDBreader files : readers) {
            ArrayList<PolyPeptideList> fileMolecules =
                                                files.getEntireProteinComplex();
            for (int i = 0; i < fileMolecules.size(); i++) {
                PolyPeptideList proteinComplex = fileMolecules.get(i);

                AtomList complexInterfaceAtoms = new AtomList();
                String complexInterfaceAtomsIds = "";

                for (int j = 0; j < proteinComplex.size(); j++) {
                    for (int k = j + 1; k < proteinComplex.size(); k++) {
                        BindingInterface bi = new BindingInterface(
                                                          proteinComplex.get(j),
                                                          proteinComplex.get(k)
                                                                  );
                        ArrayList<ArrayList<AminoAcid>> complexInterface =
                                                              bi.getInterface();

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
                allInterfacesIds += complexInterfaceAtomsIds;
                this.allInterfacesAtoms.addAll(complexInterfaceAtoms);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a list of all binding site atoms found in the union of the PDB
     * files.
     * @return List of AminoAcid objects that are found at the interface of the
     *         dimer.
     */
    public final AtomList getInterfacesAtoms() {
        return this.allInterfacesAtoms;
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
        String nL = Constants.LINE_SEPERATOR;
        if (args.length < 2) {
            System.err.print(nL
                          + "Usage: " + nL
                          + nL
                          +  "java " + BindingInterface.class.getName()
                          + " decoys.tar.gz map.pdb " + nL
                          + nL);
            System.exit(1);
        }

        ArrayList<PDBreader> readers = null;
        try {
            readers = PDBreader.createPDBreaders(args[0]);

        } catch (Exception e) {
            System.err.print(e + nL);
        }

        int complexCount = 0;
        for (PDBreader files : readers) {
            for (int i = 0; i < files.getEntireProteinComplex().size(); i++) {
                complexCount++;
             }
        }

        AvgInterface avgInterface = new AvgInterface(readers);
        AtomList allInterfacesAtoms = avgInterface.getInterfacesAtoms();

        PDBreader map = null;
        try {
            map = new PDBreader(args[1]);
        } catch (Exception e) {
            System.err.print(e + nL);
        }
        // init occupancy and temperature factor values to 0.0.
        ArrayList<PolyPeptideList> complexes = map.getEntireProteinComplex();
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

                        for (Atom atom2 : allInterfacesAtoms) {
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
