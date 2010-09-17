import java.util.ArrayList;

import structure.constants.Constants;
import structure.io.pdb.PDBreader;
import structure.matter.AtomList;
import structure.sas.BindingInterface;

/**
 * Class holding a main method to calculate the binding interfaces for
 * a protein complex.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Interface {
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     */
    protected Interface() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
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
        if (args.length == 0) {
            System.err.print(nL
                          + "Usage: " + nL
                          + nL
                          + "java " + BindingInterface.class.getName()
                          + " protein.pdb " + nL
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

        System.out.print(allInterfacesAtoms);
        
    }
}
