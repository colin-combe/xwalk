package external;

import java.io.FileNotFoundException;
import java.io.IOException;

import structure.io.ReadFile;

/**
 * Class making use of functionality provided by Simon Hubbard's NACCESS
 * application.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Naccess {

    /**
     * Path to the NACCESS application.
     */
    private String naccess;

    /**
     * Constructor.
     * @param naccessPath
     *        String holding the path to the NACCESS application.
     * @throws FileNotFoundException
     *         if NACCESS application can not be found at given path.
     */
    public Naccess(final String naccessPath) throws FileNotFoundException {
        if (!ReadFile.exists(naccessPath)) {
            throw new FileNotFoundException("ERROR: NACCESS application does "
                                          + "not reside at " + naccessPath);
        } else {
            this.naccess = naccessPath;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Runs NACCESS and returns total solvent accessible surface area.
     * @param pdbFileName
     *        String holding the path to the PDB file of the protein of which
     *        the total SAS is to be calculated.
     * @return double value representing the total SAS.
     * @throws IOException
     *         if error occurs while reading the PDB file.
     */
    public final double getTotalSurfaceArea(final String pdbFileName)
                                                  throws IOException {
        String fileNameWithoutSuffix = pdbFileName.replaceAll("\\..*", "");
        ExternCommand.execute(this.naccess + " " + pdbFileName, false);
        if (ReadFile.exists(fileNameWithoutSuffix + ".rsa")) {
            ReadFile read = new ReadFile(fileNameWithoutSuffix + ".rsa");
            for (String line : read) {
                if (line.startsWith("TOTAL")) {
                    String[] array = line.split("\\s+");
                    return Double.parseDouble(array[1]);
                }
            }
        } else {
            throw new FileNotFoundException("ERROR: Residue Accessibility file "
                                          + fileNameWithoutSuffix + ".rsa "
                                          + "could not be found");
        }
        return 0;
    }
}
