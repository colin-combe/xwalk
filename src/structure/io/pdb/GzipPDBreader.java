package structure.io.pdb;

import java.io.IOException;
import java.util.zip.DataFormatException;

import structure.io.GzipFileReader;
/**
 * Class handles GNU zipped PDB files.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class GzipPDBreader {
    //--------------------------------------------------------------------------
    /**
     * Path to the gzipped PDB file.
     */
    private String fileName;
    /**
     * Constructor.
     * @param path
     *        - String object holding the path to the Gzipped PDB file.
     */
    //--------------------------------------------------------------------------
    public GzipPDBreader(final String path) {
        this.fileName = path;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a PDBreader object holding all atom coordinates of the PDB file
     * that was gzipped.
     * @return PDBreader object.
     * @throws IOException if an error occurs while reading in the gzipped file.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    public final PDBreader getPDBreader() throws IOException,
                                                 DataFormatException {
        PDBreader reader = new PDBreader(new GzipFileReader(
                                                            this.fileName
                                                           ).getBufferedReader()
                                        );
        reader.setFileName(this.fileName);
    return reader;
    }
    //--------------------------------------------------------------------------
}
