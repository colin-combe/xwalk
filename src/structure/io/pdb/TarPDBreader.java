package structure.io.pdb;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.zip.DataFormatException;

import com.ice.tar.TarInputStream;

import structure.io.ReadFile;
import structure.io.TarFileReader;

/**
 * Class handles Tar compressed PDB files.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class TarPDBreader {
    //--------------------------------------------------------------------------
    /**
     * Path to the tar-ed PDB file.
     */
    private String fileName = "";
    //--------------------------------------------------------------------------
    /**
     * InputStream object holding a stream to a tar-ed file.
     */
    private InputStream inputStream;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param path
     *        - String object holding the path to the Gzipped PDB file.
     */
    public TarPDBreader(final String path) {
        this.fileName = path;
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param tarStream
     *        - An InputStream
     */
    public TarPDBreader(final InputStream tarStream) {
        this.inputStream = tarStream;
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
    public final ArrayList < PDBreader > getPDBreaders() throws IOException,
                                                           DataFormatException {

        ArrayList < PDBreader > pdbReaders = new ArrayList < PDBreader > ();

        TarFileReader tarReader = null;
        if (this.inputStream == null) {
            tarReader = new TarFileReader(fileName);
        } else {
            tarReader = new TarFileReader(new TarInputStream(this.inputStream));
        }
        Hashtable < String, ReadFile>  tarContent = tarReader.unTar();
        for (String filePath : tarContent.keySet()) {
            PDBreader pdbReader = new PDBreader(tarContent.get(filePath));
            pdbReader.setFileName(filePath);
            pdbReaders.add(pdbReader);
        }
        return pdbReaders;
    }
    //--------------------------------------------------------------------------
}
