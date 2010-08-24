package xwalk.io;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;

import com.ice.tar.TarEntry;
import com.ice.tar.TarInputStream;

/**
 * Class to read tar-ed files into String objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class TarFileReader {

    /**
     * TarInputStream object from which the content of the tar-ed file
     * will be read from.
     */
    private TarInputStream tin;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param tarInputStream
     *        - TarInputStream object holding the content of the tar-ed file.
     */
    public TarFileReader(final TarInputStream tarInputStream) {
        tin = tarInputStream;
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param fileName
     *        - Path to the tarred file.
     * @throws FileNotFoundException if file could not be found.
     */
    public TarFileReader(final String fileName) throws FileNotFoundException {
        tin = new TarInputStream(new FileInputStream(fileName));
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the InputStream of the tar-ed file.
     * @return TarInputStream object holding the content of the tar-ed file.
     */
    public final TarInputStream getTarInputStream() {
        return this.tin;
    }
    //--------------------------------------------------------------------------
    /**
     * Un-tar operation.
     * @return Hashtable with content of the Tar file, where the keys are the
     *         file names and the elements are the file contents.
     * @throws IOException if an Error occurs while reading the TarInputStream
     */
    public final Hashtable < String, String > unTar() throws IOException {
        Hashtable < String, String > entries = new Hashtable < String,
                                                               String > ();
        TarEntry te = this.tin.getNextEntry();
          while (te != null) {
              int buf;
              StringBuffer buffer = new StringBuffer();
              while ((buf = this.tin.read()) != -1) {
                  buffer.append((char) buf);
              }
              entries.put(te.getName(), buffer.toString());
              te = this.tin.getNextEntry();
          }
    return entries;
    }
    //--------------------------------------------------------------------------
}

