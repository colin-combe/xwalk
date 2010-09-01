package structure.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import xwalk.crosslink.CrossLinkUtilities;

/**
 * Class to read gzipped files.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class GzipFileReader {

    /**
     * Input stream holding the stream to the gzipped file.
     */
    private GZIPInputStream gz;

    /**
     * Constructor.
     * @param fileName
     *        - String object holding the path to the gzipped file.
     * @throws IOException if reading file into GZIPInputStream fails.
     */
    public GzipFileReader(final String fileName) throws IOException {
        this.gz = new GZIPInputStream(new FileInputStream(fileName));
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param gzInputStream
     *        - GZIPInputStream object holding the stream to the gzipped file.
     */
    public GzipFileReader(final GZIPInputStream gzInputStream) {
        this.gz = gzInputStream;
    }
    //--------------------------------------------------------------------------
    /**
     * Return the GZIPInputStream object of this reader.
     * @return GZIPInputStream object holding the stream to the gzipped file.
     */
    public final GZIPInputStream getGZIPInputStream() {
        return this.gz;
    }
    //--------------------------------------------------------------------------
    /**
     * Converts the GzipInputStream object to a BufferedReader object.
     * @return BufferedReader object.
     */
    public final BufferedReader getBufferedReader() {
        InputStreamReader in = new InputStreamReader(this.gz);
        return new BufferedReader(in);

    }
}

