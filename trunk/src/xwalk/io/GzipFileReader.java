package xwalk.io;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;

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

//return new BufferedInputStream(new GZIPInputStream(
//                                           new FileInputStream(tarGzFilePath))
//                                                   );


}

