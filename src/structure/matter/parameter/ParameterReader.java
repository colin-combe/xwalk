package structure.matter.parameter;

import java.io.IOException;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.ReadFile;


/**
 * Class for reading in and handling parameter files, e.g. atom radius parameter
 * files.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class ParameterReader {

    /**
     * Holds van der Waals radii from the various supported force fields.
     */
    private Hashtable < Element, Double > vdWradii =
                                           new Hashtable < Element, Double >();
    /**
     * MMRR94 radius identifier.
     */
    private static final String MMFF94_FILENAME = "MMFF94_radii.txt";
    /**
     * PARSE radius identifier.
     */
    private static final String PARSE_FILENAME = "PARSE_radii.txt";
    /**
     * SURFNET radius identifier.
     */
    private static final String SURFNET_FILENAME = "SURFNET_radii.txt";
    /**
     * RASMOL radius identifier.
     */
    private static final String RASMOL_FILENAME = "RASMOL_radii.txt";
    /**
     * CHARMM radius identifier.
     */
    private static final String CHARMM_FILENAME = "CHARMM_radii.txt";

    /**
     * Column number of element name in parameter files. Should be 0.
     */
    private static int elementNameColumn = 0;
    /**
     * Column number of van der Waals radius in parameter files. Should be 1.
     */
    private static int vdWradiusColumn = 1;

    /**
     * default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param parameter
     *        - One of the supported ParameterSet object in Xwalk.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    public ParameterReader(final Constants.ParameterSets parameter)
                                                            throws IOException {
        if (parameter == Constants.ParameterSets.RASMOL) {
            this.readParameterSet(ParameterReader.RASMOL_FILENAME);
        }
        if (parameter == Constants.ParameterSets.SURFNET) {
            this.readParameterSet(ParameterReader.SURFNET_FILENAME);
        }
        if (parameter == Constants.ParameterSets.MMFF94) {
            this.readParameterSet(ParameterReader.MMFF94_FILENAME);
        }
        if (parameter == Constants.ParameterSets.PARSE) {
            this.readParameterSet(ParameterReader.PARSE_FILENAME);
        }
        if (parameter == Constants.ParameterSets.CHARMM) {
            this.readParameterSet(ParameterReader.CHARMM_FILENAME);
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Reads in the information from a parameter file into a Hashtable.
     * @param parameterFileName
     *        - String object holding the path to the desired parameter file.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    private void readParameterSet(final String parameterFileName)
                                                            throws IOException {
        ReadFile read = new ReadFile(parameterFileName);
        for (String line : read) {
            if (!line.startsWith("#")) {
                String[] column = line.split("\t");
                String elementName = column[ParameterReader.elementNameColumn];
                Double radius = Double.parseDouble(
                                         column[ParameterReader.vdWradiusColumn]
                                                  );
                for (Element e : Element.values()) {
                    if (e.toString().equals(elementName)) {
                        this.vdWradii.put(e, radius);
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------

    /**
     * Returns parameter set of atom van der Waals radii.
     * @return Hashtable with Element keys and double elements as radii.
     */
    public final Hashtable < Element, Double > getVdwRadiusParameterSet() {
        return this.vdWradii;
    }
}
