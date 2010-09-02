package structure.constants;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Various constants that are used by various classes.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Constants {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Constants() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    // CONSTANT NUMERIC VALUES
    //--------------------------------------------------------------------------
    /**
     * The average coordinate uncertainty in PDB crystal structures is about
     * 0.28 Å. See Laskowski, R.A. (2003) "Structural Quality Assurance",
     * Structural Bioinformatics, 273-303.
     */
    public static final double COORDINATE_UNCERTAINTY = 0.28;
    //--------------------------------------------------------------------------
    /**
     * Bond length between a hydrogen and a non-hydrogen.
     */
    public static final double BOND_TO_HYDROGEN = 1.2;
    //--------------------------------------------------------------------------
    /**
     * Default van der Waals radius for any atom type.
     */
    public static final double DEFAULT_ATOM_RADIUS = 1.5;
    //--------------------------------------------------------------------------
    /**
     * Default grid size.
     */
    public static final double DEFAULT_GRID_SIZE = 1.0;
    //--------------------------------------------------------------------------
    /**
     * Default grid size.
     */
    public static final double SOLVENT_RADIUS = 1.4;
    //--------------------------------------------------------------------------
    /**
     * Maximum value that fits into the occupancy and temperature factor
     * column.
     */
    public static final double MIN_OCCUPANCY_TEMPERATURE_VALUE = -99.99;
    //--------------------------------------------------------------------------
    /**
     * Maximum value that fits into the occupancy and temperature factor
     * column.
     */
    public static final double MAX_OCCUPANCY_TEMPERATURE_VALUE = 999.99;
    //--------------------------------------------------------------------------
    /**
     * Minimum XYZ coordinate values in PDB files.
     */
    public static final double MIN_XYZ = -999.999;
    //--------------------------------------------------------------------------
    /**
     * Maximum XYZ coordinate values in PDB files.
     */
    public static final double MAX_XYZ = 9999.999;
    //--------------------------------------------------------------------------
    /**
     * Maximum serial number in PDB files.
     */
    public static final int MAX_SERIAL = 99999;
    //--------------------------------------------------------------------------
    // CONSTANT STRING VALUES
    //--------------------------------------------------------------------------
    /**
     * New line character dependent on the operating system.
     */
    public static final String LINE_SEPERATOR =
                                           System.getProperty("line.separator");
    //--------------------------------------------------------------------------
    /**
     * File path separator character dependent on the operating system.
     */
    public static final String FILE_SEPERATOR = (File.separator.equals("\\")) ?
                                          "\\\\" : File.separator;
    //--------------------------------------------------------------------------
    /**
     * All capital letters in the English alphabet+all decimal numbers+space
     * character. This String object can be useful to select e.g. all chains in
     * a protein molecule.
     */
    public static final String ALPHANUMERIC =
                                         " ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";

    //--------------------------------------------------------------------------
    /**
     * Cartesian coordinates are given up to a three digit after comma.
     */
    public static final NumberFormat CARTESIAN_DEC_FORMAT =
                                                     new DecimalFormat("0.000");
    //--------------------------------------------------------------------------
    /**
     * Occupancy and temperature factor values are given up to a two digit
     * after comma.
     */
    public static final NumberFormat OCCUPANCY_TEMP_DEC_FORMAT =
                                                     new DecimalFormat("0.00");
    //--------------------------------------------------------------------------
    // CONSTANT ENUM SETS
    //--------------------------------------------------------------------------
    /**
     * Supported atom parameter sets.
     */
    public enum ParameterSets { RASMOL, SURFNET, MMFF94, PARSE, CHARMM };

    /**
     * Supported element types.
     */
    public enum ElementTypes { ORGANIC, METAL, METALLOID, NON_METAL };

    /**
     * Supported bond types.
     */
    public enum BondTypes { SINGLE_BOND, DOUBLE_BOND, TRIPLE_BOND,
                            AROMATIC_BOND, CROSS_LINK};
}
