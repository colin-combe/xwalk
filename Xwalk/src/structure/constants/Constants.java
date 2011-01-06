/*
 * (C) 2010 Abdullah Kahraman
 *
 * This software is part of the open-source project "Xwalk". You can use this
 * software under the terms of the
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 * (http://creativecommons.org/licenses/by-nc-sa/3.0/).
 * This means that you
 * 1.) can copy, modify, distribute the software
 * 2.) must give credit to the author
 * 3.) must not use this work for commercial purposes
 * 4.) must license derivative works under the same or a similar license.
 */

package structure.constants;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import structure.matter.Atom;
import structure.matter.AtomList;

/**
 * Various constants that are used by various classes.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
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
     * 0.28 Angstroem. See Laskowski, R.A. (2003) "Structural Quality Assurance",
     * Structural Bioinformatics, 273-303.
     */
    public static final double COORDINATE_UNCERTAINTY = 0.28;
    //--------------------------------------------------------------------------
    /**
     * Returns the uncertainty of the atom's position in Angstroem, calculated
     * from the temperature factor and the average X-ray coordinate uncertainty.
     * @param atom
     *        Atom object from which the uncertainty shall be calculated.
     * @return double value representing the atom's uncertainty
     * @see #getCoordinateUncertainty(AtomList)
     */
    public static double getCoordinateUncertainty(final Atom atom) {
        double errorRange = atom.getAverageDislocation();
        //  + Constants.COORDINATE_UNCERTAINTY
    return errorRange;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the maximum uncertainty of from a list of atom positions in
     * Angstroem, calculated from the temperature factor and the average X-ray
     * coordinate uncertainty.
     * @param atoms
     *        List of Atom object from which the maximum uncertainty will be
     *        calculated.
     * @return double value representing the maximum uncertainty found within
     *         all atoms of the list.
     * @see #getCoordinateUncertainty(Atom)
     */
    public static double getCoordinateUncertainty(final AtomList atoms) {
        double maxUncertainty = Integer.MIN_VALUE;
        for (Atom atom : atoms) {
            maxUncertainty = Math.max(atom.getAverageDislocation()
                                    + Constants.COORDINATE_UNCERTAINTY,
                                      maxUncertainty);
        }
    return maxUncertainty;
    }
    //--------------------------------------------------------------------------
    /**
     * Bond length between a hydrogen and a non-hydrogen.
     */
    public static final double BOND_LENGTH_TO_HYDROGEN = 1.2;
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
     * Minimum difference in the solvent accessible surface area of an atom
     * upon protein complexation in order to be counted at a binding interface.
     */
    public static final double MINUMUM_SASA_DIFFERECE_FOR_INTERFACE = 0.1;
    //--------------------------------------------------------------------------
    /**
     * Minimum relative solvent accessibility of an amino acid in order to be
     * considered to be on the surface of a protein molecule.
     */
    public static final double MINIMUM_REL_SASA_ON_SURFACE = 5;
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
     * Name of temporary directory which is used for example by the
     * BindingInterface class.
     */
    public static final String TEMPORARY_DIR = "Trash";

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
    public enum ParameterSets { RASMOL, SURFNET, MMFF94, PARSE, CHARMM, XLOGP };

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
