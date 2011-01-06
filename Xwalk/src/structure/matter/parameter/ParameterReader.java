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
 *
 */

package structure.matter.parameter;

import java.io.File;
import java.io.IOException;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.io.ReadFile;


/**
 * Class for reading in and handling parameter files, e.g. atom radius parameter
 * files.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class ParameterReader {

    /**
     * Holds van der Waals radii from the various supported force fields.
     */
    private Hashtable < Element, Double > vdWradii =
                                           new Hashtable < Element, Double >();
    /**
     * Holds atomic XlogP values for all standard PDB amino acid atom types.
     */
    private Hashtable <AminoAcidType, Hashtable <AtomType, Double >> xlogP =
                 new Hashtable <AminoAcidType, Hashtable <AtomType, Double >>();
    /**
     * Location of parameter files
     */
    private static final String PARAMETER_DIR = "mm"
                                              + File.separatorChar
                                              + "parameter"
                                              + File.separatorChar;
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
     * XLOGP radius identifier.
     */
    private static final String XLOGP_FILENAME = "atomic_XlogP.txt";

    /**
     * Column number of element name in two column parameter files. Should be 0.
     */
    private static int elementNameColumn = 0;
    /**
     * Column number of van der Waals radius in two column parameter files.
     * Should be 1.
     */
    private static int vdWradiusColumn = 1;
    /**
     * Column number of amino acid name in three column parameter files.
     * Should be 0.
     */
    private static int aminoAcidNameColumn = 0;
    /**
     * Column number of atom name in three column parameter files. Should be 1.
     */
    private static int atomNameColumn = 1;
    /**
     * Column number of xlogP value in a three column parameter files.
     * Should be 2.
     */
    private static int xlogPcolumn = 2;

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
            this.readParameterSet(ParameterReader.PARAMETER_DIR
                                + ParameterReader.RASMOL_FILENAME);
        }
        if (parameter == Constants.ParameterSets.SURFNET) {
            this.readParameterSet(ParameterReader.PARAMETER_DIR
                                + ParameterReader.SURFNET_FILENAME);
        }
        if (parameter == Constants.ParameterSets.MMFF94) {
            this.readParameterSet(ParameterReader.PARAMETER_DIR
                                + ParameterReader.MMFF94_FILENAME);
        }
        if (parameter == Constants.ParameterSets.PARSE) {
            this.readParameterSet(ParameterReader.PARAMETER_DIR
                                + ParameterReader.PARSE_FILENAME);
        }
        if (parameter == Constants.ParameterSets.CHARMM) {
            this.readParameterSet(ParameterReader.PARAMETER_DIR
                                + ParameterReader.CHARMM_FILENAME);
        }
        if (parameter == Constants.ParameterSets.XLOGP) {
            this.readParameterSet(ParameterReader.PARAMETER_DIR
                                + ParameterReader.XLOGP_FILENAME);
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
                if (column.length == 2) {
                    String elementName =
                                      column[ParameterReader.elementNameColumn];
                    double radius = Double.parseDouble(
                                         column[ParameterReader.vdWradiusColumn]
                                                      );
                    for (Element e : Element.values()) {
                        if (e.toString().equals(elementName)) {
                            this.vdWradii.put(e, radius);
                        }
                    }
                } else if (column.length == 3) {
                    String aminoAcidName =
                                    column[ParameterReader.aminoAcidNameColumn];
                    String atomName = column[ParameterReader.atomNameColumn];
                    double xlogPvalue = Double.parseDouble(
                                             column[ParameterReader.xlogPcolumn]
                                                          );

                    for (AminoAcidType aat : AminoAcidType.values()) {
                        if (aat.getThreeLetterCode().equals(aminoAcidName)) {
                            for (AtomType at : AtomType.values()) {
                                if (at.getAbbreviation().equals(atomName)) {
                                    if (this.xlogP.get(aat) == null) {
                                        Hashtable<AtomType, Double> atomXlogP =
                                              new Hashtable<AtomType, Double>();
                                        atomXlogP.put(at, xlogPvalue);
                                        this.xlogP.put(aat, atomXlogP);
                                    } else {
                                        this.xlogP.get(aat).put(at, xlogPvalue);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Returns the set of atom van der Waals radius parameter.
     * @return Hashtable with Element keys and double elements as radii.
     */
    public final Hashtable < Element, Double > getVdwRadiusParameterSet() {
        return this.vdWradii;
    }
    /**
     * Returns the set of atomic XlogP values for all protein amino acid atoms.
     * @return Hashtable with XlogP values for all protein amino acid atoms.
     */
    public final Hashtable <AminoAcidType, Hashtable <AtomType, Double >>
                                                        getXlogPparameterSet() {
        return this.xlogP;
    }
}
