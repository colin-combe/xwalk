package structure.matter.pdb;

import java.io.IOException;
import java.util.ArrayList;

import structure.constants.Constants.ParameterSets;
import structure.matter.Atom;
import structure.matter.AtomList;


/**
 * Class representing protein complexes.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class ProteinComplex extends ArrayList < PolyPeptide > {

    //--------------------------------------------------------------------------
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    //--------------------------------------------------------------------------
    /**
     * Name of complex.
     */
    private String name = "";

    //--------------------------------------------------------------------------
    /**
     * Returns protein atoms found in this complex.
     * @return AtomList object holding all atom coordinates.
     */
    public final AtomList getAllAtoms() {
        AtomList complexCoordinates = new AtomList();
        for (PolyPeptide protein : this) {
             for (AminoAcid aminoacid : protein.getAminoAcids()) {
                  complexCoordinates.addAll(aminoacid.getAllAtoms());
             }
        }
        return complexCoordinates;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the atom radii according to a ParameterSet objectReturns protein
     * atoms found in this complex.
     * @param parameter
     *        - ParameterSet object holding atom radii.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    public final void setAtomRadii(final ParameterSets parameter)
                                                            throws IOException {
        for (PolyPeptide protein : this) {
             for (AminoAcid aminoacid : protein.getAminoAcids()) {
                 for (Atom atom : aminoacid.getAllAtoms()) {
                      atom.setVanDerWaalsRadius(parameter);
                 }
             }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the name of this complex
     * @param name
     *        - String object holding the name of this complex.
     */
    public final void setName(String name) {
        this.name = name;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the name of this complex
     * @return String object holding the name of this complex.
     */
    public final String getName() {
        return name;
    }
}
