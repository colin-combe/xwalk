package xwalk.matter.pdb;

import java.util.ArrayList;

import xwalk.constants.Constants.ParameterSets;
import xwalk.matter.Atom;
import xwalk.matter.AtomList;

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
     */
    public final void setAtomRadii(final ParameterSets parameter) {
        for (PolyPeptide protein : this) {
             for (AminoAcid aminoacid : protein.getAminoAcids()) {
                 for (Atom atom : aminoacid.getAllAtoms()) {
                      atom.setVanDerWaalsRadius(parameter);
                 }
             }
        }
    }
}
