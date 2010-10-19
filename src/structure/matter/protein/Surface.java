package structure.matter.protein;

import java.util.ArrayList;

import structure.io.pdb.PDBwriter;
import external.Naccess;

/**
 * Class for calculating and storing amino acids that are located on the surface
 * of protein molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Surface {
    //--------------------------------------------------------------------------
    /**
     * Minimum relative SASA of an amino acid in order to be considered as
     * a surface residue.
     */
    private static final int MIN_RELATIVE_SOLVENT_ACCESSIBILITY = 5;
    //--------------------------------------------------------------------------
    /**
     * List of amino acids that are forming the surface of a protein.
     */
    private ArrayList<AminoAcid> surfaceAminoAcids = new ArrayList<AminoAcid>();
    //--------------------------------------------------------------------------
    /**
     * Total SASA of the protein molecule.
     */
    private double sasa;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param complex
     *        List of proteins to which the surface will be calculated.
     * @param naccess
     *        The SASA will be calculated using NACCESS.
     */
    public Surface(final PolyPeptideList complex, final Naccess naccess) {
        String tempFileName = "temp.pdb";
        PDBwriter write = new PDBwriter();
        write.setFile(tempFileName);
        write.write(complex);
        naccess.run(tempFileName);

        ArrayList<AminoAcid> allAA = new ArrayList<AminoAcid>();
        for (AminoAcid aa : complex.getAllAminoAcids()) {
            allAA.add(aa.copy());
        }

        naccess.setSolventAccessibility(allAA);
        this.sasa = naccess.getTotalSolventAccessibility();

        for (AminoAcid aa : allAA) {
            if (aa.getRelativeSas()
                >
                Surface.MIN_RELATIVE_SOLVENT_ACCESSIBILITY) {
                this.surfaceAminoAcids.add(aa);
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the amino acids that are located at the surface of a protein
     * molecule.
     * @return List of surface amino acids.
     */
    public final ArrayList<AminoAcid> getSurface() {
        return this.surfaceAminoAcids;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the total SASA of the protein complex.
     * @return double value representing the total SASA.
     */
    public final double getTotalSasa() {
        return this.sasa;
    }
}
