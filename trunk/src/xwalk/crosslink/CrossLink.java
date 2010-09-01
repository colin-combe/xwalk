package xwalk.crosslink;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import structure.constants.Constants;
import structure.constants.Constants.BondTypes;
import structure.constants.Constants.Value;
import structure.grid.Path;
import structure.math.Mathematics;
import structure.matter.Atom;
import structure.matter.Bond;
import structure.matter.pdb.AminoAcid;


/**
 * Class for representing Cross-Link objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class CrossLink extends Bond {

    //--------------------------------------------------------------------------
    /**
     * Distance that the cross-link spans in sequence space.
     */
    private int seqDist;
    /**
     * Distance that the cross-link spans in Euclidean space.
     */
    private double eucDist;
    /**
     * Distance that the cross-link spans in Solvent-Path distance space.
     */
    private double solventPathDistance =
                                Double.parseDouble(Value.DISTANCE.getDefault());

    /**
     * First protein atom that is connected by the cross-link.
     */
    private Atom preAtom;
    /**
     * Second protein atom that is connected by the cross-link.
     */
    private Atom postAtom;
    /**
     * Second protein atom that is connected by the cross-link.
     */
    private Path solventDistancePath;
    /**
     * Ranking index of this cross-link within a list of cross-links.
     */
    private int index = -1;
    /**
     * String object holding the path to the PDB file in which the cross-link
     * has been found.
     */
    private String filePath = "";
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atom2
     *        - Second protein atom to be connected by the virtual cross-linker.
     */
    public CrossLink(final Atom atom1, final Atom atom2) {
        super(atom1, atom2, BondTypes.CROSS_LINK);
        ArrayList < Atom > list = new ArrayList < Atom > ();
        list.add(atom1);
        list.add(atom2);
        // sorting atom pair by chain id.
        Collections.sort(list, new Comparator < Atom > () {
                                  public int compare(final Atom atom1,
                                                     final Atom atom2) {
                                     String chainId1 = atom1.getChainId() + "";
                                     String chainId2 = atom2.getChainId() + "";
                                     return chainId1.compareTo(chainId2);
                                  }
                              }
                        );
        this.preAtom = list.get(0);
        this.postAtom = list.get(1);

        this.setSequenceDistance();
        this.setEuclideanDistance();
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor. No Euclidean distance calculations are performed.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atom2
     *        - Second protein atom to be connected by the virtual cross-linker.
     * @param euclideanDistance
     *        - double value representing the Euclidean distance between both 
     *          atoms.
     * @param sequenceDistance
     *        - integer value representing the distance in sequence space of  
     *          both atoms.
     */
    public CrossLink(
                     final Atom atom1,
                     final Atom atom2,
                     final int sequenceDistance,
                     final double euclideanDistance
                    ) {
        super(atom1, atom2, BondTypes.CROSS_LINK);
        ArrayList < Atom > list = new ArrayList < Atom > ();
        list.add(atom1);
        list.add(atom2);
        // sorting atom pair by chain id.
        Collections.sort(list, new Comparator < Atom > () {
                                  public int compare(final Atom atom1,
                                                     final Atom atom2) {
                                     String chainId1 = atom1.getChainId() + "";
                                     String chainId2 = atom2.getChainId() + "";
                                     return chainId1.compareTo(chainId2);
                                  }
                              }
                        );
        this.preAtom = list.get(0);
        this.postAtom = list.get(1);

        this.seqDist = sequenceDistance; 
        this.eucDist = euclideanDistance;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the distance in sequence space, which corresponds to the
     * difference of the amino acid rank to which both atoms are associated.
     */
    private void setSequenceDistance() {
        this.seqDist = Math.abs(
                                         postAtom.getRank() - preAtom.getRank()
                                        );
    }
    //--------------------------------------------------------------------------

    /**
     * Calculates the distance in Euclidean space.
     */
    private void setEuclideanDistance() {
        this.eucDist = Mathematics.distance(preAtom.getPoint3d(),
                                                      postAtom.getPoint3d());
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the Solvent-Path distance. In order to calculate the Solvent-Path
     * distance please see in the Class CrossLinkList.
     * @param dist
     *        - Distance in the Solvent-Path space.
     */
    public final void setSolventPathDistance(final double dist) {
        this.solventPathDistance = dist;
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the grid path of the Solvent-Path distance.
     * @param path -
     *        Path object holding the list of GridCell object between source
     *        and target grid cells within a Grid object.
     */
    public final void setPath(final Path path) {
        this.solventDistancePath = path;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the distance in sequence space.
     * @return integer number representing the distance in sequence space.
     */
    public final int getSequenceDistance() {
        return this.seqDist;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the distance in Euclidean space.
     * @return double number representing the distance in Euclidean space.
     */
    public final double getEuclideanDistance() {
        return Double.parseDouble(Constants.DISTANCE_DEC_FORMAT.format(
                                                                    this.eucDist
                                                                      ));
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the distance in Solvent-Path distance space.
     * @return double number representing the distance in Solvent-Path distance
     *         space.
     */
    public final double getSolventPathDistance() {
        return Double.parseDouble(Constants.DISTANCE_DEC_FORMAT.format(
                                                        this.solventPathDistance
                                                                      ));
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the grid path of the Solvent-Path distance.
     * @return Path object holding the list of GridCell object between source
     *         and target grid cells within a Grid object.
     */
    public final Path getPath() {
        return this.solventDistancePath;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether a second cross-link has the same atom identifier, i.e.
     * residue name and residue number.
     * @param crossLink
     *        - CrossLink object to be compared to this CrossLink object.
     * @return {@code TRUE} if both CrossLink object are equal in homology,
     * {@code FALSE} otherwise.
     */
    public final boolean equalsInHomolog(final CrossLink crossLink) {
        Atom preAtom1 = this.getPreAtom();
        Atom postAtom1 = this.getPostAtom();
        Atom preAtom2 = crossLink.getPreAtom();
        Atom postAtom2 = crossLink.getPostAtom();

        String residueId1 = preAtom1.getResidueName() + "#"
                            + preAtom1.getResidueNumber();
        String residueId2 = postAtom1.getResidueName() + "#"
                            + postAtom1.getResidueNumber();
        String residueId3 = preAtom2.getResidueName() + "#"
                            + preAtom2.getResidueNumber();
        String residueId4 = postAtom2.getResidueName() + "#"
                            + postAtom2.getResidueNumber();
        if ((residueId1.equals(residueId3) && residueId2.equals(residueId4))
            ||
            (residueId2.equals(residueId3) && residueId1.equals(residueId4))) {
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether a second cross-link has the same atom identifier, i.e.
     * residue name, residue number, chain Id and atom name.
     * @param crossLink
     *        - CrossLink object to be compared to this CrossLink object.
     * @return {@code TRUE} if both CrossLink object are equal in all atom
     *         identifier, {@code FALSE} otherwise.
     */
    public final boolean equals(final CrossLink crossLink) {
        Atom preAtom1 = this.getPreAtom();
        Atom postAtom1 = this.getPostAtom();
        Atom preAtom2 = crossLink.getPreAtom();
        Atom postAtom2 = crossLink.getPostAtom();

        if (this.equalsInHomolog(crossLink)) {
            String residueId1 = preAtom1.getChainId() + "#"
                                + preAtom1.getName().trim();
            String residueId2 = postAtom1.getChainId() + "#"
                                + postAtom1.getName().trim();
            String residueId3 = preAtom2.getChainId() + "#"
                                + preAtom2.getName().trim();
            String residueId4 = postAtom2.getChainId() + "#"
                                + postAtom2.getName().trim();

            if ((residueId1.equals(residueId3) && residueId2.equals(residueId4))
               ||
               (residueId2.equals(residueId3) && residueId1.equals(residueId4)))
            {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the ranking index of this cross-link.
     * @param rank
     *        - integer value representing the ranking index.
     */
    public final void setIndex(final int rank) {
        this.index = rank;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the ranking index of this cross-link.
     * @return Integer value representing the ranking index.
     */
    public final int getIndex() {
        return index;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of this cross-link in distance file
     * format.
     * @return String object holding the representation of this cross-link in
     * distance file format.
     */
    public final String toString() {
        String atomId1 = AminoAcid.getAminoAcidId(preAtom)
                         + "-"
                         + preAtom.getName().trim();
        String atomId2 = AminoAcid.getAminoAcidId(postAtom)
                         + "-"
                         + postAtom.getName().trim();

        if (this.solventPathDistance == Double.parseDouble(
                                                    Value.DISTANCE.getDefault())
                                                       ) {
            return this.filePath + "\t" + atomId1 + "\t" + atomId2 + "\t"
                   + this.seqDist + "\t"
                   + this.getEuclideanDistance() + "\t"
                   + "-"
                   + Constants.LINE_SEPERATOR;
        } else {
            return this.filePath + "\t" + atomId1 + "\t" + atomId2 + "\t"
                   + this.seqDist + "\t"
                   + this.getEuclideanDistance() + "\t"
                   + this.getSolventPathDistance()
                   + Constants.LINE_SEPERATOR;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to the file in which this cross-link has been found.
     * @return String object holding the path to the file.
     */
    public final String getFilePath() {
        return filePath;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the path to the file in which this cross-link has been found.
     * @param filePath
     *        - String object holding the path to the file.
     */
    public final void setFileName(String filePath) {
        this.filePath = filePath;
    }
    
}
