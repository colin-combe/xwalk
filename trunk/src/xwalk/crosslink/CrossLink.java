package xwalk.crosslink;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import xwalk.constants.Constants;
import xwalk.constants.Constants.BondTypes;
import xwalk.constants.Constants.Value;
import xwalk.grid.Path;
import xwalk.math.Mathematics;
import xwalk.matter.Atom;
import xwalk.matter.Bond;
import xwalk.matter.pdb.AminoAcid;

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
    private int sequenceDistance;
    /**
     * Distance that the cross-link spans in Euclidean space.
     */
    private double euclideanDistance;
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
     * Calculates the distance in sequence space, which corresponds to the
     * difference of the amino acid rank to which both atoms are associated.
     */
    private void setSequenceDistance() {
        this.sequenceDistance = Math.abs(
                                         postAtom.getRank() - preAtom.getRank()
                                        );
    }
    //--------------------------------------------------------------------------

    /**
     * Calculates the distance in Euclidean space.
     */
    private void setEuclideanDistance() {
        this.euclideanDistance = Mathematics.distance(preAtom.getPoint3d(),
                                                      postAtom.getPoint3d());
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the Solvent-Path distance. In order to calculate the Solvent-Path
     * distance please see in the Class CrossLinkSet.
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
        return this.sequenceDistance;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the distance in Euclidean space.
     * @return double number representing the distance in Euclidean space.
     */
    public final double getEuclideanDistance() {
        return this.euclideanDistance;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the distance in Solvent-Path distance space.
     * @return double number representing the distance in Solvent-Path distance
     *         space.
     */
    public final double getSolventPathDistance() {
        return this.solventPathDistance;
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
     * Checks whether a second cross-link has the same atom identifier except of
     * the chain-Id as this cross-link object.
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

        String residueId1 = preAtom1.getResidueName()
                            + preAtom1.getResidueNumber();
        String residueId2 = postAtom1.getResidueName()
                            + postAtom1.getResidueNumber();
        String residueId3 = preAtom2.getResidueName()
                            + preAtom2.getResidueNumber();
        String residueId4 = postAtom2.getResidueName()
                            + postAtom2.getResidueNumber();
        if ((residueId1.equals(residueId3) && residueId2.equals(residueId4))
            ||
            (residueId2.equals(residueId3) && residueId1.equals(residueId4))) {
            return true;
        } else {
            return false;
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Returns a String representation of this cross-link in distance file
     * format.
     * @return String object holding the representation of this cross-link in
     * distance file format.
     */
    public final String toString() {
        NumberFormat decFormat = new DecimalFormat("0.0");
        String atomId1 = AminoAcid.getAminoAcidId(preAtom)
                         + "-"
                         + preAtom.getName().trim();
        String atomId2 = AminoAcid.getAminoAcidId(postAtom)
                         + "-"
                         + postAtom.getName().trim();

        if (this.solventPathDistance == Double.parseDouble(
                                                    Value.DISTANCE.getDefault())
                                                       ) {
            return atomId1 + "\t" + atomId2 + "\t"
                   + this.sequenceDistance + "\t"
                   + decFormat.format(this.euclideanDistance)
                   + Constants.LINE_SEPERATOR;
        } else {
            return atomId1 + "\t" + atomId2 + "\t"
                   + this.sequenceDistance + "\t"
                   + decFormat.format(this.euclideanDistance) + "\t"
                   + decFormat.format(this.solventPathDistance)
                   + Constants.LINE_SEPERATOR;
        }
    }
}
