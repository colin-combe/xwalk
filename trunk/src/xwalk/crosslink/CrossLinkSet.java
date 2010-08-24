package xwalk.crosslink;

import java.util.TreeSet;

import xwalk.math.DistanceComparator;

/**
 * Container for a set of cross-link objects, where the order is determined by
 * the Solvent-Path distance or if not existent by the Euclidean distance.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class CrossLinkSet extends TreeSet < CrossLink > {

    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;

    //--------------------------------------------------------------------------
    /**
     * Constructor, which will use the DistanceComparator class for initiating
     * this CrossLinkSet.
     */
    public CrossLinkSet() {
        super(new DistanceComparator());
    }

    //--------------------------------------------------------------------------
    /**
     * Returns a String representation of all cross-links in this container.
     * @return String object holding the representation of all cross-links in
     * this container.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        int i = 1;
        for (CrossLink crossLink : this) {
            output.append(i++ + "\t" + crossLink.toString());
        }
        return output.toString();
    }
}

