package xwalk.math;

import java.util.Comparator;

import xwalk.crosslink.CrossLink;

/**
 * Distance comparator class for CrossLink object, which uses the Solvent-Path
 * distance to compare two CrossLink objects or if not existent uses the
 * Euclidean distance.
 * @author Abdullah Kahraman
 * @version 3.0
 * @version 3.0
 */
public class DistanceComparator implements Comparator < CrossLink > {
    /**
     * Compares two cross-link objects by their Solvent-Path distance or if not
     * existent Euclidean distance.
     * @param xl1
     *        - First CrossLink Object.
     * @param xl2
     *        - Second CrossLink Object.
     * @return -1 if the first has a shorter distance.<br>
     *         0 if both objects span equal distances.<br>
     *         1 if the second has a larger distance.
     */
    public final int compare(final CrossLink xl1, final CrossLink xl2) {
        if (xl1.getSolventPathDistance() < xl2.getSolventPathDistance()) {
            return -1;
        }
        if (xl1.getSolventPathDistance() > xl2.getSolventPathDistance()) {
            return 1;
        }
        if (xl1.getEuclideanDistance() < xl2.getEuclideanDistance()) {
           return -1;
        }
        if (xl1.getEuclideanDistance() > xl2.getEuclideanDistance()) {
            return 1;
        }
        return 0;
    }
}
