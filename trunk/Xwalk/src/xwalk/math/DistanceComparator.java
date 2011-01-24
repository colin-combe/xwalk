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

package xwalk.math;

import java.util.Comparator;

import xwalk.crosslink.CrossLink;

/**
 * Distance comparator class for CrossLink object, which uses the Solvent-Path
 * distance to compare two CrossLink objects or if not existent uses the
 * Euclidean distance.
 * @author Abdullah Kahraman
 * @version 0.1
 * @version 0.1
 */
public class DistanceComparator implements Comparator < CrossLink > {
    /**
     * Compares two cross-link objects by their Probabilities to occur in
     * experiments or if not assessed by their Solvent-Path distance or if not
     * existent Euclidean distance.
     * @param xl1
     *        - First CrossLink Object.
     * @param xl2
     *        - Second CrossLink Object.
     * @return -1 if the first has a shorter distance/higher probability.<br>
     *         0 if both objects span equal distances/probabilities.<br>
     *         1 if the second has a larger distance/lower probability
     */
    public final int compare(final CrossLink xl1, final CrossLink xl2) {
        if (xl1.getSolventPathDistanceProbability()
                <
            xl2.getSolventPathDistanceProbability()) {
            return 1;
        }
        if (xl1.getSolventPathDistanceProbability()
                >
            xl2.getSolventPathDistanceProbability()) {
            return -1;
        }
        if (xl1.getEuclideanDistanceProbability()
                <
            xl2.getEuclideanDistanceProbability()) {
            return 1;
        }
        if (xl1.getEuclideanDistanceProbability()
                >
            xl2.getEuclideanDistanceProbability()) {
            return -1;
        }
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
