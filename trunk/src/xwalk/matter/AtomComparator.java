package xwalk.matter;

import java.util.Comparator;

/**
 * Comparator class that returns lexicograhical comparison between two atom
 * identifiers.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class AtomComparator implements Comparator < Atom > {

    /**
     * Compares two atom identifier lexicographically.
     * @param atom1
     *        - First Atom object.
     * @param atom2
     *        - Second Atom object.
     * @return  the value {@code 0} if the argument string is equal to
     *          this string; a value less than {@code 0} if this string
     *          is lexicographically less than the string argument; and a
     *          value greater than {@code 0} if this string is
     *          lexicographically greater than the string argument.
     */
    public final int compare(final Atom atom1, final Atom atom2) {
        String residueId1 = atom1.getResidueName() + "-"
                          + atom1.getResidueNumber() + "-"
                          + atom1.getChainId() + "-"
                          + atom1.getName();
        String residueId2 = atom2.getResidueName() + "-"
                          + atom2.getResidueNumber() + "-"
                          + atom2.getChainId() + "-"
                          + atom2.getName();
        return residueId1.compareTo(residueId2);
    }
}
