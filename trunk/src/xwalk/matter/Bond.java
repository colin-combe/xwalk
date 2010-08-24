package xwalk.matter;

import xwalk.constants.Constants.BondTypes;

/**
 * Class handles chemical bonds between Atom objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Bond {
    /**
     * Atom at pre position within the bond.
     */
    private Atom pre;
    /**
     * Atom at post position within the bond.
     */
    private Atom post;
    /**
     * Type of bond.
     */
    private BondTypes bondType;
    /**
     * Constructor.
     * @param pre
     *        - Atom object of first atom.
     * @param post
     *        - Atom object of second atom.
     * @param bondType
     *        - BondType of bond.
     */
    public Bond(final Atom pre, final Atom post, final BondTypes bondType) {
        this.pre = pre;
        this.post = post;
        this.bondType = bondType;
    }
    /**
     * Returns the type of this bond.
     * @return BondType object
     */
    public final BondTypes getBondType() {
        return this.bondType;
    }

    /**
     * Returns the first atom in the bond.
     * @return Atom object.
     */
    public final Atom getPreAtom() {
        return this.pre;
    }
    /**
     * Returns the second atom in the bond.
     * @return Atom object.
     */
    public final Atom getPostAtom() {
        return this.post;
    }
    /**
     * Checks whether {@code atom} is part of the bond.
     * @return {@code TRUE} if atom is part of this bond, {@code FALSE}
     *         otherwise.
     */
    public final boolean isInBond(final Atom atom) {
        if (pre.equals(atom) || post.equals(atom)) {
            return true;
        }
        return false;
    }
}
