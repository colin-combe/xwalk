package structure.math;

import structure.constants.Constants;

/**
 * Simple class for handling Cartesian coordinates in 3 dimensional space with
 * integer numbers.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 * @see Point3d
 */
public class Point3i {
    /**
     * Cartesian X-coordinates.
     */
    private int iIndex;
    //--------------------------------------------------------------------------
    /**
     * Cartesian Y-coordinates.
     */
    private int jIndex;
    //--------------------------------------------------------------------------
    /**
     * Cartesian Z-coordinates.
     */
    private int kIndex;
    //--------------------------------------------------------------------------
    /**
     * Size of point.
     */
    private double pointRadius;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param i
     *        - Cartesian X-coordinates.
     * @param j
     *        - Cartesian Y-coordinates.
     * @param k
     *        - Cartesian Z-coordinates.
     */
    public Point3i(final int i, final int j, final int k) {
        this(i, j, k, 0);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param i
     *        - Cartesian X-coordinates.
     * @param j
     *        - Cartesian Y-coordinates.
     * @param k
     *        - Cartesian Z-coordinates.
     * @param radius
     *        - Radius of this point.
     */
    public Point3i(final int i, final int j, final int k, final double radius) {
        this.iIndex = i;
        this.jIndex = j;
        this.kIndex = k;
        this.pointRadius = radius;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the Cartesian X-coordinate of this point.
     * @return double number representing the Cartesian X-coordinate of this
     *         point.
     */
    public final int getI() {
        return iIndex;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the Cartesian Y-coordinate of this point.
     * @return double number representing the Cartesian Y-coordinate of this
     *         point.
     */
    public final int getJ() {
        return jIndex;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the Cartesian Z-coordinate of this point.
     * @return double number representing the Cartesian Z-coordinate of this
     *         point.
     */
    public final int getK() {
        return kIndex;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the radius of this point.
     * @return double number representing the radius of this point.
     */
    public final double getRadius() {
        return this.pointRadius;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether two indices have the same index values.
     * @param  point3iObject
     *              - Point3i object
     * @return {@code TRUE} if both points have the same index values,
     *         {@code FALSE} otherwise.
     * @throws ClassCastException if point3iObject is not of type Point3i.
     */
    public final boolean equals(final Object point3iObject)
                                                     throws ClassCastException {
        if (point3iObject instanceof Point3i) {
            Point3i point3i = (Point3i) point3iObject;
            if (this.iIndex == point3i.getI()
                &&
                this.jIndex == point3i.getJ()
                &&
                this.kIndex == point3i.getK()) {
                return true;
            }
        } else {
            throw new ClassCastException("ERROR: Object " + point3iObject
                                         + " is not of type Point3i"
                                         + Constants.LINE_SEPERATOR);
        }
    return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Creates a copy of this Point3i object.
     * @return Copy of this Point3i object
     */
    public final Point3i copy() {
        return new Point3i(this.getI(),
                           this.getJ(),
                           this.getK(),
                           this.getRadius());
    }
}
