package xwalk.math.pdb;

import xwalk.math.Point3d;
import xwalk.matter.Atom;
import xwalk.matter.AtomList;

/**
 * Class holding various functions for transforming, i.e. translating and
 * rotating protein atoms.
 * @author Abdullah Kahraman
 * @version 3.0
 * @version 3.0
 */
public class Transformation {
    //--------------------------------------------------------------------------
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected Transformation() {
        throw new UnsupportedOperationException();
    }
    /**
     * Returns the dimension of the AtomList object on the X,Y,Z axis.
     * @param atomList
     *           - AtomList object
     * @return A Point3d object that holds the dimensions in X,Y,Z axis
     *         in the X,Y,Z fields of the Point3d object.
     * @see #min(AtomList)
     * @see #max(AtomList)
     */
    public static Point3d dimension(final AtomList atomList) {
        Point3d min = Transformation.min(atomList);
        Point3d max = Transformation.max(atomList);
        Point3d dim = new Point3d(min.getX() - max.getX(),
                                  min.getY() - max.getY(),
                                  min.getZ() - max.getZ());
    return dim;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the minimum X,Y,Z coordinates within a list of atom coordinates.
     * @param atomList
     *           - AtomList object
     * @return A Point3d object that holds the Cartesian coordinates of the
     *         minimum
     *         X,Y,Z coordinates.
     * @see #max(AtomList)
     * @see #dimension(AtomList)
     */
    public static Point3d min(final AtomList atomList) {
        double xMin = Integer.MAX_VALUE;
        double yMin = Integer.MAX_VALUE;
        double zMin = Integer.MAX_VALUE;
        for (Atom atom : atomList) {
             double x = atom.getPoint3d().getX();
             double y = atom.getPoint3d().getY();
             double z = atom.getPoint3d().getZ();
             double r = atom.getVanDerWaalsRadius();
             xMin = Math.min(xMin, x - r);
             yMin = Math.min(yMin, y - r);
             zMin = Math.min(zMin, z - r);
        }
        return new Point3d(xMin, yMin, zMin);
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the maximum X,Y,Z coordinates within a list of atom coordinates.
     * @param atomList
     *           - AtomList object
     * @return A Point3d object that holds the Cartesian coordinates of the
     *         maximum
     *         X,Y,Z coordinates.
     * @see #min(AtomList)
     * @see #dimension(AtomList)
     */
    public static Point3d max(final AtomList atomList) {
        double xMax = Integer.MIN_VALUE;
        double yMax = Integer.MIN_VALUE;
        double zMax = Integer.MIN_VALUE;
        for (Atom atom : atomList) {
             double x = atom.getPoint3d().getX();
             double y = atom.getPoint3d().getY();
             double z = atom.getPoint3d().getZ();
            double r = atom.getVanDerWaalsRadius();
            xMax = Math.max(xMax, x + r);
            yMax = Math.max(yMax, y + r);
            zMax = Math.max(zMax, z + r);
        }
        return new Point3d(xMax, yMax, zMax);
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the center of geometry of any AtomList object.
     * @param atomList
     *        - AtomList object
     * @return Point3d object that holds the Cartesian coordinates of the
     *         center of geometry of the AtomList object.
     * @see #centerOfMass(AtomList)
     */
    public static Point3d centerOfGeometry(final AtomList atomList) {
        Point3d min = Transformation.min(atomList);
        Point3d max = Transformation.max(atomList);
        return new Point3d(min.getX() + (
                                         (max.getX() - min.getX())
                                                     /
                                                     2
                                         ),
                           min.getY() + (
                                         (max.getY() - min.getY())
                                                     /
                                                     2),
                           min.getZ() + (
                                         (max.getZ() - min.getZ())
                                                     /
                                                     2
                                         )
                                        );
    }
    //--------------------------------------------------------------------------

    /**
     * Returns the center of mass of any AtomList object.
     * @param atomList
     *           - AtomList object
     * @return A Point3d object that holds the Cartesian coordinates of the
     *         center of mass of the AtomList object.
     * @see #centerOfGeometry(AtomList)
     */
    public static Point3d centerOfMass(final AtomList atomList) {
        double[] center = new double[3];
        for (Atom atom : atomList) {
             center[0] += atom.getWeight() * atom.getPoint3d().getX();
             center[1] += atom.getWeight() * atom.getPoint3d().getY();
             center[2] += atom.getWeight() * atom.getPoint3d().getZ();
        }
        center[0] /= atomList.size();
        center[1] /= atomList.size();
        center[2] /= atomList.size();
    return new Point3d(center[0], center[1], center[2]);
    }
    //--------------------------------------------------------------------------

    /**
     * Moves all atoms in an AtomList object by X,Y,Z given by the newPosition
     * Point3d object.
     * @param atomList
     *           - AtomList object
     * @param newPosition
     *           - Point3d object
     */
    public static void move(final AtomList atomList,
                            final Point3d newPosition) {
        for (Atom atom : atomList) {
             atom.setPoint3d(new Point3d(
                                  atom.getPoint3d().getX() + newPosition.getX(),
                                  atom.getPoint3d().getY() + newPosition.getY(),
                                  atom.getPoint3d().getZ() + newPosition.getZ()
                                        )
                            );
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Moves all atoms in an AtomList object such that their center of geometry
     * is aligned with the origin of the Cartesian coordinate system.
     * @param atomList
     *           - AtomList object
     * @return Point3d object holding the old center of geometry coordinates.
     *         Note that the new center of geometry is at [0,0,0].
     * @see #centerOfMass(AtomList)
     */
    public static Point3d move2centreOfGeometry(final AtomList atomList) {
        Point3d center = Transformation.centerOfGeometry(atomList);
        Transformation.move(atomList,
                            new Point3d(
                                    -center.getX(),
                                    -center.getY(),
                                    -center.getZ()
                                       )
                           );
    return center;
    }
    //--------------------------------------------------------------------------

    /**
     * Moves all atoms in an AtomList object such that their center of mass is
     * aligned with the origin of the Cartesian coordinate system.
     * @param atomList
     *           - AtomList object
     * @return Point3d object holding the old center of mass coordinates. Note
     *         that the new center of geometry is at [0,0,0].
     * @see #centerOfGeometry(AtomList)
     */
    public static Point3d move2centreOfMass(final AtomList atomList) {
        Point3d center = Transformation.centerOfMass(atomList);
        Transformation.move(atomList, new Point3d(-center.getX(),
                                                  -center.getY(),
                                                  -center.getZ()
                                                 )
                           );
    return center;
    }
}
