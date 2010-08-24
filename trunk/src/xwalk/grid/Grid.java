package xwalk.grid;

import xwalk.math.Point3d;
import xwalk.math.Point3i;
import xwalk.matter.Atom;

/**
 * Class for handling general grid objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Grid {

    /**
     * Grid cells that constitute the grid.
     */
    private GridCell[][][] gridCells;

    /**
     * Stores the maximum Cartesian coordinates of the grid.
     */
    private Point3d max;

    /**
     * Stores the minimum Cartesian coordinates of the grid.
     */
    private Point3d min;

    /**
     * Stores the number of cells in the X dimension.
     */
    private int noOfxCells = -1;
    /**
     * Stores the number of cells in the Y dimension.
     */
    private int noOfyCells = -1;
    /**
     * Stores the number of cells in the ~ dimension.
     */
    private int noOfzCells = -1;

    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param minimum
     *        - Point3d object holding the minimum X,Y,Z Cartesian coordinates
     *          of this grid.
     * @param maximum
     *        - Point3d object holding the maximum X,Y,Z Cartesian coordinates
     *          of this grid.
     * @param gridCellSize
     *        - Double value representing the size of all grid cells i.e. their
     *          cell edge length.
     */
    public Grid(final Point3d minimum,
                final Point3d maximum,
                final double gridCellSize) {
        this.min = minimum;
        this.max = maximum;
        this.setNumberOfCells(gridCellSize);
        this.gridCells = new GridCell[this.noOfxCells]
                                     [this.noOfyCells]
                                     [this.noOfzCells];

        for (int i = 0; i < this.noOfxCells; i++) {
            for (int j = 0; j < this.noOfyCells; j++) {
                for (int k = 0; k < this.noOfzCells; k++) {
                     double x = this.min.getX() + (i * gridCellSize)
                                                           + (gridCellSize / 2);
                     double y = this.min.getY() + (j * gridCellSize)
                                                           + (gridCellSize / 2);
                     double z = this.min.getZ() + (k * gridCellSize)
                                                           + (gridCellSize / 2);

                    this.gridCells[i][j][k] = new GridCell(new Point3d(x, y, z),
                                                           gridCellSize);
                    this.gridCells[i][j][k].setPoint3i(new Point3i(i, j, k));
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    /**
     * Returns the minimum Cartesian coordinate of this grid.
     * @return Point3d
     *        - Point3d object, which stores the minimum Cartesian coordinate.
     */
    protected final Point3d getMin() {
        return this.min;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the maximum Cartesian coordinate of this grid.
     * @return Point3d
     *        - Point3d object, which stores the maximum Cartesian coordinate.
     */
    protected final Point3d getMax() {
        return this.max;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates the number of cells in all three XYZ Cartesian dimensions.
     * @param gridCellSize
     *        - Double value representing the size of all grid cells i.e. their
     *          cell edge length.
     */
    private void setNumberOfCells(final double gridCellSize) {
       this.noOfxCells = Math.round((float) ((this.max.getX() - this.min.getX())
                                             /
                                             gridCellSize
                                            )
                                   );
       this.noOfyCells = Math.round((float) ((this.max.getY() - this.min.getY())
                                             /
                                             gridCellSize
                                            )
                                   );
       this.noOfzCells = Math.round((float) ((this.max.getZ() - this.min.getZ())
                                             /
                                             gridCellSize
                                            )
                                   );
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the number of cells for each Cartesian dimension.
     * @return Point3i object holding the number of cells for each Cartesian
     *         dimension.
     */
    public final Point3i getNumberOfCells() {
        return new Point3i(this.noOfxCells, this.noOfyCells, this.noOfzCells);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a grid cell with the grid.
     * @param i
     *        - integer value representing the index on the X dimension.
     * @param j
     *        - integer value representing the index on the Y dimension.
     * @param k
     *        - integer value representing the index on the Z dimension.
     * @return GridCell object with indices i,j,k. If indices extend over grid
     *         borders, then {@code NULL} is returned.
     */
    public final GridCell get(final int i, final int j, final int k) {
        boolean indicesLarger0 = i >= 0 && j >= 0 && k >= 0;
        boolean indicesSmallerMax = i < this.noOfxCells
                                    &&
                                    j < this.noOfyCells
                                    &&
                                    k < this.noOfzCells;
        if (indicesLarger0 && indicesSmallerMax) {
            return this.gridCells[i][j][k];
        } else {
            return null;
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Resets the value, visited and occupied status of ALL grid cells in the
     * grid.
     */
    public final void reset() {

        for (int i = 0; i < noOfxCells; i++) {
            for (int j = 0; j < noOfyCells; j++) {
                for (int k = 0; k < noOfzCells; k++) {
                    this.gridCells[i][j][k].reset();
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Resets only value and visited, leaving occupied status as it is.
     */
    public final void resetSoft() {
        for (int i = 0; i < noOfxCells; i++) {
            for (int j = 0; j < noOfyCells; j++) {
                for (int k = 0; k < noOfzCells; k++) {
                    this.gridCells[i][j][k].resetSoft();
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the grid cell points in PDB format.
     * @return String object holding the grid in PDB format.
     */
    public final String toString() {
        StringBuffer output = new StringBuffer();
        for (int i = 0; i < noOfxCells; i++) {
             for (int j = 0; j < noOfyCells; j++) {
                  for (int k = 0; k < noOfzCells; k++) {
                       GridCell cell = this.get(i, j, k);
                       Atom dummy = cell.toAtom();
                       dummy.setName("C");
                       dummy.setResidueName("GRD");
                       dummy.setResidueNumber(1);
                       output.append(dummy.toString());
                }
            }
        }
        return output.toString();
    }
}
