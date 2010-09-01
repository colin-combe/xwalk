package structure.grid;

import java.util.ArrayList;

import structure.constants.Constants.Value;
import structure.io.pdb.PDBreader;
import structure.math.Mathematics;
import structure.math.Point3d;
import structure.math.Point3i;
import structure.math.algorithms.BoundarySearch;
import structure.matter.Atom;
import structure.matter.AtomList;


/**
 * Class to create and handle grid objects for AtomList objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 *
 */
public class GridUtilities {
    /**
     * Constructor with prevention against calls from subclass.
     */
    protected GridUtilities() {
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the number of grid cells that fit within the hemisphere.
     * @param radius
     *        - double value representing the radius of the hemisphere.
     * @param gridCellSize
     *        - double value representing the size of the grid cell.
     * @return integer value representing the number of grid cells.
     */
    public static int getNumberOfGridCellsFittingIntoHemisphere(
                                                       final double radius,
                                                       final double gridCellSize
                                                               ) {
        final double minRadius = 0.01;
        double numberOfCellsInAtom;
        if (radius == 0.0) {
            numberOfCellsInAtom = minRadius;
        } else {
            numberOfCellsInAtom = radius
                                  /
                                  gridCellSize;
        }
        return Math.round((float) numberOfCellsInAtom);
    }
    //--------------------------------------------------------------------------
    /**
     * Assigns all grid cells that same unique number if they can be connected
     * via other grid cells that have all the isOccupied() state set to
     * {@code TRUE}.
     * @param grid
     *        - Grid object which has been set up for an AtomList object.
     */
    public static void cluster(final Grid grid) {
        int clusterNo = 0;
        Point3i noOfCells = grid.getNumberOfCells();
        for (int i = 0; i < noOfCells.getI(); i++) {
            for (int j = 0; j < noOfCells.getJ(); j++) {
                for (int k = 0; k < noOfCells.getK(); k++) {
                    if (!grid.get(i, j, k).isVisited()) {
                        clusterNo++;
                        grid.get(i, j, k).setValue(Value.CLUSTER_NO,
                                                   Integer.toString(clusterNo)
                                                  );
                        grid.get(i, j, k).setVisitStatus();
                        GridUtilities.search4neighbours(i, j, k, grid,
                                                        clusterNo
                                                       );
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Recursive function to scan through the grid and assign the same cluster
     * number to all grid cells that can be connected to each other via the
     * isOccupied() state.
     * @param i
     *        - Integer value representing the X-axis position of the currently
     *          scanned grid cell within the grid.
     * @param j
     *        - Integer value representing the Y-axis position of the grid cell.
     * @param k
     *        - Integer value representing the X-axis position of the grid cell.
     * @param grid
     *        - Grid object which has been set up for an AtomList object.
     * @param clusterNo
     *        - Integer value representing the clusterNo to be assigned to this
     *          grid cell.
     */
    private static void search4neighbours(final int i, final int j, final int k,
                                          final Grid grid, final int clusterNo)
    {
        Point3i noOfCells = grid.getNumberOfCells();
        for (int m = i - 1;
             m >= 0 && m <= i + 1 && m < noOfCells.getI();
             m++) {
             for (int n = j - 1;
                  n >= 0 && n <= j + 1 && n < noOfCells.getJ();
                  n++) {
                  for (int o = k - 1;
                       o >= 0 && o <= k + 1 && o < noOfCells.getK();
                       o++) {
                       if (grid.get(m, n, o).isOccupied()) {
                           if (!grid.get(m, n, o).isVisited()) {
                               grid.get(m, n, o).setValue(Value.CLUSTER_NO,
                                                     Integer.toString(clusterNo)
                                                         );
                               grid.get(m, n, o).setVisitStatus();
                               GridUtilities.search4neighbours(m, n, o, grid,
                                                               clusterNo
                                                              );
                        }
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the volume information of a grid.
     * @param grid
     *        - Grid object from which the volume is calculated from.
     * @return Point3d object with volume X-coordinate = occupied volume,
     *         Y-coordinate = unoccupied volume, Z-coordinate = total volume.
     */
    public static Point3d getVolume(final Grid grid) {
        int hit = 0;
        int nonHit = 0;
        int count = 0;
        for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
            for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {
                    count++;
                    if (grid.get(i, j, k).isOccupied()) {
                        hit++;
                    } else {
                        nonHit++;
                    }
                }
            }
        }
        // double check that everything works fine.
        if (hit + nonHit != count) {
            System.err.println("WARNING: hit (" + hit + ") nonHit (" + nonHit
                               + ") != count(" + count + ")");
        }
        double gridCellSize = grid.get(0, 0, 0).getSize();
        double gridCellVolume = Math.pow(gridCellSize, 3);
    return new Point3d(gridCellVolume * hit,
                       gridCellVolume * nonHit,
                       gridCellVolume * count);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all cells that are unoccupied.
     * @param grid
     *        - GridObject from which all unoccupied cells will be extracted.
     * @return ArrayList of GridCells being unoccupied.
     */
    public static ArrayList < GridCell > getUnoccupiedCells(final Grid grid) {
        ArrayList < GridCell > candidates = new ArrayList < GridCell > ();
        for (int i = 0; i < grid.getNumberOfCells().getI(); i++) {
            for (int j = 0; j < grid.getNumberOfCells().getJ(); j++) {
                for (int k = 0; k < grid.getNumberOfCells().getK(); k++) {
                    GridCell cell = grid.get(i, j, k);
                    if (!cell.isOccupied()) {
                        candidates.add(cell);
                    }
                }
            }
        }
        return candidates;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns all direct neighboring grid cells.
     * @param cell
     *        - GridCell object whose neighborhood is to be determined.
     * @param grid
     *        - Grid object that holds all grid cells including the {@code cell}
     *          and all neighboring grid cells.
     * @param cellSize
     *        - integer value representing the number of shells from which
     *          neighborhood is to be extracted.
     * @return ArrayList of neighboring grid cells.
     */
    public static ArrayList < GridCell > getNeighbouringCells(
                                                            final GridCell cell,
                                                            final Grid grid,
                                                            final int cellSize
                                                             ) {
        ArrayList < GridCell > neighbours = new ArrayList < GridCell > ();
        Point3i ijk = cell.getPoint3i();

        for (int m = -cellSize;
             m <= cellSize && ijk.getI() + m >= 0
             &&
             ijk.getI() + m < grid.getNumberOfCells().getI();
             m++
             ) {
             for (int n = -cellSize;
                  n <= cellSize && ijk.getJ() + n >= 0
                  &&
                  ijk.getJ() + n < grid.getNumberOfCells().getJ();
                  n++
                  ) {
                  for (int o = -cellSize;
                       o <= cellSize && ijk.getK() + o >= 0
                       &&
                       ijk.getK() + o < grid.getNumberOfCells().getK();
                       o++
                       ) {
                       // do not add cell itself into the list of neighbours.
                       if (m == 0 && n == 0 && o == 0) {
                           continue;
                       }
                       neighbours.add(grid.get(ijk.getI() + m,
                                               ijk.getJ() + n,
                                               ijk.getK() + o)
                                              );
                }
            }
        }
        return neighbours;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether any neighbouring cells in a shell contain occupied labeled
     * cells.
     * @param atom
     *        - Atom object to be checked whether it is solvent accessible.
     * @param grid
     *        - Grid object that holds all grid cells including the {@code cell}
     *          and all neighboring grid cells. Prior to the execution of this
     *          method, the Grid object must have been searched for the boundary
     *          cells with the BoundarySearch class.
     * @return {@code TRUE} if cell is accessible, {@code FALSE} otherwise.
     * @see structure.math.algorithms.BoundarySearch
     */
    public static boolean isAccessible(final Atom atom,
                                       final AtomGrid grid) {
        ArrayList<GridCell> neighbours = grid.getAllGridCells(atom);
        for (GridCell neighbour : neighbours) {
             if (neighbour.isAtBoundary()) {
                 return true;
             }
        }
        return false;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the GridCell element that is equal to the target GridCell object.
     * @param target
     *        - GridCell object against which {@code cells} will be matched.
     * @param cells
     *        - List of GridCell objects to be matched against {@code target}.
     * @return GridCell object if found in the list of GridCell object,
     *         otherwise {@code NULL}.
     */
    public static GridCell equals(final GridCell target,
                                  final ArrayList < GridCell > cells) {
        for (GridCell cell : cells) {
            if (cell.equals(target)) {
                return cell;
            }
        }
        return null;
    }
    //--------------------------------------------------------------------------
}
