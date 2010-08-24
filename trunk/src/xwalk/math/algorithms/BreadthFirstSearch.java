package xwalk.math.algorithms;

import java.util.ArrayList;
import java.util.Hashtable;

import xwalk.constants.Constants.Value;
import xwalk.grid.Grid;
import xwalk.grid.GridCell;
import xwalk.grid.GridUtilities;
import xwalk.grid.Path;
import xwalk.math.Mathematics;
/**
 * Set the distances in the surroundings of a source cell with a grid.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class BreadthFirstSearch {
    /**
     * Shortest path as calculated between a source and a target cell.
     */
    private Path path = new Path();

    /**
     * Constant, indicating that the target cell is the first in the path.
     */
    public static final int CELL_NO_OF_TARGET_CELL_IN_PATH = 0;
    /**
     *
     */
    private Hashtable < GridCell, String > targetsFoundInSearch =
                                          new Hashtable < GridCell, String > ();
    //--------------------------------------------------------------------------

    /**
     * Perform breadth-first search on the grid to find the shortest path
     * between a single grid cell and a list of other grid cells.
     * @param source
     *        - Source grid cell, which represents the starting point for the
     *          distance calculation.
     * @param targets
     *        - List of target grid cells, which represent the end point in the
     *        distance calcuation
     * @param grid
     *        - Grid object in which the entire search is done.
     * @param maxDist
     *        - double value representing the maximum distance to search for in
     *          the grid
     * @return List of path objects, holding each the path between the source
     *         and one target cell.
     */
    public final ArrayList < Path > findShortestPath(
                                           final GridCell source,
                                           final ArrayList < GridCell > targets,
                                           final Grid grid,
                                           final double maxDist
                                                    ) {
        ArrayList < GridCell > actives = new ArrayList < GridCell > ();

        // set value of source cell to 0.0
        source.setValue(Value.DISTANCE, "0.0");
        actives.add(source);
        
        // start breadth-first search from grid cell.
        this.setDistanceRecursively(actives, targets, grid, maxDist);

        // trace back the path
        ArrayList < Path > paths = new ArrayList < Path > ();
        for (GridCell target : targets) {
            // first element in the path is the target cell itself.
            // Last element will be the source cell.
            this.path.add(BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH,
                          target.copy());
            if (target.getValue(Value.DISTANCE).equals(
                                                     Value.DISTANCE.getDefault()
                                                      )
                                     ) {
                paths.add(this.path);
                path = new Path();
                continue;
            }
            this.backtrackPath(target, source, grid);
            paths.add(this.path);
            path = new Path();
        }
        return paths;
    }
    //--------------------------------------------------------------------------

    /**
     * Sets the distances for all grid cells that lie in-between the source cell
     * and a list of target cells.
     * @param actives
     *        - List of GridCells for which neighbouring cells have to be
     *          determined and distances calculated.
     * @param targets
     *        - List of GridCell objects that represent the target grid cell to
     *          be reached.
     * @param grid
     *        - Grid object in which the entire search is done.
     * @param maxDist
     *        - double value representing the maximum distance to search for in
     *          the grid
     */
    private void setDistanceRecursively(final ArrayList < GridCell > actives,
                                        final ArrayList < GridCell > targets,
                                        final Grid grid,
                                        final double maxDist) {
        // neighbouring grid cells become new actives for the next round of
        // breadth-first-search.
        ArrayList < GridCell > newActives = new ArrayList < GridCell >();
        for (GridCell active : actives) {
            
             ArrayList < GridCell > neighbours =
                            GridUtilities.getNeighbouringCells(active, grid, 1);

             for (GridCell neighbour : neighbours) {
                  double currentDist = Integer.MIN_VALUE;
                  double newDist = Integer.MIN_VALUE;

                  if (!neighbour.isOccupied()) {
                      // The distance of the neighbouring grid cell is the
                      // distance of the current active cell + the distance
                      // between the active and the neighbouring cell.
                      currentDist = Double.parseDouble(
                                              neighbour.getValue(Value.DISTANCE)
                                                      );
                      newDist = Double.parseDouble(active.getValue(
                                                                 Value.DISTANCE)
                                                                  )
                              + Mathematics.distance(active.getPoint3d(),
                                                     neighbour.getPoint3d()
                                                    );
                      // distance from other active cell might be shorter.
                      if (newDist < currentDist) {
                          neighbour.setValue(Value.DISTANCE, newDist + "");
                          if (!newActives.contains(neighbour)) {
                              newActives.add(neighbour);
                          }
                      }
                  }
                  // remove target from list of targets if is has been found.
                  GridCell equal = GridUtilities.equals(neighbour, targets);
                  if (equal != null) {
                      this.targetsFoundInSearch.put(neighbour, "");
                      if (targets.size() == this.targetsFoundInSearch.size()) {
                          return;
                      }
                  }
             }
        }
        // set the visit flag in all new active cells to true to avoid
        // recalculation of distances for these cells.
        // Check furthermore whether it is necessary to continue distance
        // calculation, as all newActives might have already distances larger
        // than maxDist.
        int maxDistCount = 0;
        for (GridCell neighbour : newActives) {
             neighbour.setVisitStatus();
             double neighbourDistance = Double.parseDouble(
                                               neighbour.getValue(Value.DISTANCE
                                                                 )
                                                           );
             if (neighbourDistance > maxDist) {
                 maxDistCount++;
             }
        }
        // break up recursive loop.
        if (maxDistCount == newActives.size()) {
            return;
        }
        this.setDistanceRecursively(newActives, targets, grid, maxDist);
    }
    //--------------------------------------------------------------------------
    /**
     * Backtraces the path starting between a target cell and source cell.
     * @param source
     *        - Source grid cell, which represents the starting point for the
     *          distance calculation.
     * @param target
     *        - Target grid cell, which represent the end point in the
     *        distance calculation.
     * @param grid
     *        - Grid object in which the entire search is done.
     */
    private void backtrackPath(final GridCell target,
                               final GridCell source,
                               final Grid grid) {
        ArrayList < GridCell > neighbours =
                            GridUtilities.getNeighbouringCells(target, grid, 1);
        double minDist = Double.parseDouble(target.getValue(Value.DISTANCE));
        GridCell minGridCell = target;
        for (GridCell neighbour : neighbours) {
             double dist = Double.parseDouble(
                                              neighbour.getValue(Value.DISTANCE)
                                             );
            if (dist < minDist) {
                minDist = dist;
                minGridCell = neighbour;
            }
        }
        if (minGridCell == source) {
            this.path.add(minGridCell.copy());
            return;
        }
        this.path.add(minGridCell.copy());
        this.backtrackPath(minGridCell, source, grid);
    }
}
