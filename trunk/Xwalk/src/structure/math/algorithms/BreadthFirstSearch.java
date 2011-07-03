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

package structure.math.algorithms;

import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.grid.Grid;
import structure.grid.GridCell;
import structure.grid.GridUtilities;
import structure.grid.Path;
import structure.math.Mathematics;

/**
 * Set the distances in the surroundings of a source cell with a grid.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
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
     * Boolean indicating whether distance assignment had to end prematurely,
     * due to solvent inaccessibility of cross-link.
     */
    private boolean hasFinished = false;

    /**
     * Global variable to keep track on targets that have already been assigned
     * a distance.
     */
    private Hashtable < GridCell, String > targetsFoundInSearch =
                                          new Hashtable < GridCell, String >();
    //--------------------------------------------------------------------------

    /**
     * Perform breadth-first search on the grid to find the shortest path
     * between a single grid cell and a list of other grid cells.
     * @param source
     *        - Source grid cell, which represents the starting point for the
     *          distance calculation.
     * @param targets
     *        - List of target grid cells, which represent the end point in the
     *        distance calculation
     * @param grid
     *        - Grid object in which the entire search is done.
     * @param maxDist
     *        - float value representing the maximum distance to search for in
     *          the grid
     * @return List of path objects, holding each the path between the source
     *         and one target cell. If no path could be found, than each
     *         path object holds only the target cell.
     */
    public final ArrayList < Path > findShortestPath(
                                           final GridCell source,
                                           final ArrayList < GridCell > targets,
                                           final Grid grid,
                                           final float maxDist
                                                    ) {
        ArrayList < GridCell > actives = new ArrayList < GridCell >();

        // set value of source cell to 0.0
        source.setDistance(0.0f);
        actives.add(source);

        // start breadth-first search from grid cell.
        this.setDistanceRecursively(actives, targets, grid, maxDist);

        // trace back the path
        ArrayList < Path > paths = new ArrayList < Path >();
        for (GridCell target : targets) {
            // first element in the path is the target cell itself.
            // Last element will be the source cell.
            this.path.add(BreadthFirstSearch.CELL_NO_OF_TARGET_CELL_IN_PATH,
                          target.copy());
            if (target.getDistance() == Constants.DEFAULT_GRID_DISTANCE) {
                paths.add(this.path);
                path = new Path();
                continue;
            } else {
                this.backtrackPath(target, source, grid);
                paths.add(this.path);
                path = new Path();
            }
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
     *        - float value representing the maximum distance to search for in
     *          the grid
     */
    private void setDistanceRecursively(final ArrayList < GridCell > actives,
                                        final ArrayList < GridCell > targets,
                                        final Grid grid,
                                        final float maxDist) {
        // neighbouring grid cells become new actives for the next round of
        // breadth-first-search.
        ArrayList < GridCell > newActives = new ArrayList < GridCell >();
        for (GridCell active : actives) {

             ArrayList < GridCell > neighbours =
                            GridUtilities.getNeighbouringCells(active, grid, 1);

             for (GridCell neighbour : neighbours) {
                  float currentDist = Integer.MIN_VALUE;
                  float newDist = Integer.MIN_VALUE;

                  if (!neighbour.isOccupied()) {
                      currentDist = neighbour.getDistance();
                      // The distance of the neighbouring grid cell is the
                      // distance of the current active cell + the distance
                      // between the active and the neighbouring cell.
                      newDist = active.getDistance()
                              + (float) Mathematics.distance(   active.getXYZ(),
                                                             neighbour.getXYZ()
                                                           );
                      // distance from other active cell might be shorter.
                      if (newDist < currentDist) {
                          neighbour.setDistance(newDist);
                          if (!newActives.contains(neighbour)) {
                              newActives.add(neighbour);
                          }
                      }
                  }
                  // check whether we have reached the target cell already.
                  // If so, remove target from list of targets.
                  GridCell equal = GridUtilities.equals(neighbour, targets);
                  if (equal != null) {
                      this.targetsFoundInSearch.put(neighbour, "");
                      if (targets.size() == this.targetsFoundInSearch.size()) {
                          this.hasFinished = true;
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
             if (neighbour.getDistance() > maxDist) {
                 maxDistCount++;
                 this.hasFinished = true;
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
        float minDist = target.getDistance();
        GridCell minGridCell = target;
        for (GridCell neighbour : neighbours) {
             float dist = neighbour.getDistance();
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
    //--------------------------------------------------------------------------
    /**
     * Returns whether the search for a target was successful, hence could
     * be finished.
     * @return {@code TRUE} if search found target cell, {@code FALSE}
     * otherwise.
     */
    public final boolean hasFinished() {
        return this.hasFinished;
    }
}
