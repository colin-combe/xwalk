package structure.grid;

import java.util.Hashtable;

import structure.constants.Constants;
import structure.constants.Constants.Value;
import structure.math.Point3d;
import structure.math.Point3i;
import structure.matter.Atom;



/**
 * @author Abdullah Kahraman
 * @version 3.0
 */
public class GridCell {

    //-------------------------------------------------------------------------
    /**
     * Length of the grid cell edge.
     * Default {@code size = 0.5}.
     */
    private double size = Constants.DEFAULT_GRID_SIZE;

    /**
     * Diagonal of the grid cell.
     * Default {@code diagonal = Math.sqrt(0.5)}
     */
    private double diagonal = Math.sqrt(2 * Math.pow(this.size, 2));

    /**
     * Cartesian coordinates of the grid cell center.
     */
    private Point3d xyzCoordinates;
    /**
     * Indices if grid cell is located in a 3D space.
     */
    private Point3i ijkIndices;
    /**
     * Stores boolean information about whether grid cell is occupied by
     * e.g. a molecule.
     * Default is {@code FALSE}.
     */
    private boolean isOccupied = false;
    /**
     * Stores boolean information about whether grid cell lies on the 
     * boundary between occupied and unoccupied grid cells.
     * Default is {@code FALSE}.
     */
    private boolean isBoundary = false;
    /**
     * Stores boolean information about whether grid cell has been visited
     * e.g. while searching through the parent grid.
     * Default is {@code FALSE}.
     */
    private boolean hasBeenVisited = false;
    /**
     * Stores any value information as a String object.
     */
    private Hashtable < Value, String > value =
            new Hashtable < Value, String > ();

    //-------------------------------------------------------------------------
    /**
     * Constructor.
     * @param xyz
     *        - Point3d object with XYZ coordinates of the grid cell.
     * @param edgeLength
     *        - double value representing the edge length of the gridCell cube.
     */
    public GridCell(final Point3d xyz, final double edgeLength) {
        this.xyzCoordinates = xyz;
        this.setSize(edgeLength);
        this.setValue(Value.CLUSTER_NO, Value.CLUSTER_NO.getDefault());
        this.setValue(Value.DISTANCE, Value.DISTANCE.getDefault());
        this.setValue(Value.GENERAL, Value.GENERAL.getDefault());
    }
    //-------------------------------------------------------------------------
   /**
     * Sets new Cartesian coordinates for the center of this grid cell.
     * @param xyz
     *        - Point3d object holding the new Cartesian coordinates.
     */
    public final void setPoint3d(final Point3d xyz) {
        this.xyzCoordinates = xyz;
    }
    //-------------------------------------------------------------------------
   /**
     * Returns the Cartesian coordinates of the cell center of this grid cell.
     * @return Point3d object holding the Cartesian coordinates.
     */
    public final Point3d getPoint3d() {
        return this.xyzCoordinates;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets new indices for the grid cell.
     * @param ijk
     *        - Point3i object holding the indices.
     */
    public final void setPoint3i(final Point3i ijk) {
        this.ijkIndices = ijk;
    }
    //-------------------------------------------------------------------------
    /**
     * Returns the indices of this grid cell.
     * @return Point3i object holding the indices.
     */
    public final Point3i getPoint3i() {
        return this.ijkIndices;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the size, i.e. edge length of this grid cell and its diagonal.
     * @param edgeLength
     *        - double value representing the size of the grid cell.
     */
    private void setSize(final double edgeLength) {
        this.size = edgeLength;
        this.diagonal = Math.sqrt(2 * Math.pow(this.size, 2));
    }
    //-------------------------------------------------------------------------
    /**
     * Returns the size, i.e. edge length of this grid cell.
     * @return double value representing the size of the grid cell.
     */
    public final double getSize() {
        return this.size;
    }
    //-------------------------------------------------------------------------
    /**
     * Returns the diagonal of this grid cell.
     * @return double value representing the diagonal of the grid cell.
     */
    public final double getDiagonalLength() {
        return this.diagonal;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the occupation status of this grid cell to {@code TRUE}.
     */
    public final void setOccupation() {
        this.isOccupied = true;
    }
    //-------------------------------------------------------------------------

    /**
     * Sets the occupation status of this grid cell to {@code FALSE}.
     */
    public final void unsetOccupation() {
        this.isOccupied = false;
    }
    //-------------------------------------------------------------------------

    /**
     * Returns the occupation status of this grid cell.
     * @return {@code TRUE} if cell has been labeled as occupied,
     *         {@code FALSE} otherwise.
     */
    public final boolean isOccupied() {
        return this.isOccupied;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the visited status of this grid cell to {@code TRUE}.
     */
    public final void setVisitStatus() {
        this.hasBeenVisited = true;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the visited status of this grid cell to {@code FALSE}.
     */
    public final void unsetVisitStatus() {
        this.hasBeenVisited = false;
    }
    //-------------------------------------------------------------------------

    /**
     * Returns the visited status of this grid cell.
     * @return {@code TRUE} if cell has been labeled as visited,
     *         {@code FALSE} otherwise.
     */
    public final boolean isVisited() {
        return this.hasBeenVisited;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the boundary status of this grid cell to {@code TRUE}.
     */
    public final void setBoundaryStatus() {
        this.isBoundary = true;
    }
    //-------------------------------------------------------------------------
    /**
     * Sets the boundary status of this grid cell to {@code FALSE}.
     */
    public final void unsetBoundaryStatus() {
        this.isBoundary = false;
    }
    //-------------------------------------------------------------------------

    /**
     * Returns the boundary status of this grid cell.
     * @return {@code TRUE} if cell has been labeled as lying at the boundary
     *          between occupied and unoccupied cells,
     *         {@code FALSE} otherwise.
     */
    public final boolean isAtBoundary() {
        return this.isBoundary;
    }
    //-------------------------------------------------------------------------

    /**
     * Sets a String value to this grid cell, which should represent a type of
     * Value.
     * @param type
     *        - Value enumerator representing the possible value types.
     * @param valueText
     *        - String object holding the information about the value.
     * @see #getValue(Value)
     */
    public final void setValue(final Value type, final String valueText) {
        this.value.put(type, valueText);
    }
    //-------------------------------------------------------------------------

    /**
     * Returns the value of a specific type associated to this grid cell.
     * @param type
     *        - Value enumerator representing the possible value types.
     * @return String object holding the information about the value.
     * @see #setValue(Value, String)
     */
    public final String getValue(final Value type) {
        return this.value.get(type);
    }
    //-------------------------------------------------------------------------

    /**
     * Resets the distance and cluster_no value, visited and occupied status of 
     * this grid cells.
     * @see #resetSoft()
     */
    public final void reset() {
        this.setValue(Value.CLUSTER_NO, Value.CLUSTER_NO.getDefault());
        this.setValue(Value.DISTANCE, Value.DISTANCE.getDefault());
        this.unsetVisitStatus();
        this.unsetOccupation();
    }
    //-------------------------------------------------------------------------

    /**
     * Resets only the distance value and visited status of this grid cells, 
     * but not the occupied status.
     * @see #reset()
     */
    public final void resetSoft() {
        this.setValue(Value.DISTANCE, Value.DISTANCE.getDefault());
        this.unsetVisitStatus();
    }
    //-------------------------------------------------------------------------

    /**
     * Creates a copy of this GridCell object.
     * @return A new copy of this GridCell object.
     */
    public final GridCell copy() {
        GridCell copy = new GridCell(this.getPoint3d().copy(), this.getSize());
        copy.setPoint3i(this.getPoint3i().copy());
        copy.setValue(Value.CLUSTER_NO, this.getValue(Value.CLUSTER_NO));
        copy.setValue(Value.GENERAL, this.getValue(Value.GENERAL));
        copy.setValue(Value.DISTANCE, this.getValue(Value.DISTANCE));

        if (this.isOccupied) {
            copy.setOccupation();
        }
        if (this.isVisited()) {
            copy.setVisitStatus();
        }
        return copy;
    }
    //-------------------------------------------------------------------------

    /**
     * Checks whether a second GridCell object is the same as this GridCell
     * object.
     * @param cell
     *        - GridCell object to be compared to this GridCell object.
     * @return {@code TRUE} if both GridCell objects have the same XYZ
     *         Cartesian coordinates, IJK indices and sizes, {@code FALSE}
     *         otherwise.
     */
    public final boolean equals(final GridCell cell) {
        return  this.getPoint3d().equals(cell.getPoint3d())
                &&
                this.getPoint3i().equals(cell.getPoint3i())
                &&
                this.getSize() == cell.getSize() ? true : false;
    }
    //-------------------------------------------------------------------------
    /**
     * Returns the grid cell as an Atom Object.
     * @return Atom object with chainID=Y or N for occupied or unoccupied grid
     *         cells respectively and with distance information in the
     *         temperature factor column.
     */
    public final Atom toAtom() {
        Atom atom = new Atom();
        double maxTempFactorValue = Constants.MAX_OCCUPANCY_TEMPERATURE_VALUE;

        atom.setFlag("HETATM");
        
        atom.setPoint3d(new Point3d(this.getPoint3d().getX(),
                                    this.getPoint3d().getY(),
                                    this.getPoint3d().getZ()));
        try {
            double tempValue = Double.parseDouble(this.getValue(
                                                                Value.DISTANCE
                                                               )
                                                 );
            if (tempValue > maxTempFactorValue) {
                tempValue = maxTempFactorValue;
            }
            atom.setTemperatureFactor(tempValue);
        } catch (Exception e) {
              atom.setTemperatureFactor(0.0);
        }

        if (this.isOccupied()) {
           atom.setSerialNumber(Integer.parseInt(this.getValue(Value.GENERAL)));
           atom.setChainId('Y');
        } else {
            atom.setSerialNumber(0);
            atom.setChainId('N');
        }
        if (this.isAtBoundary()) {
            atom.setChainId('B');
         }
        return atom;
    }

    public String toString(){
        return this.toAtom().toString();
    }
}
