package structure.matter.hetgroups;

import structure.matter.AtomList;
import structure.matter.Molecule;
import xwalk.crosslink.CrossLinkerType;

/**
 * Class representing Cross-linker molecules.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class CrossLinker extends Molecule {
    /**
     * Default serialVersionUID.
     */
    private static final long serialVersionUID = 1L;
    /**
     * Object to store the type of this CrossLinker object.
     */
    private CrossLinkerType type;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param xlType
     *        - Type of cross-linker as defined by the CrossLinkerType enum.
     */
    public CrossLinker(final CrossLinkerType xlType) {
        super(new AtomList());
        this.type = xlType;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the type of this cross-linker.
     * @return CrossLinkerType object.
     */
    public final CrossLinkerType getType() {
        return type;
    }
}
