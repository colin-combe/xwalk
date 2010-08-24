package xwalk.matter.parameter;

/**
 * Supported atom enum types for cross-linking.
 * @author Abdullah
 * @version 3.0
 * @since 3.0
 *
 */
public enum AtomType {
    CARBON_ALPHA(Element.CARBON, "CA"),
    NITROGEN_EPSILON(Element.NITROGEN, "NZ");
    //--------------------------------------------------------------------------
    /**
     * Element to which atom belongs.
     */
    private Element elementType;
    /**
     * Abbreviation, e.g. PDB short-name of atom.
     */
    private String shortName;
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param element
     *        - Element object holding the element information to which atom
     *          belongs.
     * @param abbreviation
     *        - String object holding the short name of the atom, e.g. PDB
     *          short name
     */
    AtomType(final Element element, final String abbreviation) {
        this.elementType = element;
        this.shortName = abbreviation;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the associated element type of this atom type.
     * @return Element object representing element type.
     */
    public Element getElement() {
        return this.elementType;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the short name of a PDB atom of this atom type.
     * @return String object holding the short name.
     */
    public String getAbbreviation() {
        return this.shortName;
    }
}
