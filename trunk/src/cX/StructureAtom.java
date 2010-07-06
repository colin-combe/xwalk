package cX;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.Locale;

/**
 * Class that defines all the model fields of a standard PDB file.
 * Functions are provided to check and set every entry individually. 
 * The only checking done concerns the allowed format and range of 
 * the numbers, chemical validity must be checked elsewhere.
 * The function toString returns the entry of one atom correctly 
 * formatted for output (multiple conformations must however be taken 
 * care of in the calling routines and are not checked for here).
 * PDB reading (and model writing) is handling in StructureCoords.
 */
public class StructureAtom {

	/*--------------------------------------------------------------------------*/
	// Classobjects
	/*--------------------------------------------------------------------------*/    

    private String flag;       //  1 -  6
    private int serial;        //  7 - 11
    private String name;       // 13 - 16
    private char altLoc;       // 17
    private String resName;    // 18 - 20
    private char chainID;      // 22  
    private int resNo;        // 23 - 26 
    private char iCode;        // 27
    private double x;          // 31 - 38
    private double y;          // 39 - 46
    private double z;          // 47 - 54
    private double occupancy;  // 55 - 60
    private double tempFactor; // 61 - 66
    private String segID;      // 73 - 76
    private String element;    // 77 - 78
    private double potential;  // 79 - 80
    private String remark;      // any kind of information about atom

    static private String ff	= "mmff";   	   // which forcefield to use to obtain atomic radii 

    private StructureAtom covalentBondAtom;
	// To account for resolution error distance is taken to be less than 1.8 (0.28  estimated standard error for X-ray structures)
    private static double covalentRadiusErrorRange = 0.28;
    
    // ANISOU
    private int u00;           // 2nd Line 29 - 35
    private int u11;           // 2nd Line 36 - 42 
    private int u22;           // 2nd Line 43 - 49
    private int u01;           // 2nd Line 50 - 56
    private int u02;           // 2nd Line 57 - 63
    private int u12;           // 2nd Line 64 - 70
    
	// SPHERICAL COORDINATES
	private double r;
	private double phi;
	private double theta;
	
	// WEIGHT OF ATOM
	private double weight;
	
	// RADIUS OF ATOM
	private double radius;

	// FOR CONSURF
	private int colour;
	
	// FOR ORDERING ATOM GROUPS TO CLUSTER, EG. SIZE CLUSTERING OF CLEFTS
	private int clusterPos;

	// FOR ATOMIC HYDROPHOBICITY CALCULATION OF PROTEINS
	private double solvationParameter;

	// FOR VAN DER WAALS CALCULATION
//	private String eccpAtomType;
	
	// FOR VAN DER WAALS CALCULATION
	private String mmffAtomType;

	boolean isAromatic = false;
		
	// FOR HYDROGEN BOND POTENTIAL CALCULATION
	private int hybridizationState;

	// for physicochemical property calculation
	private double partialCharge;
	private double xlogP;
	
	/*--------------------------------------------------------------------------*/
	// Constructors
	/*--------------------------------------------------------------------------*/    

    /**
     * Constructor sets all the fields to defaults (mainly 0 and "") 
     * but enough to enable printing by later just parsing x,y,z
     */ 
    public StructureAtom () {
		// a few defauilts to speed up parsing
		flag = "ATOM  "; 
		serial = 0;      
		name = "";    
		altLoc = ' ';     
		resName = ""; 
		chainID = ' ';   
		resNo = 0;     
		iCode = ' ';     
		x = 0.0;        
		y = 0.0;        
		z = 0.0;
		occupancy = 0.0; 
		tempFactor = 0.0;
		segID = "";     
		element = "";  
		potential = 0.0; 

		remark = "";   
		r = 0.0;
		theta = 0.0;         
		phi = 0.0;
	
		weight = 1.0;
		radius = 1.5;
		colour = 0; 
		clusterPos = -1;
	
		solvationParameter = 0.0;
		mmffAtomType = "";
		hybridizationState = 0;
		
		partialCharge = 0;
		xlogP = 0;
			
		// ANISOU
		u00 = 0;        
		u11 = 0;        
		u22 = 0;        
		u01 = 0;        
		u02 = 0;        
		u12 = 0;
    }

	/*--------------------------------------------------------------------------*/
	// Classmethods
	/*--------------------------------------------------------------------------*/    

	/**
	 * clone atom
	 * @param -none 
	 * @return - none
	 */
	public StructureAtom copy(){
		StructureAtom atom = new StructureAtom();
		atom.setFlag(new String(this.getFlag()));
		atom.setSerial((new Integer(this.getSerial())).intValue());
		atom.setName(new String(this.getName()));
		atom.setAltLoc(new Character(this.getAltLoc()).charValue());
		atom.setResName(new String(this.getResName()));
		atom.setChainId(new Character(this.getChainId()).charValue());
		atom.setResNo((new Integer(this.getResNo())).intValue());
		atom.setICode(new Character(this.getICode()).charValue());
		atom.setX((new Double(this.getX())).doubleValue());
		atom.setY((new Double(this.getY())).doubleValue());
		atom.setZ((new Double(this.getZ())).doubleValue());
		atom.setOccupancy((new Double(this.getOccupancy())).doubleValue());
		atom.setTemp((new Double(this.getTemp())).doubleValue());
		atom.setSegID(new String(this.getSegID()));
		atom.setElement(new String(this.getElement()));
		atom.setPotential((new Double(this.getPotential())).doubleValue());
	
		atom.setRemark(new String(this.getRemark()));
		atom.setR((new Double(this.getR())).doubleValue());
		atom.setTheta((new Double(this.getTheta())).doubleValue());
		atom.setPhi(new Double(this.getPhi()).doubleValue());

		atom.setWeight(new Double(this.getWeight()).doubleValue());
		atom.setColour(new Integer(this.getColour()).intValue());
		atom.setRadius((new Double(this.getRadius())).doubleValue());
		atom.setClusterPos(new Integer(this.getClusterPos()).intValue());
				
		atom.setSolvParam(new Double(this.getSolvParam()).doubleValue());
		atom.setMMFFatomType(new String(this.getMMFFatomType()));
		
		atom.setPartialCharge(new Double(this.getPartialCharge()));
		atom.setXlogP(new Double(this.getXlogP()));
		
		if(this.isAromatic())atom.setAromatic();

		atom.setHybridizationState(new Integer(this.getHybridizationState()));
		
		return atom;
	}

	/*--------------------------------------------------------------------------*/
	public void setCovalentBondAtom(StructureAtom atom){
		this.covalentBondAtom = atom;
	}
	/*--------------------------------------------------------------------------*/
	public StructureAtom getCovalentBondAtom(){
		return this.covalentBondAtom;
	}
	/*--------------------------------------------------------------------------*/
	private void setCHARMMradius(){
		// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
		if (this.getFlag().equals("ATOM  ") && this.getName().length() > 0){
			char atomSort = this.getName().trim().charAt(0);
			switch(atomSort){

			case 'C': {
				this.setRadius(2.175); 
				break;
			}
			case 'N': {
				this.setRadius(1.850);
				break;
			}
			case 'O': {
				this.setRadius(1.700);
				break;
			}
			case 'S': {
				this.setRadius(2.000);
				break;
			}
			case 'H': {
				this.setRadius(1.100);
				break;
			}
			
			// in some NMR structures a hydrogen is represented by 
				// the letter "Q"
				case 'Q': {
					this.setRadius(1.100);
					break;
				}
				default: this.setRadius(1.500);
			}
		}
		
		if (this.getFlag().equals("HETATM") && this.getName().length() > 0){
			String atomName = this.getName().trim();
			if(atomName.equalsIgnoreCase("CA")){ 
				this.setRadius(1.367);
			}
			else if(atomName.equalsIgnoreCase("FE")){ 
				this.setRadius(0.65);
			}
			else if(atomName.equalsIgnoreCase("ZN")){ 
				this.setRadius(1.090);
			}
			else if(atomName.equalsIgnoreCase("MG")){ 
				this.setRadius(1.185);
			}
			else if(atomName.equalsIgnoreCase("MN")){ 
				this.setRadius(1.367);
			}
			else if(atomName.equalsIgnoreCase("LU")){ 
				this.setRadius(1.367);
			}
			else if(atomName.equalsIgnoreCase("K")){ 
				this.setRadius(1.764);
			}
			else if(atomName.equalsIgnoreCase("CD")){ 
				this.setRadius(1.367);
			}
			else if(atomName.equalsIgnoreCase("NA")){ 
				this.setRadius(1.364);
			}
			else if(atomName.equalsIgnoreCase("I")){ 
				this.setRadius(2.270);
			}
			else if(atomName.equalsIgnoreCase("U")){ 
				this.setRadius(1.364);
			}
			else if(atomName.startsWith("C")) {
				this.setRadius(1.90); 
				}
			else if(atomName.startsWith("N")) {
				this.setRadius(1.850);
				}
			else if(atomName.startsWith("O")) {
				this.setRadius(1.700);
			}
			else if(atomName.startsWith("S")) {
				this.setRadius(2.000);
			}
			else if(atomName.startsWith("P")) {
				this.setRadius(2.150);
			}
			else if(atomName.startsWith("H")) {  // taken from NAD molecule
				this.setRadius(1.100);
			}
			else{
				this.setRadius(1.500);
			}

			if(this.covalentBondAtom!=null){ 
				if(atomName.startsWith("H") && this.covalentBondAtom.getName().trim().startsWith("C")) {  
					this.setRadius(1.320);
				}
				else if(atomName.startsWith("H") && this.covalentBondAtom.getName().trim().startsWith("N")) {  
					this.setRadius(0.225);
				}
				else if(atomName.startsWith("H") && this.covalentBondAtom.getName().trim().startsWith("O")) {  
					this.setRadius(0.225);
				}
				else if(atomName.startsWith("H") && this.covalentBondAtom.getName().trim().startsWith("S")) {  
					this.setRadius(0.450);
				}
			}
		}
	}

	/*--------------------------------------------------------------------------*/
	private void setMMFFradius(){
		// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
		if (this.getFlag().equals("ATOM  ") && this.getName().length() > 0){
			char atomSort = this.getName().trim().charAt(0);
			switch(atomSort){

			case 'C': {
				this.setRadius(1.969); 
				break;
			}
			case 'N': {
				this.setRadius(2.014);
				break;
			}
			case 'O': {
				this.setRadius(1.779);
				break;
			}
			case 'S': {
				this.setRadius(2.185);
				break;
			}
			case 'H': {
				this.setRadius(1.307);
				break;
			}
			
			// in some NMR structures a hydrogen is represented by 
				// the letter "Q"
				case 'Q': {
					this.setRadius(1.307);
					break;
				}
				default: this.setRadius(1.500);
			}
		}
		
		if (this.getFlag().equals("HETATM") && this.getName().length() > 0){
			String atomName = this.getName().trim();
			if(atomName.equalsIgnoreCase("CA")){ 
				this.setRadius(1.948);
			}
			else if(atomName.equalsIgnoreCase("FE")){ 
				this.setRadius(1.638);
			}
			else if(atomName.equalsIgnoreCase("ZN")){ 
				this.setRadius(1.620);
			}
			else if(atomName.equalsIgnoreCase("MG")){ 
				this.setRadius(1.538);
			}
			else if(atomName.equalsIgnoreCase("MN")){ 
				this.setRadius(1.638);
			}
			else if(atomName.equalsIgnoreCase("LU")){ 
				this.setRadius(1.948);
			}
			else if(atomName.equalsIgnoreCase("K")){ 
				this.setRadius(2.000);
			}
			else if(atomName.equalsIgnoreCase("CD")){ 
				this.setRadius(1.620);
			}
			else if(atomName.equalsIgnoreCase("NA")){ 
				this.setRadius(1.591);
			}
			else if(atomName.equalsIgnoreCase("I")){ 
				this.setRadius(2.358);
			}
			else if(atomName.equalsIgnoreCase("U")){ 
				this.setRadius(1.538);
			}
			else if(atomName.startsWith("C")) {
				this.setRadius(1.969); 
				}
			else if(atomName.startsWith("N")) {
				this.setRadius(2.014);
				}
			else if(atomName.startsWith("O")) {
				this.setRadius(1.779);
			}
			else if(atomName.startsWith("S")) {
				this.setRadius(2.185);
			}
			else if(atomName.startsWith("P")) {
				this.setRadius(2.287);
			}
			else if(atomName.startsWith("H")) {  // taken from NAD molecule
				this.setRadius(1.307);
			}
			else{
				this.setRadius(1.500);
			}
		}
	}
	
	/*--------------------------------------------------------------------------*/
	public void setVdwPlusHydrogenRadius(){
		// taken from SURFNET
		// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
		if (this.getFlag().equals("ATOM  ") && this.getName().length() > 0){
			char atomSort = this.getName().trim().charAt(0);
			switch(atomSort){
				case 'C': {
					this.setRadius(1.87); 
					break;
				}
				case 'N': {
					this.setRadius(1.65);
					break;
				}
				case 'O': {
					this.setRadius(1.40);
					break;
				}
				case 'S': {
					this.setRadius(1.85);
					break;
				}
				case 'H': {
					this.setRadius(1.20);
					break;
				}
	
				// in some NMR structures a hydrogen is represented by 
				// the letter "Q"
				case 'Q': {
					this.setRadius(1.200);
					break;
				}
				default: this.setRadius(1.500);
			}
		}
		
		if (this.getFlag().equals("HETATM") && this.getName().length() > 0){
			String atomName = this.getName().trim();
			if(atomName.equalsIgnoreCase("FE")){
				this.setRadius(1.10);
			}
			else if(atomName.startsWith("C")) {
				this.setRadius(1.87); 
			}
			else if(atomName.startsWith("N")) {
				this.setRadius(1.65);
			}
			else if(atomName.startsWith("O")) {
				this.setRadius(1.40);
			}
			else if(atomName.startsWith("S")) {
				this.setRadius(1.85);
			}
			else if(atomName.startsWith("P")) {
				this.setRadius(1.90);
			}
			else if(atomName.startsWith("H")) {
				this.setRadius(1.200);
			}
			
			else{
				this.setRadius(1.500);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------//
	private void setRasmolRadius(){
		// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
		if (this.getFlag().equals("ATOM  ") && this.getName().length() > 0){
			char atomSort = this.getName().trim().charAt(0);
			switch(atomSort){
 				case 'C': {
					this.setRadius(1.548); 
					break;
				}
				case 'N': {
					this.setRadius(1.400);
					break;
				}
				case 'O': {
					this.setRadius(1.348);
					break;
				}
				case 'S': {
					this.setRadius(1.808);
					break;
				}
				case 'H': {
					this.setRadius(1.100);
					break;
				}
			
				// in some NMR structures a hydrogen is represented by 
				// the letter "Q"
				case 'Q': {
					this.setRadius(1.100);
					break;
				}
				default: this.setRadius(1.500);
			}
		}
		
		if (this.getFlag().equals("HETATM") && this.getName().length() > 0){
			String atomName = this.getName().trim();
			if(atomName.equalsIgnoreCase("CA")){
				this.setRadius(1.948);
			}
			else if(atomName.equalsIgnoreCase("FE")){
				this.setRadius(1.948);
			}
			else if(atomName.equalsIgnoreCase("ZN")){
				this.setRadius(1.148);
			}
			else if(atomName.equalsIgnoreCase("CD")){
				this.setRadius(1.748);
			}
			else if(atomName.equalsIgnoreCase("I")){
				this.setRadius(1.748);
			}
			else if(atomName.startsWith("C")) {
				this.setRadius(1.548); 
			}
			else if(atomName.startsWith("N")) {
				this.setRadius(1.400);
			}
			else if(atomName.startsWith("O")) {
				this.setRadius(1.348);
			}
			else if(atomName.startsWith("S")) {
				this.setRadius(1.808);
			}
			else if(atomName.startsWith("P")) {
				this.setRadius(1.880);
			}
			else if(atomName.startsWith("H")) {
				this.setRadius(1.100);
			}
		}
	}

	//------------------------------------------------------------------------------------------------------------------//
	private void setPARSEradius(){
		// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
		if (this.getFlag().equals("ATOM  ") && this.getName().length() > 0){
			char atomSort = this.getName().trim().charAt(0);
			switch(atomSort){
 				case 'C': {
					this.setRadius(2.000); 
					break;
				}
				case 'N': {
					this.setRadius(1.500);
					break;
				}
				case 'O': {
					this.setRadius(1.400);
					break;
				}
				case 'S': {
					this.setRadius(1.850);
					break;
				}
				case 'H': {
					this.setRadius(1.000);
					break;
				}
				// in some NMR structures a hydrogen is represented by 
				// the letter "Q"
				case 'Q': {
					this.setRadius(1.000);
					break;
				}
				default: this.setRadius(1.500);
			}
		}
		
		if (this.getFlag().equals("HETATM") && this.getName().length() > 0){
			String atomName = this.getName().trim();
			if(atomName.equalsIgnoreCase("CA")){	// PAULING univalent radius
				this.setRadius(1.18);
			}
			else if(atomName.equalsIgnoreCase("ZN")){ // PAULING univalent radius
				this.setRadius(0.88);
			}
			else if(atomName.equalsIgnoreCase("MG")){ // PAULING univalent radius
				this.setRadius(0.82);
			}
			else if(atomName.equalsIgnoreCase("K")){ // PAULING univalent radius
				this.setRadius(1.33);
			}
			else if(atomName.equalsIgnoreCase("CD")){ // PAULING univalent radius
				this.setRadius(1.14);
			}
			else if(atomName.equalsIgnoreCase("NA")){ // PAULING univalent radius
				this.setRadius(0.95);
			}
			else if(atomName.equalsIgnoreCase("LI")){ // PAULING univalent radius
				this.setRadius(0.60);
			}
			else if(atomName.equalsIgnoreCase("SI")){ // PAULING univalent radius
				this.setRadius(0.65);
			}
			else if(atomName.equalsIgnoreCase("CU")){ // PAULING univalent radius
				this.setRadius(0.96);
			}
			else if(atomName.equalsIgnoreCase("SE")){ // PAULING univalent radius
				this.setRadius(0.66);
			}
			else if(atomName.equalsIgnoreCase("FE")){ // PAULING empirical crystal radius
				this.setRadius(0.65);
			}
			else if(atomName.equalsIgnoreCase("MN")){ // PAULING empirical crystal radius
				this.setRadius(0.80);
			}
			else if(atomName.equalsIgnoreCase("LU")){ // PAULING empirical crystal radius
				this.setRadius(0.93);
			}
			else if(atomName.equalsIgnoreCase("U")){ // PAULING empirical crystal radius
				this.setRadius(1.11);
			}
			else if(atomName.equalsIgnoreCase("I")){ // PAULING vdW radius
				this.setRadius(2.16);
			}
			else if(atomName.equalsIgnoreCase("F")){ // PAULING vdW radius
				this.setRadius(1.36);
			}
			else if(atomName.equalsIgnoreCase("CL")){ // PAULING vdW radius
				this.setRadius(1.81);
			}
			else if(atomName.equalsIgnoreCase("BR")){ // PAULING vdW radius
				this.setRadius(1.95);
			}
			else if(atomName.startsWith("P")) { // PAULING vdW radius
				this.setRadius(1.9);
			}			
			else if(atomName.startsWith("C")) { // PARSE
				this.setRadius(2.000); 
				}
			else if(atomName.startsWith("N")) { // PARSE
				this.setRadius(1.500);
			}
			else if(atomName.startsWith("O")) { // PARSE
				this.setRadius(1.400);
			}
			else if(atomName.startsWith("S")) { // PARSE
				this.setRadius(1.850);
			}
			else if(atomName.startsWith("H")) { // PARSE 
				this.setRadius(1.000);
			}
			else{
				this.setRadius(1.500);
			}
		}
		
		if(this.isAromatic() && this.getName().trim().startsWith("C")){
			this.setRadius(1.700);
		}
		// The vdW radius of aliphatic hydrogens is included in the vdW radius of carbon atoms 
		if(this.covalentBondAtom!=null){
			String atomName = this.getName().trim();
			if(atomName.startsWith("H") && 
			   this.covalentBondAtom.getName().trim().startsWith("C") &&
			   !this.covalentBondAtom.isAromatic) {  
				this.setRadius(0.000);
			}
		}
	}

	//-----------------------------------------------------------------------------------------------------------//
	
    /**
     * Set the flag of an atom for a PDB file
     * Gives a warning if not "ATOM" or "HETATM" and an error message 
     * if the flag is longer than six characters, defaults to ATOM.
     * @param text String
     * @return - none
     */ 
    public void setFlag (String text) {
		if (text.length() > 6) {
		    System.err.println ("WARNING PDB flag \""+text+"\" > six characters!");
	    	flag = "ATOM  ";
	    	return;
		}
		if (!(text.equalsIgnoreCase("ATOM  ")) && !(text.equals("HETATM"))) {
		    System.err.println ("Warning PDB flag \""+text+"\" neither ATOM nor HETATM!");
		}
		flag = text;
	return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the PDB atom serial number
     * @param n int 
     * @return - none
     */
    public void setSerial (int n) {
		if ( (n > 99999) || (n < -9999) ) {
		    System.err.println ("WARNING PDB serial \""+n+"\" number out of range!");
	    	return;
		}
		serial = n;
	return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the aromaticity flag to true
     * @param - none
     * @return - none
     */ 
    public void setAromatic () {
    	isAromatic = true;
    }
	/*--------------------------------------------------------------------------*/
    /**
     * returns true if aromaticity flag was set previously
     * @param - none
     * @return - boolean
     */ 
    public boolean isAromatic() {
    	return isAromatic;
    }
	/*--------------------------------------------------------------------------*/

    /**
     * Set the name of an atom for a PDB file
     * @param text String
     * @return - none
     */ 
    public void setName (String text) {
		if (text.length() > 4) {
		    System.err.println ("WARNING PDB name > \""+text+"\" four characters!");
	    	System.err.println (text.length() );
		    return;
		}
		name = text;
	return;
    }

	/*--------------------------------------------------------------------------*/

	public void setRemark (String text) {
		remark = text;
	}

	/*--------------------------------------------------------------------------*/
    /**
     * Set the proteinAltLoc of an atom for a PDB file
     * @param text char
     * @return - none
     */ 
    public void setAltLoc (char text) {
		altLoc = text;
		return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the resName of an atom for a PDB file
     * @param text String
     * @return - none
     */ 
    public void setResName (String text) {
	if (text.length() > 3) {
	    System.err.println ("WARNING PDB resName >= \""+text+"\" four characters!");
	    return;
	}
		resName = text;
	return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the chainID of an atom for a PDB file
     * @param text char
     * @return - none
     */ 
    public void setChainId (char text) {
		chainID = text;
		return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the resSeq of an atom for a PDB file
     * @param n int
     * @return - none
     */ 
    public void setResNo ( int n) {
		if ( (n > 9999) || (n < -999) ) {
	    	System.err.println ("WARNING PDB resSeq number \""+n+"\" out of range!");
	    	return;
		}
		resNo = n;
	return;
    };    

	/*--------------------------------------------------------------------------*/

    /**
     * Set the icode of an atom for a PDB file
     * @param text char
     * @return - none
     */ 
    public void setICode (char text) {
		iCode = text;
		return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the atom's x coordinate
     * @param a double
     * @return - none
     */
    public void setX (double a) {
		if ((a > 99999.999) || (a < -9999.999)) {
		    System.err.println ("WARNING PDB x out of bounds!: "+a);
	    	return;
		}
		this.x = a;
		
	return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * Set the atom's y coordinate
     * @param a double
     * @return - none
     */
    public void setY (double a) {
		if ((a > 99999.999) || (a < -9999.999)) {
	    	System.err.println ("WARNING PDB y out of bounds!: "+a);
	    	return;
		}
		this.y = a;
	return;
    }
    
	/*--------------------------------------------------------------------------*/

    /**
     * Set the atom's z coordinate
     * @param a double
     * @return - none
     */
    public void setZ (double a) {
		if ((a > 99999.999) || (a < -9999.999)) {
	    	System.err.println ("WARNING PDB z out of bounds!: "+a);
	    	return;
		}
		this.z = a;
	return;
    }

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the radius for spherical coordinate
	 * @param double r - radius 
	 * @return - none
	 */
	public void setR (double radius) {
		
		this.r = radius;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the phi degree for spherical coordinate
	 * @param double phiDegree - phi degree 
	 * @return - none
	 */
	public void setPhi (double phiDegree) {		
		phi = phiDegree;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the theta degree for spherical coordinate
	 * @param double thetaDegree - theta degree 
	 * @return - none
	 */
	public void setTheta (double thetaDegree) {
		
		theta = thetaDegree;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * sets the weight of the atom
	 * @param double weight - weight of the atom
	 * @return - none
	 */
	public void setWeight (double weight) {
		
		this.weight = weight;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * sets the radius of the atom
	 * @param double radius - radius of the atom
	 * @return - none
	 */
	public void setRadius (double radius) {
		
		this.radius = radius;
	}
	/*--------------------------------------------------------------------------*/

	/**
	 * increases the radius of the atom
	 * @param double increase - value to increase radius
	 * @return - none
	 */
	public void increaseRadiusBy (double increase) {
		
		this.radius += increase;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * sets the radius of the atom
	 * @param - none
	 * @return - none
	 */
	
	public void setRadius () {
		
		if(ff.equalsIgnoreCase("charmm")){
			this.setCHARMMradius();
		}
		else if(ff.equalsIgnoreCase("surfnet")){
			this.setVdwPlusHydrogenRadius();
		}
		else if(ff.equalsIgnoreCase("mmff")){
			this.setMMFFradius();
//			this.setPARSEradius();
		}
		
		else{
			this.setPARSEradius();
		}
	}
	
	
	/*--------------------------------------------------------------------------*/

	/**
	 * returns the ff type for atomic radii
	 * @param - none
	 * @return - String: parse, charmm, amber
	 */
	static public String getFF () {
		return ff;
	}
	
	/*--------------------------------------------------------------------------*/

	/**
	 * sets the ff type for atomic radii
	 * @param - String: parse, charmm, mmff, surfnet
	 * @return - none
	 */
	static public void setFF (String parseMmffSurfnetCharmm) {
		ff = parseMmffSurfnetCharmm;
	}
	
	/*--------------------------------------------------------------------------*/

	/**
	 * sets the radius of the atom to the SURFNET specified radius
	 * @param - none
	 * @return - none
	 */
	public void setSurfnetRadius () {
		this.setVdwPlusHydrogenRadius();
	}
	/*--------------------------------------------------------------------------*/

	/**
	 * sets the atomic solvation parameter for hydrophobicity calculation
	 * @param double solvParam - solvation parameter of atom
	 * @return - none
	 */
	public void setSolvParam (double solvParam) {
		this.solvationParameter = solvParam;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * returns the atomic solvation parameter for hydrophobicity calculation
	 * @param -none
	 * @return - double - solvation parameter of atom
	 */
	public double getSolvParam () {
		return this.solvationParameter;
	}
	/*--------------------------------------------------------------------------*

	/**
	 * sets the atomic logP for hydrophobicity calculation
	 * @param double atomicLogP - atomic logP
	 * @return - none
	 */
//	public void setLogP (double atomicLogP) {
//		this.atomicLogP = atomicLogP;
//	}

	/*--------------------------------------------------------------------------*

	/**
	 * returns the atomic logP for hydrophobicity calculation
	 * @param -none
	 * @return - double atomicLogP - atomic logP
	 */
//	public double getLogP () {
//		return this.atomicLogP;
//	}
	/*--------------------------------------------------------------------------*/
	/**
	 * sets the colour grade of the atom
	 * @param dint colour - colour grade of the atom
	 * @return - none
	 */
	public void setColour (int colourGrade) {
		colour = colourGrade;
	}
	/*--------------------------------------------------------------------------*/
	/**
	 * gets the colour grade of the atom
	 * @param -none
	 * @return double colour - colour grade of the atom
	 */
	public int getColour () {
		return colour;
	}

	/*--------------------------------------------------------------------------*/
	/**
	 * gets the weight of the atom
	 * @param -none
	 * @return double weight - weight of the atom
	 */
	public double getWeight () {
		return weight;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the radius of atom
	 * @param - none
	 * @return double radius - radius 
	 */
	public double getRadius () {
		return radius;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the radius of spherical coordinate
	 * @param - none
	 * @return double r - radius 
	 */
	public double getR () {
		return r;
	}

	/*--------------------------------------------------------------------------*/
	/**
	 * get the phi degree of spherical coordinate
	 * @param - none
	 * @return double phiDegree - phi degree
	 */
	public double getPhi () {
		return phi;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the theta degree of spherical coordinate
	 * @param - none 
	 * @return double thetaDegree - theta degree
	 */
	public double getTheta () {
		return theta;
	}


	/*--------------------------------------------------------------------------*/

    /**
     * Set the atom's occupancy 
     * @param a double
     * @return - none
     */
    public void setOccupancy (double a) {
		if (a > 999.99){
			occupancy = 999.99;
			System.err.println ("WARNING PDB occupancy > 999.99 ("+a+"). Sets to 999.99!");
	    	return;
		}
		else if(a < -99.99){
			occupancy = -99.99;
			System.err.println ("WARNING PDB occupancy > -99.99 ("+a+"). Sets to 99.99!");
	    	return;			
		}
		else{
			occupancy = a;
		}
	return;
    }
    
	/*--------------------------------------------------------------------------*/

    /**
     * Set the atom's occupancy 
     * @param a double
     * @return - none
     */
    public void setTemp (double temp) {
		if (temp > 999.99) {
			tempFactor = 999.99;
			System.err.println ("WARNING PDB tempFactor > 999.99 ("+temp+"). Sets to 999.99 for atom: "+this);
			return;
		}
		else if (temp < -99.99){
			tempFactor = 99.99;
			System.err.println ("WARNING PDB tempFactor < -99.99 ("+temp+"). Sets to -99.99 for atom: "+this);
			return;
		}
		else{
			tempFactor = temp;
		}
	return;
    }
	/*--------------------------------------------------------------------------*/

    /**
     * Set the atom's occupancy 
     * @param a double
     * @return - none
     */
    public void setQuickTemp (double temp) {
		if (temp > 999.99) {
			tempFactor = 999.99;
			return;
		}
		else if (temp < -99.99){
			tempFactor = 99.99;
			return;
		}
		else{
			tempFactor = temp;
		}
	return;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * get the temperature factor of an atom for a PDB file
     * @param - none
     * @return - double tempFactor = temperature factor of an atom as a String
     */ 
    public double getTemp () {
		return tempFactor;
    }

	/*--------------------------------------------------------------------------*/

	/**
	 * get the segID of an atom for a PDB file
	 * @param - none
	 * @return - segment ID of an atom as a String
	 */ 
	public String getSegID () {
		return segID;
	}


	/*--------------------------------------------------------------------------*/

   	/**
     * get the element 
     * @param - none 
     * @return - element as a String
     */ 
    public String getElement () {
		return element;
    }

	/*--------------------------------------------------------------------------*/

    /**
     * get the potential
     * @param - none
     * @return - potential as an String
     */ 
    public double getPotential () {
		return potential;
    }
    
	/*--------------------------------------------------------------------------*/

    /**
     * get the anisotropic temperature factors
     * @param - none 
     * @return - none
     */
    public int[] getAniso () {
		int[] aniso = {u00, u11, u22, u01, u02, u12};
		return aniso;
    }

	/*--------------------------------------------------------------------------*/

	/**
	 * get the flag of an atom for a PDB file
	 * @param - none
	 * @return - flag as a String
	 */ 
	public String getFlag () {
		return flag;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the PDB atom serial number
	 * @param - none
	 * @return - serial number as a int
	 */
	public int getSerial () {
		return serial;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the name of an atom for a PDB file
	 * @param -none 
	 * @return - name of atom as a String
	 */ 
	public String getName () {
		return name;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the proteinAltLoc of an atom for a PDB file
	 * @param - none
	 * @return - proteinAltLoc of an atom as a String
	 */ 
	public char getAltLoc () {
		return altLoc;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the resName of an atom for a PDB file
	 * @param - none 
	 * @return - residue name of a atom as a String
	 */ 
	public String getResName () {
		return resName;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the chainID of an atom for a PDB file
	 * @param - none
	 * @return - chain ID of an atom as a String
	 */ 
	public char getChainId () {
		return chainID;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the resSeq of an atom for a PDB file
	 * @param - none
	 * @return - residue number of a sequence as a int
	 */ 
	public int getResNo () {
		return resNo;
	};    

	/*--------------------------------------------------------------------------*/

	/**
	 * get the iCode of an atom for a PDB file
	 * @param - none 
	 * @return - get iCode of an atom as a String
	 */ 
	public char getICode () {
		return iCode;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the atom's x coordinate
	 * @param - none 
	 * @return - x coordinate as a double
	 */
	public double getX () {
		return x;
	}
	
	/*--------------------------------------------------------------------------*/

	/**
	 * get the atom's y coordinate
	 * @param - none 
	 * @return - y coordinate as a double
	 */
	public double getY () {
		return y;
	}
    
	/*--------------------------------------------------------------------------*/

	/**
	 * get the atom's z coordinate
	 * @param - none
	 * @return - z coordinate as a double
	 */
	public double getZ () {
		return z;
	}

	/*--------------------------------------------------------------------------*/
	public String getRemark () {
		return remark;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * get the atom's occupancy 
	 * @param - none
	 * @return - occupancy of an atom as a double
	 */
	public double getOccupancy () {
		return occupancy;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * get position in cluster 
	 * @param - none
	 * @return - cluster position as integer number
	 */
	public int getClusterPos () {
		return clusterPos;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * get atom type according to MMFF94 
	 * @param - none
	 * @return - mmff94 - atom type
	 */
	public String getMMFFatomType () {
		return mmffAtomType;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * set type of atom accoring to MMFF94 
	 * @param - mmff94 - atom type
	 * @return - none
	 */
	public void setMMFFatomType (String mmffAtomType) {
		this.mmffAtomType = mmffAtomType;
	}
	/*--------------------------------------------------------------------------*/

	/**
	 * get atom type of eccp/3 force field 
	 * @param - none
	 * @return - eccp/3 - atom type
	 */
//	public String getEccpAtomType () {
//		return eccpAtomType;
//	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * set atom type of eccp/3 force field 
	 * @param - eccp/3 - atom type
	 * @return - none
	 */
//	public void setEccpAtomType (String eccpAtomType) {
//		this.eccpAtomType = eccpAtomType;
//	}
	/*--------------------------------------------------------------------------*/

	/**
	 * set atom type of eccp/3 force field 
	 * @param - none
	 * @return - none
	 */
/*	public void setEccpAtomType () {
		switch(this.getName().trim().charAt(0)){
		case('C'):{
			this.setEccpAtomType("c6");
			break;
		}
		case('O'):{
			this.setEccpAtomType("o18");			
			break;
		}
		case('N'):{
			this.setEccpAtomType("n14");
			break;
		}
		case('H'):{
			this.setEccpAtomType("h1");
			break;
		}
		case('s'):{
			this.setEccpAtomType("s20");
			break;
		}
		default:{
			this.setEccpAtomType("c8");
		}
		}
	}
  */      
	/*--------------------------------------------------------------------------*/

	/**
	 * get partial charge 
	 * @param - none
	 * @return - partial charge of atom
	 */
	public double getPartialCharge () {
		return partialCharge;
	}
	/*--------------------------------------------------------------------------*/

	/**
	 * set partial charge 
	 * @param - partial charge of atom
	 * @return - none
	 */
	public void setPartialCharge (double partialCharge) {
		this.partialCharge = partialCharge;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * get atomic logP  
	 * @param - none
	 * @return - atomic logP value
	 */
	public double getXlogP() {
		return xlogP;
	}
	/*--------------------------------------------------------------------------*/

	/**
	 * set atomic logP 
	 * @param - atomic logP value
	 * @return - none
	 */
	public void setXlogP (double logP) {
		this.xlogP = logP;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * get hybridization state 
	 * @param - none
	 * @return - hybridization state as string: 1,2,3
	 */
	public int getHybridizationState () {
		return hybridizationState;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * set hybridization state 
	 * @param - hybridization state as string: sp,sp2,sp3
	 * @return - none
	 */
	public void setHybridizationState (int hybridizationState) {
		this.hybridizationState = hybridizationState;
	}
        
	/*--------------------------------------------------------------------------*/

	/**
	 * Set the segID of an atom for a PDB file
	 * @param text String
	 * @return - none
	 */ 
	public void setSegID (String text) {
		if (text.length() > 4) {
			System.err.println ("WARNING PDB segID > \""+text+"\" four characters!");
			return;
		}
		segID = text;
	return;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the element 
	 * @param text String
	 * @return - none
	 */ 
	public void setElement (String text) {
		if (text.length() > 2) {
			System.err.println ("WARNING PDB element > \""+text+"\" two characters!");
			return;
		}
		element = text;
	return;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the potential
	 * @param text int
	 * @return - none
	 */ 
	public void setPotential (double n) {
		potential = n;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the anisotropic temperature factors
	 * @param aniso00 int 
	 * @param aniso11 int  
	 * @param aniso22 int  
	 * @param aniso01 int  
	 * @param aniso02 int  
	 * @param aniso12 int 
	 * @return - none
	 */
	public void setAniso (int aniso00, int aniso11,int aniso22, int aniso01, int aniso02, 
						  int aniso12) {
		u00 = aniso00;      
		u11 = aniso11;      
		u22 = aniso22;      
		u01 = aniso01;      
		u02 = aniso02;      
		u12 = aniso12;
		return;
	}

	/*--------------------------------------------------------------------------*/

	/**
	 * Set the position in a cluster
	 * @param clusterPosition - int 
	 * @return - none
	 */
	public void setClusterPos (int clusterPosition) {
		clusterPos = clusterPosition;
	}

	/*--------------------------------------------------------------------------*/

    /**
     * Rather ugly method of formatting the output with String functions 
     * to standard PDB format (but it seems to work alright)
     * @param - none
     * @return - output string
     */
    public String toString(double occupancy, double tempFactor) {

		Locale.setDefault(Locale.US);
		NumberFormat decFormat = new DecimalFormat("0.000");

		StringBuffer output = new StringBuffer();

//		try{
		for (int j=this.getFlag().length();j<6; j++) {
			output.append(" ");
		}
		output.append(this.getFlag());
		
		for (int j=Integer.toString(this.getSerial()).length();j<5; j++) {
			output.append(" ");
		}
		output.append(this.getSerial());

		output.append(" ");
					
		for (int j=this.getName().length();j<4; j++) {
			output.append(" ");
		}
		output.append(this.getName());
	
		output.append(this.getAltLoc());
		
		for (int j=this.getResName().length();j<3; j++) {
			output.append(" ");
		}
		output.append(this.getResName());
		
		output.append(" ");

		output.append(this.getChainId());

		for (int j=Integer.toString(this.getResNo()).length();j<4; j++) {
			output.append(" ");
		}
		output.append(this.getResNo());
		
		output.append(this.getICode());
		
		output.append("   ");

		for (int j=decFormat.format(this.getX()).length();j<8; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getX()));
			
			
		for (int j=decFormat.format(this.getY()).length();j<8; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getY()));


		for (int j=decFormat.format(this.getZ()).length();j<8; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getZ()));

		decFormat = new DecimalFormat("0.00");

		for (int j=decFormat.format(occupancy).length();j<6; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(occupancy));

		for (int j=decFormat.format(tempFactor).length();j<6; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(tempFactor));

		output.append("      ");

		for (int j=this.getSegID().length();j<4; j++) {
			output.append(" ");
		}
		output.append(this.getSegID());

		for (int j=this.getElement().length();j<2; j++) {
			output.append(" ");
		}
		output.append(this.getElement());
//		output.append(" ");
//		output.append(decFormat.format(this.getPotential()));

		output.append("\n");
//	}
//	catch(IOException e){
//		System.err.print("Could not print out Atom line\n\n");
//	}

	return output.toString();
    }
    
	/*--------------------------------------------------------------------------*/

    /**
     * Rather ugly method of formatting the output with String functions 
     * to standard PDB format (but it seems to work alright)
     * @param - none
     * @return - output string
     */
    public String toString() {

		Locale.setDefault(Locale.US);
		NumberFormat decFormat = new DecimalFormat("0.000");

		StringBuffer output = new StringBuffer();

//		try{
		for (int j=this.getFlag().length();j<6; j++) {
			output.append(" ");
		}
		output.append(this.getFlag());
		
		for (int j=Integer.toString(this.getSerial()).length();j<5; j++) {
			output.append(" ");
		}
		output.append(this.getSerial());

		output.append(" ");
					
		for (int j=this.getName().length();j<4; j++) {
			output.append(" ");
		}
		output.append(this.getName());
	
		output.append(this.getAltLoc());
		
		for (int j=this.getResName().length();j<3; j++) {
			output.append(" ");
		}
		output.append(this.getResName());
		
		output.append(" ");

		output.append(this.getChainId());

		for (int j=Integer.toString(this.getResNo()).length();j<4; j++) {
			output.append(" ");
		}
		output.append(this.getResNo());
		
		output.append(this.getICode());
		
		output.append("   ");

		for (int j=decFormat.format(this.getX()).length();j<8; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getX()));
			
			
		for (int j=decFormat.format(this.getY()).length();j<8; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getY()));


		for (int j=decFormat.format(this.getZ()).length();j<8; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getZ()));

		decFormat = new DecimalFormat("0.00");

		for (int j=decFormat.format(this.getOccupancy()).length();j<6; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getOccupancy()));

		for (int j=decFormat.format(this.getTemp()).length();j<6; j++) {
			output.append(" ");
		}
		output.append(decFormat.format(this.getTemp()));

		output.append("      ");

		for (int j=this.getSegID().length();j<4; j++) {
			output.append(" ");
		}
		output.append(this.getSegID());

		for (int j=this.getElement().length();j<2; j++) {
			output.append(" ");
		}
		output.append(this.getElement());
		output.append(" ");		
//		output.append(decFormat.format(this.getPotential()));
		output.append(System.getProperty("line.separator"));
//	}
//	catch(IOException e){
//		System.err.print("Could not print out Atom line\n\n");
//	}

	return output.toString();
    }
	/*--------------------------------------------------------------------------*/

	public boolean equalsAtom(Object o){
		StructureAtom atom = (StructureAtom)o;
		if((this.getName().trim().equals(atom.getName().trim())) && 
		   (this.getChainId() == atom.getChainId()) &&
		   (this.getAltLoc() == atom.getAltLoc()) && 
		   (this.getResName().trim().equals(atom.getResName().trim())) &&
		   this.getResNo() == atom.getResNo() &&
		   Math.abs(this.getX()-atom.getX())<0.0001 &&
		   Math.abs(this.getY()-atom.getY())<0.0001 &&
		   Math.abs(this.getZ()-atom.getZ())<0.0001){
		   	return true;
		   }
		else {
			return false;
		}
	}

	/*--------------------------------------------------------------------------*/

	public boolean isEqualType(StructureAtom atom){
		if((this.getFlag().equals(atom.getFlag())) &&
		   (this.getName().trim().equals(atom.getName().trim())) &&
		   (this.getResName().trim().equals(atom.getResName().trim())) ){
		   	return true;
		}	
		else {
			return false;
		}
	}
	/*--------------------------------------------------------------------------*/

	public boolean isEqualResidue(StructureAtom atom){
		if( (this.getResName().trim().equals(atom.getResName().trim()) ) &&
		    (this.getResNo()   == atom.getResNo()) &&
		    (this.getChainId() == atom.getChainId()) ){
		   	return true;
		}	
		else {
			return false;
		}
	}
	/*--------------------------------------------------------------------------*/

	public boolean isEqualPos(StructureAtom atom){
		if( Math.abs(this.getX()-atom.getX())<=0.001 &&
			Math.abs(this.getY()-atom.getY())<=0.001 &&
			Math.abs(this.getZ()-atom.getZ())<=0.001 ){
		   	return true;
		}	
		else {
			return false;
		}
	}
	/*--------------------------------------------------------------------------*/

	public double distance(StructureAtom atom){
		
		return Mathematics.distance(atom, this) - this.getRadius() - atom.getRadius();
	}

	/*--------------------------------------------------------------------------*/

	public void move(double[] move){
		this.setX(this.getX()+move[0]);
		this.setY(this.getY()+move[1]);
		this.setZ(this.getZ()+move[2]);
	}

	//----------------------------------------------------------------------------//
	public StructureCoords getCovalentBondAtoms(StructureCoords coords, boolean withoutHydrogen){

		// Covalent bonds in organic molecules should have a distance less than 1.54  between both atom centres. (see http://en.wikipedia.org/wiki/Bond_length)
		// To account for resolution error distance is taken to be less than 1.8 (0.28  estimated standard error for X-ray structures)
		StructureCoords covalentAtoms = new StructureCoords();
		for(Iterator i=coords.iterator(); i.hasNext();){			
			StructureAtom atom2 = (StructureAtom)i.next();
			
			double maxDist = (StructureAtom.getVdWradius(this.getElement())+StructureAtom.getVdWradius(atom2.getElement()));
			maxDist = (maxDist/2)+this.covalentRadiusErrorRange;

			if(withoutHydrogen){
				if(atom2.getName().trim().startsWith("H")) continue;
			}
			double dist = Mathematics.distance(this, atom2);
			if(dist <=maxDist && !this.equalsAtom(atom2) && !atom2.isMetal(coords)){
				covalentAtoms.add(atom2);
			}
		}
/*		if(covalentAtoms.size()==0 && increaseMaxDist){
			double steps = 0.05;
			for(double dist=maxDist+steps;dist<2.0;dist+=steps){
				for(Iterator i=coords.iterator(); i.hasNext();){			
					StructureAtom atom2 = (StructureAtom)i.next();
					if(withoutHydrogen){
						if(atom2.getName().trim().startsWith("H")) continue;
					}
					if(Mathematics.distance(this, atom2)<=dist && !this.equals(atom2) && !atom2.isMetal(coords)){
						covalentAtoms.add(atom2);
					}
				}

				if(covalentAtoms.size()>0){
					//					System.err.print("WARNING: Couldn't find covalent atoms for dist \""+maxDist+"\". Had to increase distance to \""+dist+"\"\n"+this+"\n");					
					break;
				}
			}
		}
		*/
	return covalentAtoms;
	}
	//----------------------------------------------------------------------------//
	public StructureCoords getCovalentBondHydrogens(StructureCoords coords){
		StructureCoords hydrogens = new StructureCoords();
		for(Iterator i=coords.iterator(); i.hasNext();){			
			StructureAtom atom = (StructureAtom)i.next();
			if(!atom.getElement().equals("H")) continue;
			if(Mathematics.distance(this, atom)<=1.2 && !this.equalsAtom(atom)){
				hydrogens.add(atom);
			}
		}
	return hydrogens;
	}
	//----------------------------------------------------------------------------//
	
	static public double getVdWradius (String atomName){
		StructureAtom atom = new StructureAtom();
		atom.setName(atomName);
		atom.setRadius();
	return atom.getRadius();
	}
	
	//----------------------------------------------------------------------------//
	
	static public double getCovalentRadius (StructureAtom atom1, StructureAtom atom2){
		double maxDist = StructureAtom.getVdWradius(atom1.getElement())+StructureAtom.getVdWradius(atom2.getElement());
		return (maxDist/2)+StructureAtom.covalentRadiusErrorRange;
	}

	//----------------------------------------------------------------------------//
	// check whether coords are HETATM and then their atom number and then their name to assess whether atom is metal or not 
	public boolean isMetal(StructureCoords coords){
		
		if(this.getFlag().equals("ATOM  ")) 
			return false;

		if(coords.size()==1){
			String atomName = this.getName().trim();		
			if(atomName.equals("CA")){ 
				for(Iterator i=coords.iterator(); i.hasNext();){
					StructureAtom atom = (StructureAtom)i.next();
					double dist = this.distance(atom);
					if(dist<=0){
						return false;
					}
				}
				return true;
			}
			else if(atomName.equals("LI")){ 
				return true;
			}
			else if(atomName.equals("BE")){ 
				return true;
			}
			else if(atomName.equals("NA")){ 
				return true;
			}
			else if(atomName.equals("MG")){ 
				return true;
			}
			else if(atomName.equals("K")){ 
				return true;
			}
			else if(atomName.equals("AL")){ 
				return true;
			}
			else if(atomName.equals("SE")){ 
				return true;
			}
			else if(atomName.equals("SI")){ 
				return true;
			}
			else if(atomName.equals("FE")){ 
				return true;
			}
			else if(atomName.equals("ZN")){ 
				return true;
			}
			else if(atomName.equals("MN")){ 
				return true;
			}
			else if(atomName.equals("LU")){ 
				return true;
			}
			else if(atomName.equals("CD")){ 
				return true;
			}
			else if(atomName.equals("I")){ 
				return true;
			}
			else if(atomName.equals("U")){ 
				return true;
			}
			else if(atomName.equals("HG")){ 
				return true;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	}

	//----------------------------------------------------------------------------//
	// check whether atom is C, O, H, N 
	public boolean isOrganic(StructureCoords coords){
		
		if(this.isMetal(coords))
			return false;
		
		// if molecule consists of only one atom, it is most likely a metal ion
		String atomName = this.getName().trim();		
		if(coords.size()==1){
			return false;
		}		
		else if(atomName.startsWith("C")){ 
			return true;
		}
		else if(atomName.startsWith("O")){ 
			return true;
		}
		else if(atomName.startsWith("N")){ 
			return true;
		}
		else if(atomName.startsWith("H")){ 
			return true;
		}
		else if(atomName.startsWith("S")){ 
			return true;
		}
		else if(atomName.startsWith("P")){ 
			return true;
		}
		else{
			return false;
		}
	}

	//----------------------------------------------------------------------------//
	// check whether atom is H 
	public boolean isHydrogen(){
		
		// if molecule consists of only one atom, it is most likely a metal ion
		String atomName = this.getName().trim();
		if(atomName.startsWith("H")){ 
			return true;
		}
		else if(atomName.matches("\\d+H.*")){ 
			return true;
		}
		return false; 
	}	
	
} // End of class StructureAtom








