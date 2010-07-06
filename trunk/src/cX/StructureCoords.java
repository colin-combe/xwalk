package cX;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Collections;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.util.Locale;


/**
 * Class with functions for handling input and output of PDB format files
 */
public class StructureCoords extends ArrayList {

	/*--------------------------------------------------------------------------*/
	// Classobjects
	/*--------------------------------------------------------------------------*/    

	protected String volume = "";
	protected String fileName = "";
	protected String remark = "";
	
	protected double gridLen = 0.5;
	
	protected double[] originalCentre1 = new double[3];
	protected double[] originalCentre2 = new double[3];

	/*--------------------------------------------------------------------------*/
	// Constructors
	/*--------------------------------------------------------------------------*/    

	/**
	 * Constructor - model is read from given file name  
	 * @param - none
	 * @return - none
	 */
	public StructureCoords () {
	}

	/**
	 * Constructor - model is read from given file name  
	 * @param fileName - path of file to read
	 * @return - none
	 */
	public StructureCoords (String fileName) {
		this.fileName = fileName;
	}

	/*--------------------------------------------------------------------------*/
	// Classmethods
	/*--------------------------------------------------------------------------*/    

    /**
     * Read the model from a given PDB file
     * @param fileName String
     */ 
	public void extractAtoms (String flags, String chainIds, String alternativeLocs, 
							  boolean getH2O, boolean getHydrogen) {		
   		try {

			InputStream ins = this.getClass().getClassLoader().getResourceAsStream(fileName);
   			BufferedReader f;
   			if(ins==null){
   		   		f = new BufferedReader (new FileReader (fileName));
   			}
   			else{
   				f = new BufferedReader(new InputStreamReader(ins));
   			}

	   		String s;
	   		String flag;
	   		String name;
	   		char chain;
	   		char altLoc;
	   		String resName;
	   		String resNo;

	   		// PDB2PQR can't handle ICode columns. So those atoms that have an ICode, those have to get assigned a new
	   		// ResNo
	   		StructureCoords v = new StructureCoords();
	   		int maxResNo = 0;
			while ((s = f.readLine()) != null) {
		   		if (s.length() > 26) {
		   			flag 	= s.substring(0, 6); 
		   			name 	= s.substring(12, 16); 
		   			altLoc 	= s.substring(16, 17).charAt(0);
		   			chain	= s.substring(21, 22).charAt(0);
		   			resName = s.substring(17, 20);

		   			if( (flags.indexOf(flag) != -1) &&
		   				(chainIds.indexOf(chain) != -1) && 
		   			  	(alternativeLocs.indexOf(altLoc) != -1) && 
		   			  	( (!resName.equals("HOH") || getH2O) ) &&
		   			  	( (!name.trim().startsWith("H") || getHydrogen ) )
		   			) 
		   			{
		   				StructureAtom atom = this.parseAtom(s);

			   			if(atom.getICode()==' '){
			   				this.add(atom);			   				
			   			}
		   				maxResNo = Math.max(maxResNo, atom.getResNo());
		   			}
		   		}

				// NMR structures have a lot of Modelstructure in one PDB file. It should
				// be enough to read just the first model.
//				if (s.startsWith("ENDMDL") && (this.size() > 10)){
//					break;
//				}
	   		}
			if(v.size()>0){
				StructureResidue[] resi = v.getResidues();
				for(int i=0; i<resi.length; i++){
					for(Iterator j=resi[i].iterator(); j.hasNext();){
						StructureAtom atom = (StructureAtom) j.next();
						atom.setResNo(maxResNo+i);
						atom.setICode(' ');
						System.err.print("WARNING: Atom has ICode which PDB2PQR can't handle, " +
										 "so delete ICode and assign new residue number\n"+atom);
					}
				}
			}
       }
       catch (IOException e){
	       System.err.println ("Problems reading model from \""+ fileName+"\" "+e.getMessage());
	   }
	   // before parsing the coordinates, the ligand object has to initialized. This could not
	   // be done in the constructures because the new StructureCoords() would call itself recursivly
	   // till a StackOverFlow Error massage
	   
		// transform the PDB coordinates to the coordinate origin, for visualisation purposes.
		// All shapes showing in a visualisation program are then centered and not spread
		// over whole 3D space
       
       // assign element names. Has to be done here to check properly for atom being a metal
       StructureCoords atomCoords = new StructureCoords();
       StructureCoords hetatmCoords = new StructureCoords();       

       for(Iterator k=this.iterator(); k.hasNext();){
    	   StructureAtom atom = (StructureAtom)k.next();
    	   if(atom.getFlag().equals("ATOM  "))
    		   atomCoords.add((Object)atom);
    	   else if(atom.getFlag().equals("HETATM"))
    		   hetatmCoords.add((Object)atom);
       }
       
       if(atomCoords.size()>0){
           for(Iterator j=atomCoords.iterator(); j.hasNext();){
        	   StructureAtom atom = (StructureAtom)j.next();
        	   // assign element
        	   String name = atom.getName().trim();
        	   String nameWithoutDigits = name.replaceAll("[0-9]+","");
        	   atom.setElement(String.valueOf(nameWithoutDigits.charAt(0)));
           }    	   
       }
	}
	/*--------------------------------------------------------------------------*/    

    /**
     * Read the model from a given PDB file
     * @param fileName String
     */ 
	public void extractAtoms (String flags, String chainIds, String alternativeLocs, 
							  boolean getH2O, boolean getHydrogen, boolean verbose) {		
   		try {

			InputStream ins = this.getClass().getClassLoader().getResourceAsStream(fileName);
   			BufferedReader f;
   			if(ins==null){
   		   		f = new BufferedReader (new FileReader (fileName));
   			}
   			else{
   				f = new BufferedReader(new InputStreamReader(ins));
   			}

	   		String s;
	   		String flag;
	   		String name;
	   		char chain;
	   		char altLoc;
	   		String resName;
	   		String resNo;

	   		// PDB2PQR can't handle ICode columns. So those atoms that have an ICode, those have to get assigned a new
	   		// ResNo
	   		StructureCoords v = new StructureCoords();
	   		int maxResNo = 0;
			while ((s = f.readLine()) != null) {
		   		if (s.length() > 26) {
		   			flag 	= s.substring(0, 6); 
		   			name 	= s.substring(12, 16); 
		   			altLoc 	= s.substring(16, 17).charAt(0);
		   			chain	= s.substring(21, 22).charAt(0);
		   			resName = s.substring(17, 20);

		   			if( (flags.indexOf(flag) != -1) &&
		   				(chainIds.indexOf(chain) != -1) && 
		   			  	(alternativeLocs.indexOf(altLoc) != -1) && 
		   			  	( (!resName.equals("HOH") || getH2O) ) &&
		   			  	( (!name.trim().startsWith("H") || getHydrogen ) )
		   			) 
		   			{
		   				StructureAtom atom = this.parseAtom(s, verbose ? 1:0);
		   				// Reject any inserted amino acids
//		   				if(atom.getICode()==' '){
			   				this.add(atom);
//			   			}

		   				maxResNo = Math.max(maxResNo, atom.getResNo());

		   			}
		   		}

				// NMR structures have a lot of Modelstructure in one PDB file. It should
				// be enough to read just the first model.
//				if (s.startsWith("ENDMDL") && (this.size() > 10)){
//					break;
//				}
	   		}
			if(v.size()>0){
				StructureResidue[] resi = v.getResidues();
				for(int i=0; i<resi.length; i++){
					for(Iterator j=resi[i].iterator(); j.hasNext();){
						StructureAtom atom = (StructureAtom) j.next();
						atom.setResNo(maxResNo+i);
						atom.setICode(' ');
						if(verbose)
							System.err.print("WARNING: Atom has ICode which PDB2PQR can't handle, " +
											 "so delete ICode and assign new residue number\n"+atom);
					}
				}
			}
       }
       catch (IOException e){
    	   if(verbose)
    		   System.err.println ("Problems reading model from \""+ fileName+"\" "+e.getMessage());
	   }
	   // before parsing the coordinates, the ligand object has to initialized. This could not
	   // be done in the constructures because the new StructureCoords() would call itself recursivly
	   // till a StackOverFlow Error massage
	   
		// transform the PDB coordinates to the coordinate origin, for visualisation purposes.
		// All shapes showing in a visualisation program are then centered and not spread
		// over whole 3D space
       
       // assign element names. Has to be done here to check properly for atom being a metal
       StructureCoords atomCoords = new StructureCoords();
       StructureCoords hetatmCoords = new StructureCoords();       

       for(Iterator k=this.iterator(); k.hasNext();){
    	   StructureAtom atom = (StructureAtom)k.next();
    	   if(atom.getFlag().equals("ATOM  "))
    		   atomCoords.add((Object)atom);
    	   else if(atom.getFlag().equals("HETATM"))
    		   hetatmCoords.add((Object)atom);
       }
       
       if(atomCoords.size()>0){
           for(Iterator j=atomCoords.iterator(); j.hasNext();){
        	   StructureAtom atom = (StructureAtom)j.next();
        	   // assign element
        	   String name = atom.getName().trim();
        	   String nameWithoutDigits = name.replaceAll("[0-9]+","");
        	   atom.setElement(String.valueOf(nameWithoutDigits.charAt(0)));
           }    	   
       }
	}
	/*--------------------------------------------------------------------------*/		

	/**
	 * Function for extracting the ATOM fields from a string
	 * @param oneline String
	 */
	public StructureAtom parseAtom (String oneline) {

		int length = oneline.length();
		StructureAtom atom = new StructureAtom();

		try {
			
			atom.setFlag(oneline.substring(0,6));
			atom.setSerial(Integer.parseInt (oneline.substring(6,11).trim())); 
			atom.setName(oneline.substring(12,16));
			atom.setAltLoc(oneline.charAt(16));
			atom.setResName(oneline.substring(17,20));
			atom.setChainId(oneline.charAt(21));    
			atom.setResNo(Integer.parseInt(oneline.substring(22,26).trim()));
			atom.setICode(oneline.charAt(26));
			atom.setX(Double.parseDouble(oneline.substring(30,38).trim()));
			atom.setY(Double.parseDouble(oneline.substring(38,46).trim()));
			atom.setZ(Double.parseDouble(oneline.substring(46,54).trim()));
			atom.setOccupancy(Double.parseDouble(oneline.substring(54,60).trim()));
			atom.setTemp(Double.parseDouble(oneline.substring(60,66).trim()));
// no need to parse in further beyond the temp column
/*			if (length >= 75) {
				atom.setSegID(oneline.substring(72,76));
				if (length >= 78) {
					atom.setElement(oneline.substring(76,78));
					if (length >= 80) {
						if(!oneline.substring(78,80).equals("  ")){
							atom.setCharge(Double.parseDouble(oneline.substring(78,80)));
						}
					}
				}
			}
*/
			// for cases where no information at the rightest columns are 
			// provided, just fill with empty space
			
			if(atom.getSegID().equals("")){
				atom.setSegID("    ");
			}

			if(length<80){
				atom.setPartialCharge(0.0);
			}
			// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
			if (atom.getFlag().equals("ATOM  ") && atom.getName().length() > 0){
//				atom.setRadius();
			}
			// if file is a sample point file than assign them radius = 0
			if( atom.getFlag().equals("ATOM  ") && atom.getResName().equals("SAM")&& atom.getResName().equals("SUR")){
					atom.setRadius(0.0);
			}

			if (atom.getFlag().equals("HETATM") && atom.getName().length() > 0){
//				atom.setRadius();
			}			
		}
		catch (Exception e) {
			System.err.println ("ATOM: Does not seem to be PDB Format: "+e+"\n");
		}
	
	return atom;
	}

	/*--------------------------------------------------------------------------*/		

	/**
	 * Function for extracting the ATOM fields from a string
	 * @param oneline String
	 * @param int 0 if no output should be made, any other number will cause output 
	 */
	public StructureAtom parseAtom (String oneline, int output) {

		int length = oneline.length();
		StructureAtom atom = new StructureAtom();

		try {
			
			atom.setFlag(oneline.substring(0,6));
			atom.setSerial(Integer.parseInt (oneline.substring(6,11).trim())); 
			atom.setName(oneline.substring(12,16));
			atom.setAltLoc(oneline.charAt(16));
			atom.setResName(oneline.substring(17,20));
			atom.setChainId(oneline.charAt(21));    
			atom.setResNo(Integer.parseInt(oneline.substring(22,26).trim()));
			atom.setICode(oneline.charAt(26));
			atom.setX(Double.parseDouble(oneline.substring(30,38).trim()));
			atom.setY(Double.parseDouble(oneline.substring(38,46).trim()));
			atom.setZ(Double.parseDouble(oneline.substring(46,54).trim()));
			atom.setOccupancy(Double.parseDouble(oneline.substring(54,60).trim()));
			atom.setTemp(Double.parseDouble(oneline.substring(60,66).trim()));
// no need to parse in further beyond the temp column
/*			if (length >= 75) {
				atom.setSegID(oneline.substring(72,76));
				if (length >= 78) {
					atom.setElement(oneline.substring(76,78));
					if (length >= 80) {
						if(!oneline.substring(78,80).equals("  ")){
							atom.setCharge(Double.parseDouble(oneline.substring(78,80)));
						}
					}
				}
			}
*/
			// for cases where no information at the rightest columns are 
			// provided, just fill with empty space
			
			if(atom.getSegID().equals("")){
				atom.setSegID("    ");
			}

			if(length<80){
				atom.setPartialCharge(0.0);
			}
			// determing the sort of atom (C,O,N,H,S) to assign the radius to the atom
			if (atom.getFlag().equals("ATOM  ") && atom.getName().length() > 0){
//				atom.setRadius();
			}
			// if file is a sample point file than assign them radius = 0
			if( atom.getFlag().equals("ATOM  ") && atom.getResName().equals("SAM")&& atom.getResName().equals("SUR")){
					atom.setRadius(0.0);
			}

			if (atom.getFlag().equals("HETATM") && atom.getName().length() > 0){
//				atom.setRadius();
			}			
		}
		catch (Exception e) {
			if(output!=0)
				System.err.println ("ATOM: Does not seem to be PDB Format: "+e+"\n");
		}
	
	return atom;
	}
	/*--------------------------------------------------------------------------*/

	/**
	 * Function for extracting the ANISOU fields from a string
	 * @param oneline String
	 */
	private void parseANISOU (ArrayList anisou) {
		for(Iterator i=anisou.iterator();i.hasNext();) {
			String oneline = (String)i.next();
			StructureAtom atom = (StructureAtom)this.get(this.size()-1);
			try {
				int a00 = Integer.parseInt(oneline.substring(28,35).trim());
				int a11 = Integer.parseInt(oneline.substring(35,42).trim());
				int a22 = Integer.parseInt(oneline.substring(42,49).trim());
				int a01 = Integer.parseInt(oneline.substring(49,56).trim());
				int a02 = Integer.parseInt(oneline.substring(56,63).trim());
				int a12 = Integer.parseInt(oneline.substring(63,70).trim());
				atom.setAniso(a00,a11,a22,a01,a02,a12);
			}
			catch (Exception e) {
				System.err.println ("ANISOU: Does not seem to be PDB Format");
			}
		}
	}

	/*--------------------------------------------------------------------------*/

	private double[] centerOfGeometry () {
		double[] center = {0,0,0};
		double xMax = -99999999;//this.getAtom(0).getX()+this.getAtom(0).getRadius();
		double xMin = 999999999;//this.getAtom(0).getX()-this.getAtom(0).getRadius();
		double yMax = -99999999;//this.getAtom(0).getY()+this.getAtom(0).getRadius();
		double yMin = 999999999;//this.getAtom(0).getY()-this.getAtom(0).getRadius();
		double zMax = -99999999;//this.getAtom(0).getZ()+this.getAtom(0).getRadius();
		double zMin = 999999999;//this.getAtom(0).getZ()-this.getAtom(0).getRadius();

		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			double x = atom.getX();
			double y = atom.getY();
			double z = atom.getZ();
			double r = atom.getRadius();
		
			xMax = Math.max(xMax, x+r);
			xMin = Math.min(xMin, x-r);
			yMax = Math.max(yMax, y+r);
			yMin = Math.min(yMin, y-r);
			zMax = Math.max(zMax, z+r);
			zMin = Math.min(zMin, z-r);
		}
		center[0] = xMin+((xMax-xMin)/2);
		center[1] = yMin+((yMax-yMin)/2);
		center[2] = zMin+((zMax-zMin)/2);
		
	return center;
	}
	/*--------------------------------------------------------------------------*/

	public double[] dimension() {
		double[] dim = {0,0,0};
		double xMax = -99999999;//this.getAtom(0).getX()+this.getAtom(0).getRadius();
		double xMin = 999999999;//this.getAtom(0).getX()-this.getAtom(0).getRadius();
		double yMax = -99999999;//this.getAtom(0).getY()+this.getAtom(0).getRadius();
		double yMin = 999999999;//this.getAtom(0).getY()-this.getAtom(0).getRadius();
		double zMax = -99999999;//this.getAtom(0).getZ()+this.getAtom(0).getRadius();
		double zMin = 999999999;//this.getAtom(0).getZ()-this.getAtom(0).getRadius();

		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			double x = atom.getX();
			double y = atom.getY();
			double z = atom.getZ();
			double r = atom.getRadius();
		
			xMax = Math.max(xMax, x+r);
			xMin = Math.min(xMin, x-r);
			yMax = Math.max(yMax, y+r);
			yMin = Math.min(yMin, y-r);
			zMax = Math.max(zMax, z+r);
			zMin = Math.min(zMin, z-r);
		}
		dim[0] = xMax-xMin;
		dim[1] = yMax-yMin;
		dim[2] = zMax-zMin;
	return dim;
	}
	/*--------------------------------------------------------------------------*/

	private double[] centerOfGravity () {
		double[] center = {0,0,0};	
		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			center[0] += atom.getWeight() * atom.getX();
			center[1] += atom.getWeight() * atom.getY();
			center[2] += atom.getWeight() * atom.getZ();
		}
		center[0] /= this.size();
		center[1] /= this.size();
		center[2] /= this.size();
	return center;
	}

	/*--------------------------------------------------------------------------*/

	public double[] move2centreOfGeometry() {
		double[] centre = this.centerOfGeometry();
		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			atom.setX(atom.getX()-centre[0]);
			atom.setY(atom.getY()-centre[1]);
			atom.setZ(atom.getZ()-centre[2]);
		}

/*		// test if really in center of mass
		double[] testCenter = this.centerOfGeometry();
	System.out.print(testCenter[0]+" "+testCenter[1]+" "+testCenter[2]+"\n");	
		if (Math.abs(testCenter[0])>0.00001 || Math.abs(testCenter[1])> 0.00001 ||
		    Math.abs(testCenter[0])> 0.00001){
			System.err.print("!!!Error during moving object to its center of geometry!!!. "+
							 "Object could just be moved to:\n"+
							 "x: "+testCenter[0]+"\n"+
							 "y: "+testCenter[1]+"\n"+
							 "z: "+testCenter[2]+"\n");
		}
*/
		this.calcSphereCoord();
	return centre;
	}

	/*--------------------------------------------------------------------------*/
	public double[] move2centreOfGravity() {
		double[] centre = this.centerOfGravity();
		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			atom.setX(atom.getX()-centre[0]);
			atom.setY(atom.getY()-centre[1]);
			atom.setZ(atom.getZ()-centre[2]);
		}
		this.calcSphereCoord();		
	return centre;
	}

	/*--------------------------------------------------------------------------*/

	// move whole structure to point center
	public void move2centre(double[] centre) {
		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			atom.setX(atom.getX()-centre[0]);
			atom.setY(atom.getY()-centre[1]);
			atom.setZ(atom.getZ()-centre[2]);
		}
		this.calcSphereCoord();

	}
	
	/*--------------------------------------------------------------------------*/
	// move whole structure to a certain distance
	public void move(double[] move) {
		for(Iterator i=this.iterator(); i.hasNext();)
			((StructureAtom)i.next()).move(move);
	}
	
	/*--------------------------------------------------------------------------*/

//	public void setEccpAtomType(String type){
//		for(Iterator i=this.iterator(); i.hasNext();)
//			((StructureAtom)i.next()).setEccpAtomType(type);
//	}
	/*--------------------------------------------------------------------------*/
	public void setMMFFatomType(String type){
		for(Iterator i=this.iterator(); i.hasNext();)
			((StructureAtom)i.next()).setMMFFatomType(type);
	}
	/*--------------------------------------------------------------------------*/
	
	public void calcSphereCoord() {
		
		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			double[] sphCoords = Mathematics.xyz2rtp(atom.getX(), atom.getY(), atom.getZ());
			atom.setR(sphCoords[0]);
			atom.setTheta(sphCoords[1]);
			atom.setPhi(sphCoords[2]);	

/*			double[] xyz = stat.rtp2xyz(atom.getR(), atom.getTheta(), atom.getPhi());
			if(Math.abs(atom.getX() - xyz[0]) > 0.001) {
				System.err.print("x wrong: "+atom.getX()+"\t"+xyz[0]+"\n");
			}
			if(Math.abs(atom.getY() - xyz[1]) > 0.001) {
				System.err.print("y wrong: "+atom.getY()+"\t"+xyz[1]+"\n");
			}
			if(Math.abs(atom.getZ() - xyz[2]) > 0.001) {
				System.err.print("z wrong: "+atom.getZ()+"\t"+xyz[2]+"\n");
			}
*/		}		
	}

	/*--------------------------------------------------------------------------*/
	/*
	 * calculates the minimum distance between all atoms in this StructureCoords
	 * object and all atoms in coords2
	 */
	public double distance(StructureCoords coords2, boolean includeAtomRadius) {
		double minDist = 99999.9;
		for(Iterator<StructureAtom> i=this.iterator(); i.hasNext();) {
			StructureAtom atom1 = i.next();
			for(Iterator<StructureAtom> j=coords2.iterator(); j.hasNext();) {
				StructureAtom atom2 = j.next();
				double dist = 0; 
				if(includeAtomRadius){
					dist = atom1.distance(atom2);
				}
				else{
					dist = Mathematics.distance(atom1, atom2);
				}
				
				if(minDist > dist){
					minDist = dist;
				}
			}
		}
		return minDist;
	}
	/*--------------------------------------------------------------------------*/

	public void recalcXYZ() {
		for(Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();

			double[] xyz = Mathematics.rtp2xyz(atom.getR(), atom.getTheta(), atom.getPhi());
			
			atom.setX(xyz[0]);
			atom.setY(xyz[1]);
			atom.setZ(xyz[2]);	
		}		
	}

	/*--------------------------------------------------------------------------*/
    public String toString (boolean printRadiusPotential) {

		StringBuffer output = new StringBuffer();
		for (Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();
			if(printRadiusPotential) output.append(atom.toString(atom.getRadius(), atom.getPotential())); 
			else output.append(atom.toString(atom.getOccupancy(), atom.getTemp())); 
		}

	return output.toString();
    }
	/*--------------------------------------------------------------------------*/
    public String toString () {

		StringBuffer output = new StringBuffer();
		for (Iterator i=this.iterator(); i.hasNext();) {
			StructureAtom atom = (StructureAtom)i.next();
			output.append(atom.toString(atom.getOccupancy(), atom.getTemp())); 
		}
	return output.toString();
    }

	/*--------------------------------------------------------------------------*/
	/**
	 * sort the spherical coordinates as first phi than theta, with stable sort
	 * algorithm merge sort
	 * @param fileName String
	 */

	public void sortRPT() {

		Collections.sort(this,new Comparator() {
			public int compare(Object a, Object b) {			
				StructureAtom atom1 = (StructureAtom) a;
				StructureAtom atom2 = (StructureAtom) b;
			
				if (atom1.getTheta() < atom2.getTheta()) {
					return -1;
				}				
				else if (atom1.getTheta() > atom2.getTheta()) {
					return 1;
				}
				else {
					return 0;
				}
			}
		});

		Collections.sort(this,new Comparator() {
			public int compare(Object a, Object b) {			
				StructureAtom atom1 = (StructureAtom) a;
				StructureAtom atom2 = (StructureAtom) b;
			
				if (atom1.getPhi() < atom2.getPhi()) {
					return -1;
				}				
				else if (atom1.getPhi() > atom2.getPhi()) {
					return 1;
				}
				else {
					return 0;
				}
			}
		});				
	}

	/*--------------------------------------------------------------------------*/
	/**
	 * @param fileName String
	 */

	public StructureAtom getAtom(int index) {
			return (StructureAtom)this.get(index);

	}
	/*--------------------------------------------------------------------------*/

	/**
	 * get spherical coordinates form all cartesian coordinates as a String object
	 * @param none
	 * @return - 
	 */
	public String outputSphericalCoords () {
		Locale.setDefault(Locale.US);
		NumberFormat decFormat = new DecimalFormat("0.000");

		StringBuffer output = new StringBuffer();
		
		for (int i = 0; i < this.size(); i++) {

			output.append("SPHCOORDS                     ");
			for (int j=decFormat.format(this.getAtom(i).getR()).length();j<8; j++) {
				output.append(" ");
			}
			output.append(decFormat.format(this.getAtom(i).getR()));
			
			for (int j=decFormat.format(this.getAtom(i).getTheta()).length();j<8; j++) {
				output.append(" ");
			}
			output.append(decFormat.format(this.getAtom(i).getTheta())+"\n");
			
			for (int j=decFormat.format(this.getAtom(i).getPhi()).length();j<8; j++) {
				output.append(" ");
			}
			output.append(decFormat.format(this.getAtom(i).getPhi()));
		}
		return output.toString();
	}

	/*--------------------------------------------------------------------------*/	
	public boolean containsAtom(StructureAtom atom) {

		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom2 = (StructureAtom)i.next();
			if(atom2.equalsAtom(atom)){
				return true;
			}							
		}
		return false;
	}
	
	/*--------------------------------------------------------------------------*/	
	public String getChainIds() {
	
		StringBuffer chainIds = new StringBuffer();
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			
			if (chainIds.indexOf(Character.toString(atom.getChainId())) == -1){
				chainIds.append(atom.getChainId());
			}
		}
		
	return chainIds.toString();
	}

	/*--------------------------------------------------------------------------*/	
	public double[] getCenterOfGeometry(){
		return this.centerOfGeometry();
	}

	/*--------------------------------------------------------------------------*/	
	public double[] getCenterOfGravity(){
		return this.centerOfGravity();
	}

	/*--------------------------------------------------------------------------*/	
	public void setRadiusFromTempCol() {
	
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setRadius(atom.getTemp());
		}
	}
	
	/*--------------------------------------------------------------------------*/	
	public void setRadiusFromOccupancyCol() {
	
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setRadius(atom.getOccupancy());
		}
	}
		
	//----------------------------------------------------------------------------//
	public Hashtable<StructureAtom, StructureCoords> vibrate(boolean doVibrate){		
		// vibrating ligand atoms to calculate avg potential
		Hashtable vibratedCoords = new Hashtable();
		for(Iterator e=this.iterator(); e.hasNext();){
			StructureAtom atom = (StructureAtom)e.next();
						
			if(doVibrate){
				// infer scale of vibration from temperature column value B = 8*PI*(mean-square displacement)^2 = 79*(mean-square displacement)^2
				double vibrateScale = Math.sqrt(atom.getTemp()/79);

				// if vibration scale is set to 0.0 due to temperature factor is 0.0, than vibrate 0.25Å. Temperature factor can be 0.0 in cases
				// in which e.g. hydrogen atoms are added artifically
				if(vibrateScale <= 0.000009)
					vibrateScale = 0.25;

				StructureCoords vibrate = new StructureCoords();
				for(double x=-vibrateScale; x<=vibrateScale; x+=vibrateScale){
					for(double y=-vibrateScale; y<=vibrateScale; y+=vibrateScale){
						for(double z=-vibrateScale; z<=vibrateScale; z+=vibrateScale){
							StructureAtom atomCopy = atom.copy();
							double[] xyz = {x,y,z}; 
							atomCopy.move(xyz);
							vibrate.add(atomCopy);
						}
					}
				}
				vibratedCoords.put(atom, vibrate);
			}
			else{
				StructureCoords coords = new StructureCoords();
				coords.add(atom);
				vibratedCoords.put(atom, coords);				
			}
		}
		return vibratedCoords;
	}

	//--------------------------------------------------------------------------------//
	public void avgPotentialAfterVibration(Hashtable ensembleVibratedCoordinates){

		Hashtable pots = new Hashtable();
		Hashtable<Integer, Integer> count = new Hashtable();
		for(Enumeration j=ensembleVibratedCoordinates.keys(); j.hasMoreElements();){			
			StructureAtom coordsAtom = (StructureAtom)j.nextElement();
			StructureCoords vibratedCoords = (StructureCoords)ensembleVibratedCoordinates.get(coordsAtom);
			for(Iterator i=vibratedCoords.iterator(); i.hasNext();){
				StructureAtom coordsPoint = (StructureAtom)i.next();
					
				int serial = coordsPoint.getSerial();
				if(pots.get(serial)==null){
					pots.put(serial, coordsPoint.getPotential());
					count.put(serial, 1);
				}
				else{
					pots.put(serial, ((Double)pots.get(serial))+coordsPoint.getPotential());
					count.put(serial, count.get(serial)+1);
				}
			}
			for(Enumeration e=pots.keys(); e.hasMoreElements();){
				int serial = (Integer)e.nextElement();
				double pot = (Double)pots.get(serial)/count.get(serial);
				this.getAtomWithSerial(serial).setPotential(pot);
			}

		}
	}
	//--------------------------------------------------------------------------------//
	public void avgPotentialAfterVibration(StructureCoords ensembleVibratedCoordinates){
		// no hashtable as parameter required as avg value is sampled 
		
		Hashtable pots = new Hashtable();
		Hashtable<Integer, Integer> count = new Hashtable();
		for(Iterator i=ensembleVibratedCoordinates.iterator(); i.hasNext();){
			StructureAtom surfacePoint = (StructureAtom)i.next();
			
			int serial = surfacePoint.getSerial();
			if(pots.get(serial)==null){
				pots.put(serial, surfacePoint.getPotential());
				count.put(serial, 1);
			}
			else{
				pots.put(serial, ((Double)pots.get(serial))+surfacePoint.getPotential());
				count.put(serial, count.get(serial)+1);
			}
		}
		for(Enumeration e=pots.keys(); e.hasMoreElements();){
			int serial = (Integer)e.nextElement();
			double pot = (Double)pots.get(serial)/count.get(serial);
			
			this.getAtomWithSerial(serial).setPotential(pot);
		}
	}

	/*--------------------------------------------------------------------------*/
	public StructureAtom getAtomWithSerial(int serial) {
		for(Iterator<StructureAtom> i=this.iterator(); i.hasNext();){
			StructureAtom atom = i.next();
			if(atom.getSerial()==serial) return atom;
		}
		return null;
	}
	//--------------------------------------------------------------------------------------------------------------//

	public void setPot2Temp() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setPotential(atom.getTemp());
		}
	}
	//--------------------------------------------------------------------------------------------------------------//

	public void setPot2partialCharge() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setPotential(atom.getPartialCharge());
		}
	}
	//--------------------------------------------------------------------------------------------------------------//
	public void setPot2XlogP() {	
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setPotential(atom.getXlogP());
		}
	}
	//--------------------------------------------------------------------------------------------------------------//
	public void setTemp2ColourCode() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setTemp(atom.getColour());
		}
	}
	/*--------------------------------------------------------------------------*/	
	public void setOccupancyFromRadius() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setOccupancy(atom.getRadius());
		}
	}

	/*--------------------------------------------------------------------------*/	
	public void setSurfnetRadius() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setSurfnetRadius();
		}
	}
	/*--------------------------------------------------------------------------*/	
	public void setClusterPos(int pos) {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setClusterPos(pos);
		}
	}
	/*--------------------------------------------------------------------------*/	
	public int getClusterPos() {
		return this.getAtom(0).getClusterPos();
	}
	/*--------------------------------------------------------------------------*/	
	public void setVdwPlusHydrogenRadius() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setVdwPlusHydrogenRadius();
		}
	}
	/*--------------------------------------------------------------------------*/	
	public void setName(String name) {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setName(name);
		}
	}
	/*--------------------------------------------------------------------------*/	
	public void setFlag(String flag) {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setFlag(flag);
		}
	}
	/*--------------------------------------------------------------------------*/	
	public void setHybridizationState(int hybridizationState) {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setHybridizationState(hybridizationState);
		}
	}
	
	/*--------------------------------------------------------------------------*/	
	public void setRadius(double radius) {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setRadius(radius);
		}
	}
	/*--------------------------------------------------------------------------*/	
	public void increaseRadiusBy(double increaseValue) {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.increaseRadiusBy(increaseValue);
		}
	}
	
	/*--------------------------------------------------------------------------*/	
	public void setRadius() {
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			atom.setRadius();
		}
	}
	/*--------------------------------------------------------------------------*/	
	public StructureCoords getHydrogens() {
	
		ArrayList hydrogens = new ArrayList();
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			
			if(atom.getName().trim().startsWith("H")){
				hydrogens.add(atom);
			}
		}
		StructureCoords hCoords = new StructureCoords();
		hCoords.addAll(hydrogens);
		return hCoords;
	}
	/*--------------------------------------------------------------------------*/	
	public StructureCoords getHeavyAtoms() {
	
		StructureCoords heavyCoords = new StructureCoords();
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			
			if(!atom.getName().trim().startsWith("H")){
				heavyCoords.add(atom);
			}
		}
		return heavyCoords;
	}
	/*--------------------------------------------------------------------------*/	
	public StructureResidue getResidue(int pos){

		StructureAtom preAtom = this.getAtom(0);
		StructureCoords coords = new StructureCoords();
		
		int p = 0;
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			if(!atom.isEqualResidue(preAtom)){
				preAtom = atom;
				p++;
			}
			if(p==pos){
				coords.add(atom);
			}
		}
		if(coords.size()==0){
			System.err.print(this+"Could not find for protein residue in position "+pos+"\n");
			return null;
		}

		return new StructureResidue(coords);
	}
	/*--------------------------------------------------------------------------*/	
	public StructureResidue[] getResidues(){

		StructureAtom preAtom = this.getAtom(0);
		StructureResidue coords = new StructureResidue();
		Vector resiTmp = new Vector();
		
		for(Iterator i=this.iterator(); i.hasNext();){
			StructureAtom atom = (StructureAtom)i.next();
			if(!atom.isEqualResidue(preAtom)){
				resiTmp.add(coords);
				preAtom = atom;
				coords = new StructureResidue();
				coords.add(atom);				
			}
			else{
				coords.add(atom);
			}
		}
		resiTmp.add(coords);
		StructureResidue[] resi = new StructureResidue[resiTmp.size()];
		for(int i=0; i<resi.length; i++){
			resi[i] = (StructureResidue)resiTmp.get(i);
		}
	return resi;
	}

	/*--------------------------------------------------------------------------*/
	public void add(StructureAtom atom){
		super.add(atom);
	}
	
	//----------------------------------------------------------------------------//

	public void setOriginalCenter1(double[] originalCentre){
		this.originalCentre1 = originalCentre;
	}
	//----------------------------------------------------------------------------//

	public void setOriginalCenter2(double[] originalCentre){
		this.originalCentre2 = originalCentre;
	}
	/*--------------------------------------------------------------------------*/
	public double[] getOriginalCenter1(){
		return this.originalCentre1;
	}

	/*--------------------------------------------------------------------------*/
	public double[] getOriginalCenter2(){
		return this.originalCentre2;
	}

	//----------------------------------------------------------------------------//

	public String getXYZforMultivalue(){
		StringBuffer output = new StringBuffer(); 
		for(Iterator j=this.iterator(); j.hasNext();){			
			StructureAtom atom = (StructureAtom)j.next();
			output.append(atom.getX()+","+atom.getY()+","+atom.getZ()+"\n");
		}
	return output.toString();
	}
	
	//----------------------------------------------------------------------------//

	public void resetPotentialValues(){
		for(Iterator j=this.iterator(); j.hasNext();) ((StructureAtom)j.next()).setPotential(0.00);
	}
	//-----------------------------------------------------------------------//

	public void setCovalentBoundAtomsOfHydrogens(){
		for(Iterator j=this.iterator(); j.hasNext();){
			StructureAtom atom = (StructureAtom)j.next();
			if(atom.getName().trim().startsWith("H")){
				StructureCoords heavyAtoms = atom.getCovalentBondAtoms(this, true);
				if(heavyAtoms.size()==1){
					atom.setCovalentBondAtom(heavyAtoms.getAtom(0));
				}
				else if(heavyAtoms.size()!=1){
					double minDist = 999999;
					StructureAtom minAtom = new StructureAtom();
					for(Iterator i=heavyAtoms.iterator(); i.hasNext();){
						StructureAtom heavyAtom = (StructureAtom)i.next();
						double dist = Mathematics.distance(atom, heavyAtom);
						if(dist<minDist){
							minDist = dist;
							minAtom = heavyAtom;
						}
					}					
					atom.setCovalentBondAtom(minAtom);
				}
			}
		}
	}

	//----------------------------------------------------------------------------//

	public static void main(String args[]) {
		
	} // End of main method
} // end of class StructureCoords
/*

ATOM

Overview

The ATOM records present the atomic coordinates for standard residues. They
also present the occupancy and temperature factor for each atom. Heterogen
coordinates use the HETATM record type. The element symbol is always present
on each ATOM record; segment identifier and charge are optional.

Record Format

COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM"

 7 - 11        Integer         serial        Atom serial number.

13 - 16        Atom            name          Atom name.

17             Character       proteinAltLoc        Alternate location indicator.

18 - 20        Residue name    resName       Residue name.

22             Character       chainID       Chain identifier.

23 - 26        Integer         resSeq        Residue sequence number.

27             AChar           iCode         Code for insertion of residues.

31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
                                             Angstroms.

39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
                                             Angstroms.

47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
                                             Angstroms.

55 - 60        Real(6.2)       occupancy     Occupancy.

61 - 66        Real(6.2)       tempFactor    Temperature factor.

73 - 76        LString(4)      segID         Segment identifier, left-justified.

77 - 78        LString(2)      element       Element symbol, right-justified.

79 - 80        LString(2)      charge        Charge on the atom.


Example

         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C
ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C



ANISOU

Overview

The ANISOU records present the anisotropic temperature factors.

Record Format

COLUMNS        DATA TYPE       FIELD         DEFINITION
----------------------------------------------------------------------
 1 -  6        Record name     "ANISOU"

 7 - 11        Integer         serial        Atom serial number.

13 - 16        Atom            name          Atom name.

17             Character       proteinAltLoc        Alternate location
                                             indicator.

18 - 20        Residue name    resName       Residue name.

22             Character       chainID       Chain identifier.

23 - 26        Integer         resSeq        Residue sequence number.

27             AChar           iCode         Insertion code.

29 - 35        Integer         u[0][0]       U(1,1)

36 - 42        Integer         u[1][1]       U(2,2)

43 - 49        Integer         u[2][2]       U(3,3)

50 - 56        Integer         u[0][1]       U(1,2)

57 - 63        Integer         u[0][2]       U(1,3)

64 - 70        Integer         u[1][2]       U(2,3)

73 - 76        LString(4)      segID         Segment identifier, left-justified.

77 - 78        LString(2)      element       Element symbol, right-justified.

79 - 80        LString(2)      charge        Charge on the atom.


Example

         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    107  N   GLY    13      12.681  37.302 -25.211 1.000 15.56           N
ANISOU  107  N   GLY    13     2406   1892   1614    198    519   -328       N
ATOM    108  CA  GLY    13      11.982  37.996 -26.241 1.000 16.92           C
ANISOU  108  CA  GLY    13     2748   2004   1679    -21    155   -419       C
ATOM    109  C   GLY    13      11.678  39.447 -26.008 1.000 15.73           C
ANISOU  109  C   GLY    13     2555   1955   1468     87    357   -109       C
ATOM    110  O   GLY    13      11.444  40.201 -26.971 1.000 20.93           O
ANISOU  110  O   GLY    13     3837   2505   1611    164   -121    189       O
ATOM    111  N   ASN    14      11.608  39.863 -24.755 1.000 13.68           N
ANISOU  111  N   ASN    14     2059   1674   1462     27    244    -96       N

HETATM

Overview

The HETATM records present the atomic coordinate records for atoms within
"non-standard" groups. These records are used for water molecules and atoms
presented in HET groups.

Record Format

COLUMNS        DATA TYPE       FIELD          DEFINITION
--------------------------------------------------------------------------------
 1 -  6        Record name     "HETATM"

 7 - 11        Integer         serial         Atom serial number.

13 - 16        Atom            name           Atom name.

17             Character       proteinAltLoc         Alternate location indicator.

18 - 20        Residue name    resName        Residue name.

22             Character       chainID        Chain identifier.

23 - 26        Integer         resSeq         Residue sequence number.

27             AChar           iCode          Code for insertion of residues.

31 - 38        Real(8.3)       x              Orthogonal coordinates for X.

39 - 46        Real(8.3)       y              Orthogonal coordinates for Y.

47 - 54        Real(8.3)       z              Orthogonal coordinates for Z.

55 - 60        Real(6.2)       occupancy      Occupancy.

61 - 66        Real(6.2)       tempFactor     Temperature factor.

73 - 76        LString(4)      segID          Segment identifier;
                                              left-justified.

77 - 78        LString(2)      element        Element symbol; right-justified.

79 - 80        LString(2)      charge         Charge on the atom.


Amino Acids

RESIDUE                     ABBREVIATION                SYNONYM
-----------------------------------------------------------------------------
Alanine                     ALA                         A
Arginine                    ARG                         R
Asparagine                  ASN                         N
Aspartic acid               ASP                         D
ASP/ASN ambiguous           ASX                         B
Cysteine                    CYS                         C
Glutamine                   GLN                         Q
Glutamic acid               GLU                         E
GLU/GLN ambiguous           GLX                         Z
Glycine                     GLY                         G
Histidine                   HIS                         H
Isoleucine                  ILE                         I
Leucine                     LEU                         L
Lysine                      LYS                         K
Methionine                  MET                         M
Phenylalanine               PHE                         F
Proline                     PRO                         P
Serine                      SER                         S
Threonine                   THR                         T
Tryptophan                  TRP                         W
Tyrosine                    TYR                         Y
Unknown                     UNK
Valine                      VAL                         V

Nucleic Acids

RESIDUE                                  ABBREVIATION
-----------------------------------------------------------------------
Adenosine                                  A
Modified adenosine                        +A
Cytidine                                   C
Modified cytidine                         +C
Guanosine                                  G
Modified guanosine                        +G
Inosine                                    I
Modified inosine                          +I
Thymidine                                  T
Modified thymidine                        +T
Uridine                                    U
Modified uridine                          +U
Unknown                                  UNK


Amino Acids

NAME                    CODE           FORMULA                 MOL. WT.
-----------------------------------------------------------------------------
Alanine                 ALA            C3 H7 N1 O2             89.09
Arginine                ARG            C6 H14 N4 O2            174.20
Asparagine              ASN            C4 H8 N2 O3             132.12
Aspartic acid           ASP            C4 H7 N1 O4             133.10
ASP/ASN ambiguous       ASX            C4 H71/2 N11/2 O31/2    132.61
Cysteine                CYS            C3 H7 N1 O2 S1          121.15
Glutamine               GLN            C5 H10 N2 O3            146.15
Glutamic acid           GLU            C5 H9 N1 O4             147.13
GLU/GLN ambiguous       GLX            C5 H91/2 N11/2 O31/2    146.64
Glycine                 GLY            C2 H5 N1 O2             75.07
Histidine               HIS            C6 H9 N3 O2             155.16
Isoleucine              ILE            C6 H13 N1 O2            131.17
Leucine                 LEU            C6 H13 N1 O2            131.17
Lysine                  LYS            C6 H14 N2 O2            146.19
Methionine              MET            C5 H11 N1 O2 S1         149.21
Phenylalanine           PHE            C9 H11 N1 O2            165.19
Proline                 PRO            C5 H9 N1 O2             115.13
Serine                  SER            C3 H7 N1 O3             105.09
Threonine               THR            C4 H9 N1 O3             119.12
Tryptophan              TRP            C11 H12 N2 O2           204.23
Tyrosine                TYR            C9 H11 N1 O3            181.19
Valine                  VAL            C5 H11 N1 O2            117.15
Undetermined            UNK            C5 H6 N1 O3             128.16

Nucleotides

NAME                    CODE           FORMULA                 MOL. WT.
------------------------------------------------------------------------------
Adenosine               A              C10 H14 N5 O7 P1        347.22
Cytidine                C              C9 H14 N3 O8 P1         323.20
Guanosine               G              C10 H14 N5 O8 P1        363.22
Inosine                 I              C10 H13 N4 08 P1        348.21
Thymidine               T              C10 H15 N2 08 P1        322.21
Uridine                 U              C9 H13 N2 09 P1         324.18


SOURCE: ftp://ftp.rcsb.org/pub/pdb/doc/format_descriptions/Contents_Guide_21.txt



COLUMNS DATA TYPE FIELD DEFINITION
------------------------------------------------------------------------------
       1 - 6 Record name "CONECT"
       7 - 11 Integer serial Atom serial number
       12 - 16 Integer serial Serial number of bonded atom
       17 - 21 Integer serial Serial number of bonded atom
       22 - 26 Integer serial Serial number of bonded atom
       27 - 31 Integer serial Serial number of bonded atom
       32 - 36 Integer serial Serial number of hydrogen bonded atom
       37 - 41 Integer serial Serial number of hydrogen bonded atom
       42 - 46 Integer serial Serial number of salt bridged atom
       47 - 51 Integer serial Serial number of hydrogen bonded atom
       52 - 56 Integer serial Serial number of hydrogen bonded atom
       57 - 61 Integer serial Serial number of salt bridge atom 

*/












