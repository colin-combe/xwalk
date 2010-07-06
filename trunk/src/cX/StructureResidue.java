package cX;

import java.util.Hashtable;
import java.util.Iterator;

public class StructureResidue extends StructureCoords {
		
	//----------------------------------------------------------------------------//
	
	public StructureResidue(){
		super();

	}	
	//----------------------------------------------------------------------------//
	
	public StructureResidue(StructureCoords residue){
		
		try{
			StructureAtom preAtom = residue.getAtom(0);
			for(Iterator i=residue.iterator(); i.hasNext();){
				StructureAtom atom = (StructureAtom)i.next();
				if(!atom.getResName().equals(preAtom.getResName()) &&
						atom.getResNo() != preAtom.getResNo()){
					throw new Exception();
				}
			}
		}
		catch (Exception e){
			System.err.print("Error in reading in Residue information. Data contains more than one residue information\n");
		}

		for(int i=0; i<residue.size(); i++){
			this.add(residue.getAtom(i));
		}
	}

	//----------------------------------------------------------------------------//

	public void setId(String id){
		this.remark = id;
	}
	//----------------------------------------------------------------------------//

	public String getId(){
		return this.remark;
	}
	//----------------------------------------------------------------------------//

	public String getName(){
		return this.getAtom(0).getResName();
	}

	//----------------------------------------------------------------------------//
		
	static public char get1LetterCode(String threeLetterCode){
		Hashtable singleLetterCode = new Hashtable();
		singleLetterCode.put("ALA", 'A');
		singleLetterCode.put("ARG", 'R');
		singleLetterCode.put("ASN", 'N');
		singleLetterCode.put("ASP", 'D');
		singleLetterCode.put("ASX", 'B');
		singleLetterCode.put("CYS", 'C');
		singleLetterCode.put("GLN", 'Q');
		singleLetterCode.put("GLU", 'E');
		singleLetterCode.put("GLX", 'Z');
		singleLetterCode.put("GLY", 'G');
		singleLetterCode.put("HIS", 'H');
		singleLetterCode.put("ILE", 'I');
		singleLetterCode.put("LEU", 'L');
		singleLetterCode.put("LYS", 'K');
		singleLetterCode.put("MET", 'M');
		singleLetterCode.put("PHE", 'F');
		singleLetterCode.put("PRO", 'P');
		singleLetterCode.put("SER", 'S');
		singleLetterCode.put("THR", 'T');
		singleLetterCode.put("TRP", 'W');
		singleLetterCode.put("TYR", 'Y');
		singleLetterCode.put("VAL", 'V');
		singleLetterCode.put("UNK", '-');
		return (Character)singleLetterCode.get(threeLetterCode);
	}	
	//----------------------------------------------------------------------------//
	
	public char get1LetterCode(){
		return this.get1LetterCode(this.getAtom(0).getResName());
	}	
	//----------------------------------------------------------------------------//
	
	static public String get3LetterCode(String singleLetterCode){
		Hashtable threeLetterCode = new Hashtable();
		threeLetterCode.put("A", "ALA");
		threeLetterCode.put("R", "ARG");
		threeLetterCode.put("N", "ASN");
		threeLetterCode.put("D", "ASP");
		threeLetterCode.put("B", "ASX");
		threeLetterCode.put("C", "CYS");
		threeLetterCode.put("Q", "GLN");
		threeLetterCode.put("E", "GLU");
		threeLetterCode.put("Z", "GLX");
		threeLetterCode.put("G", "GLY");
		threeLetterCode.put("H", "HIS");
		threeLetterCode.put("I", "ILE");
		threeLetterCode.put("L", "LEU");
		threeLetterCode.put("K", "LYS");
		threeLetterCode.put("M", "MET");
		threeLetterCode.put("F", "PHE");
		threeLetterCode.put("P", "PRO");
		threeLetterCode.put("S", "SER");
		threeLetterCode.put("T", "THR");
		threeLetterCode.put("W", "TRP");
		threeLetterCode.put("Y", "TYR");
		threeLetterCode.put("V", "VAL");
		threeLetterCode.put("-", "UNK");
			
		return (String)threeLetterCode.get(threeLetterCode);
	}	
}
