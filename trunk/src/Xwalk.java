

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Locale;
import java.util.Vector;

import cX.Commandline;
import cX.Grid;
import cX.Mathematics;
import cX.ReadFile;
import cX.StructureAtom;
import cX.StructureCoords;
import cX.StructureResidue;
import cX.WriteFile;

public class Xwalk {

	Hashtable<StructureAtom,Hashtable<StructureAtom,Hashtable<String,Double>>> xlinksHash = new Hashtable();
	ArrayList<Vector> xlinksArray;
	
	static String nl = System.getProperty("line.separator");
	static String alphabet			= " ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";

	static public String sequenceDistLabel 	= "Sequence Distance";  
	static public String euclideanDistLabel = "Euclidean Distance";  
	static public String gridDistLabel 		= "Solvent Path Distance";  
	
	static public String fileSeperator;
	String infile;
	boolean verbose;
	boolean verboseGrid;
	Grid grid;
	StructureResidue[] residues;
	StructureResidue[] residuesCopy;
	double gridCellLength;
	double templateRadius;
	boolean solvDist;
	String atomType1 = "";
	String atomType2 = "";
	String residueType1 = "";
	String residueType2 = "";
	boolean onlyHomo;
	double maxDist;
	
	
	public Xwalk(String infile, boolean solvDist, double gridCellLength, double templateRadius, boolean verbose, boolean verboseGrid){
		xlinksArray = new ArrayList();
		this.verbose = verbose;
		this.verboseGrid = verboseGrid;
		this.infile = infile;
		this.gridCellLength = gridCellLength;
		this.templateRadius = templateRadius;
		this.solvDist = solvDist;
		
		fileSeperator = File.separator; if((File.separator).equals("\\")) fileSeperator = "\\\\";

		StructureCoords coords = new StructureCoords(infile);
		coords.extractAtoms("ATOM  ", alphabet, alphabet, false, false, false);
		StructureAtom.setFF("surfnet");
		coords.setRadius();

		residues = coords.getResidues();
		
		StructureCoords coordsCopy = new StructureCoords(infile);
		coordsCopy.extractAtoms("ATOM  ", alphabet, alphabet, false, false, false);
		coordsCopy.setRadius();

		residuesCopy = coordsCopy.getResidues();
		
		if(this.solvDist) 
			this.grid = new Grid(coords, gridCellLength, 0, templateRadius);
	}
	
	//--------------------------------------------------------------------------------------------------------------------------------------//
	public void resetXlinks(){
		xlinksArray = new ArrayList();
		xlinksHash = new Hashtable();
	}
	//--------------------------------------------------------------------------------------------------------------------------------------//
	private Vector<Vector<StructureResidue>> findAllRelevantResidues(String aa1resNo, String aa2resNo,
											 						 String chainIds1, String chainIds2,
											 						 String altLocs1, String altLocs2){

		
		Vector<StructureResidue> candidate1 = new Vector();
		Vector<StructureResidue> candidate2 = new Vector();

		
		for(int i=0; i<this.residues.length; i++){
			StructureResidue residue = residues[i];
			residue.setClusterPos(i+1);
			if(
				((residueType1.equals("") &&
				 aa1resNo.indexOf("#"+residue.getAtom(0).getResNo()+"#")!=-1)
				 
				 ||

			     (!residueType1.equals("") && !aa1resNo.equals("-999") &&
				(this.residueType1.indexOf("#"+residue.getAtom(0).getResName()+"#")!=-1 &&
				 aa1resNo.indexOf("#"+residue.getAtom(0).getResNo()+"#")!=-1)) 
				 
				 ||

				 (!residueType1.equals("") && aa1resNo.equals("-999") &&
				  this.residueType1.indexOf("#"+residue.getAtom(0).getResName()+"#")!=-1))
			     
			     &&
			    
				(chainIds1.indexOf(Character.toString(residue.getAtom(0).getChainId()))!=-1 &&
				 altLocs1.indexOf(Character.toString(residue.getAtom(0).getAltLoc()))!=-1)){
				if(!this.atomType1.equals("")){
					StructureResidue tmp = new StructureResidue();
					for(Iterator<StructureAtom> j=residue.iterator(); j.hasNext();){
						StructureAtom atom = j.next();
						if(this.atomType1.indexOf("#"+atom.getName().trim()+"#")!=-1){
							tmp.add(atom);
						}
					}
					if(tmp.size()!=0){
						candidate1.add(tmp);
						if(!aa1resNo.equals("-999")) candidate2.add(tmp);
					}
					else{
						if(this.verbose) 
							System.err.print("WARNING: "+this.infile.replaceAll(".*"+fileSeperator, "")+"\tAtom \"-aa1 "+this.atomType1.replaceAll("#", "")+"\" not found in residue "+residue.getAtom(0).getResName()+residue.getAtom(0).getResNo()+residue.getAtom(0).getChainId()+nl);
					}
				}
				else{
					candidate1.add(residue);
					if(!aa1resNo.equals("-999")) candidate2.add(residue);
				}
			}

			if(
				((residueType2.equals("") &&
				 aa2resNo.indexOf("#"+residue.getAtom(0).getResNo()+"#")!=-1) 
							 
				 ||

			     (!residueType2.equals("") && !aa2resNo.equals("-999") &&
				(this.residueType2.indexOf("#"+residue.getAtom(0).getResName()+"#")!=-1 &&
				 aa2resNo.indexOf("#"+residue.getAtom(0).getResNo()+"#")!=-1)) 
							 
				 ||

				 (!residueType2.equals("") && aa2resNo.equals("-999") &&
				  this.residueType2.indexOf("#"+residue.getAtom(0).getResName()+"#")!=-1))
						     
			     &&
						    
				(chainIds2.indexOf(Character.toString(residue.getAtom(0).getChainId()))!=-1 &&
				 altLocs2.indexOf(Character.toString(residue.getAtom(0).getAltLoc()))!=-1)){
				if(!this.atomType2.equals("")){
					StructureResidue tmp = new StructureResidue();
					for(Iterator<StructureAtom> j=residue.iterator(); j.hasNext();){
						StructureAtom atom = j.next();
						if(this.atomType2.indexOf("#"+atom.getName().trim()+"#")!=-1){
							tmp.add(atom);
						}
					}
					if(tmp.size()!=0){
						if(!aa2resNo.equals("-999")) candidate1.add(tmp);
						candidate2.add(tmp);
					}
					else{
						if(this.verbose) 
							System.err.print("WARNING: "+this.infile.replaceAll(".*"+fileSeperator, "")+"\tAtom \"-aa2 "+this.atomType2.replaceAll("#", "")+"\" not found in residue "+residue.getAtom(0).getResName()+residue.getAtom(0).getResNo()+residue.getAtom(0).getChainId()+nl);
					}
				}
				else{
					if(!aa2resNo.equals("-999")) candidate1.add(residue);
					candidate2.add(residue);
				}
			}
		}
		Vector<Vector<StructureResidue>> v = new Vector();
		v.add(candidate1);
		v.add(candidate2);
		return v;
	}
	
	//--------------------------------------------------------------------------------------------------------------------------------------//
	private Hashtable<StructureResidue, Vector<StructureResidue>> determinePairsToRelevantResidues(Vector<StructureResidue> candidate1, Vector<StructureResidue> candidate2, 
																					 boolean onlyIntra, boolean onlyInter, double minSasRatio){

		Hashtable<StructureResidue, Vector<StructureResidue>> hash = new Hashtable();
		// the 2nd hash serves as a means to determine whether a distance pair is redundant in a homomeric structure
		Hashtable<String,String> tmpHash = new Hashtable();
		for(Iterator<StructureResidue> i=candidate1.iterator(); i.hasNext();){
			StructureResidue residue1 = i.next();
			if(minSasRatio>0 && this.grid!=null){
				if(this.grid.getSolventAccessRatio(residue1, -1) < minSasRatio){
					continue;
				}
			}
			for(Iterator<StructureResidue> j=candidate2.iterator(); j.hasNext();){
				StructureResidue residue2 = j.next();
				if(onlyIntra && residue1.getAtom(0).getChainId()!=residue2.getAtom(0).getChainId()) continue;
				if(onlyInter && residue1.getAtom(0).getChainId()==residue2.getAtom(0).getChainId()) continue;
				if(minSasRatio>0 && this.grid!=null){
					if(this.grid.getSolventAccessRatio(residue2, -1) < minSasRatio){
						continue;
					}
				}
				
				if(residue1.distance(residue2, false)>this.maxDist) continue;
				String residueId1 = residue1.getAtom(0).getResName()+"-"+residue1.getAtom(0).getResNo();
				String residueId2 = residue2.getAtom(0).getResName()+"-"+residue2.getAtom(0).getResNo();
				if(!this.onlyHomo) residueId1 += "-"+residue1.getAtom(0).getChainId();
				if(!this.onlyHomo) residueId2 += "-"+residue2.getAtom(0).getChainId();
					
				// Xlinks to same residue makes obviously no sense
				if(residueId1.equals(residueId2)){ 
					continue;
				}
					
				String tmp1 =tmpHash.get(residue1);
				if(tmp1!=null)
					if(tmp1.indexOf("#"+residueId2+"#")==-1)
						tmpHash.put(residueId1, tmp1+"#"+residueId2+"#");							
					else
						continue;
				else
					tmpHash.put(residueId1, "#"+residueId2+"#");
				
				Vector<StructureResidue> v = hash.get(residue1);
				if(v==null){
					v = new Vector();
					v.add(residue2);
					hash.put(residue1, v);
				}
				else{
					v.add(residue2);
					hash.put(residue1, v);
				}
			}
		}
		
		return hash;
	}
	//--------------------------------------------------------------------------------------------------------------------------------------//
	private Hashtable<StructureResidue, Vector<StructureResidue>> makePairsUnique(Hashtable<StructureResidue, Vector<StructureResidue>> residuePairCandidates){

		Hashtable<String, Double> resHash = new Hashtable();
		Hashtable<StructureResidue, Vector<StructureResidue>> residuePairs = new Hashtable();
		
		for(Enumeration<StructureResidue> e = residuePairCandidates.keys(); e.hasMoreElements();){
			StructureResidue residue1 = e.nextElement();
			StructureAtom atom1 = residue1.getAtom(0);
			String residueId1 = atom1.getResName().trim()+"-"+atom1.getResNo();
			if(!this.onlyHomo) residueId1 += "-"+atom1.getChainId();
			residue1.setId(residueId1);
			
			Vector<StructureResidue> v1 = residuePairCandidates.get(residue1);
			for(Iterator<StructureResidue> i=v1.iterator(); i.hasNext();){
				StructureResidue residue2 = i.next();
				StructureAtom atom2 = residue2.getAtom(0);
				String residueId2 = atom2.getResName().trim()+"-"+atom2.getResNo();
				if(!this.onlyHomo) residueId2 += "-"+atom2.getChainId();
				residue2.setId(residueId2);
						
				double dist = Mathematics.distance(atom1, atom2);
				
				if(resHash.get("#"+residue1.getId()+"-"+residue2.getId()+"#")==null && 
				   resHash.get("#"+residue2.getId()+"-"+residue1.getId()+"#")==null){
										
					Vector<StructureResidue> v = new Vector();
					v.add(residue1);
					v.add(residue2);
					Collections.sort(v, new Comparator() {
					    public int compare(Object o1, Object o2) {
							StructureResidue residue1 = ((StructureResidue) o1);
							StructureAtom atom1 = residue1.getAtom(0);
							String residueId1 = atom1.getResName().trim()+"-"+atom1.getResNo()+"-"+atom1.getChainId();
							StructureResidue residue2 = ((StructureResidue) o2);
							StructureAtom atom2 = residue2.getAtom(0);
							String residueId2 = atom2.getResName().trim()+"-"+atom2.getResNo()+"-"+atom2.getChainId();
							return residueId1.compareTo(residueId2);
					    }
					});
					
					resHash.put("#"+v.get(0).getId()+"-"+v.get(1).getId()+"#", dist);

					Vector<StructureResidue> f = residuePairs.get(v.get(0));
					if(f==null){
						f = new Vector();
						f.add(v.get(1));
						residuePairs.put(v.get(0), f);

					}
					else{
						f.add(v.get(1));
						residuePairs.put(v.get(0), f);
					}
				}

				else{
					// replace object in residuePairs if current pair has smaller distance
					double oldDist = -1;
					if(resHash.get("#"+residueId1+"-"+residueId2+"#")!=null)
						oldDist = resHash.get("#"+residueId1+"-"+residueId2+"#");
					if(resHash.get("#"+residueId2+"-"+residueId1+"#")!=null)
						oldDist = resHash.get("#"+residueId2+"-"+residueId1+"#");
					
					if(dist < oldDist){
						// need to search for residue ID as residue objects are different for homo parameter.
						boolean found = false;
						for(Enumeration<StructureResidue> e2 = residuePairs.keys(); e2.hasMoreElements();){

							StructureResidue r3 = e2.nextElement();

							if(r3.getId().equals(residue1.getId()) || r3.getId().equals(residue2.getId())){
								Vector<StructureResidue> v = residuePairs.get(r3);
								for(int j=0; j<v.size(); j++){
									StructureResidue r4 = v.get(j);

									if(r3.getId().equals(residue1.getId()) && r4.getId().equals(residue2.getId())){
										// automatically deletes also residue from residueHash
										v.remove(r4);
										if(v.size()==0) residuePairs.remove(r3);
										resHash.put("#"+r3.getId()+"-"+r4.getId()+"#", dist);
										Vector<StructureResidue> newV = residuePairs.get(residue1);
										if(newV==null){
											newV = new Vector();
											newV.add(residue2);
											residuePairs.put(residue1, newV);
										}
										else{
											newV.add(residue2);											
											residuePairs.put(residue1,newV);
										}
										found = true;
										break;
									}
									else if(r3.getId().equals(residue2.getId()) && r4.getId().equals(residue1.getId())){
										v.remove(r4);
										if(v.size()==0)	residuePairs.remove(r3);
										resHash.put("#"+r3.getId()+"-"+r4.getId()+"#",dist);
										Vector<StructureResidue> newV = residuePairs.get(residue2);
										if(newV==null){
											newV = new Vector();
											newV.add(residue1);
											residuePairs.put(residue2, newV);
										}
										else{
											newV.add(residue1);
											residuePairs.put(residue2, newV);											
										}
										found = true;
										break;
									}
								}
							}
							if(found) {
								break;
							}
						}
						if(!found)
							System.err.print("ERROR:"+residueId1+"\t"+residueId2+"\n");
					}
				}
			}
		}
		return residuePairs;
	}

	//--------------------------------------------------------------------------------------------------------------------------------------//
	private Hashtable<StructureResidue, Vector<StructureResidue>> findRelevantResiduePairs(
										String aa1resNo, String aa2resNo,
										String chainIds1, String chainIds2,
										String altLocs1, String altLocs2,
										boolean onlyIntra, boolean onlyInter, double minSasRatio){

		Vector<Vector<StructureResidue>> v = this.findAllRelevantResidues(aa1resNo, aa2resNo, chainIds1, chainIds2, altLocs1, altLocs2);
		Hashtable<StructureResidue, Vector<StructureResidue>> hash = this.determinePairsToRelevantResidues(v.get(0), v.get(1), onlyIntra, onlyInter, minSasRatio);
//		hash = this.makePairsUnique(hash);  
		return hash;
	}
	
	
	//--------------------------------------------------------------------------------------------------------------------------------------//
	public ArrayList<Vector> findXlinks(String atomType1, String atomType2,
										String residueType1, String residueType2,
										String aa1resNo, String aa2resNo,
										String chainIds1, String chainIds2,
										String altLocs1, String altLocs2,
										boolean onlyIntra, boolean onlyInter, boolean onlyHomo,
										double maxDist, double minSasRatio){
		this.onlyHomo = onlyHomo;
		this.maxDist = maxDist;
		this.atomType1 = atomType1;
		this.atomType2 = atomType2;
		this.residueType1 = residueType1;
		this.residueType2 = residueType2;
		
		Hashtable<StructureResidue, Vector<StructureResidue>> candidates = this.findRelevantResiduePairs(aa1resNo, aa2resNo, 
																										 chainIds1, chainIds2, 
																										 altLocs1, altLocs2, 
																										 onlyIntra, onlyInter,
																										 minSasRatio);
		this.assignDistances(candidates, minSasRatio);
		this.sort();
		return this.xlinksArray;
	}
	
	//--------------------------------------------------------------------------------------------------------------------------------------//
	private void assignDistances(Hashtable<StructureResidue, Vector<StructureResidue>> candidates, double minSasRatio){
		Locale.setDefault(Locale.US);
		NumberFormat decFormat = new DecimalFormat("0.0");
		
		// reset all grid cells of residue pairs in order to allow distance calculation.
		for(Enumeration<StructureResidue> e=candidates.keys(); e.hasMoreElements();){
			StructureResidue residue1 = e.nextElement();
			Vector<StructureResidue> residues2 = candidates.get(residue1);
			if(this.solvDist){					

				this.grid.reset4DistCalc();
				this.grid.treatCoordAsSolvent(residue1);
				this.grid.treatCoordAsSolvent(residues2);
				this.grid.setSurfaceDistances(residue1, this.maxDist, this.templateRadius, minSasRatio, this.verboseGrid);
			}
			StructureAtom atom1 = residue1.getAtom(0);
				
			Hashtable<StructureAtom, Hashtable<String,Double>> residues2hash = new Hashtable();
			for(Iterator<StructureResidue> j=residues2.iterator(); j.hasNext();){
				StructureResidue residue2 = j.next();
				StructureAtom atom2 = residue2.getAtom(0);
										
				double seqDist = Math.abs(residue2.getClusterPos()-residue1.getClusterPos());
				double euclDist = residue1.distance(residue2, false);

				Hashtable<String, Double> dists = new Hashtable();
				dists.put(this.sequenceDistLabel, seqDist);					
				dists.put(this.euclideanDistLabel, Double.valueOf(decFormat.format(euclDist)));
				dists.put(this.gridDistLabel,-1.0);
				if(this.solvDist){
					double gridDist = grid.getSurfaceDistances(residue2);
					if(euclDist > gridDist){
						gridDist = euclDist;
					}
					if(gridDist>this.maxDist) {
						continue;
					}
					else
						dists.put(this.gridDistLabel, Double.valueOf(decFormat.format(gridDist)));
				}
				residues2hash.put(atom2, dists);
			}
			if(residues2hash.size()>0)
				this.xlinksHash.put(atom1, residues2hash);
		}
	}

	//-------------------------------------------------------------------------------------------------------------------------------------//
	private Hashtable<String, Double> getDistHashOfShortestDist(StructureAtom atom1, StructureAtom atom2){
		Hashtable<StructureAtom, Hashtable<String,Double>> pair1 = xlinksHash.get(atom1);
		Hashtable<StructureAtom, Hashtable<String,Double>> pair2 = xlinksHash.get(atom2);
		Hashtable<String, Double> distHash = null;
		Hashtable<String, Double> distHash2 = null;		
		if(pair1!=null && pair2!=null){
			distHash = pair1.get(atom2);
			distHash2 = pair2.get(atom1);

			if(distHash!=null && distHash2!=null){
				if(this.solvDist){
					if(distHash.get(Xwalk.gridDistLabel)>distHash2.get(Xwalk.gridDistLabel)){
						distHash = distHash2;
					}							
				}
				else{
					if(distHash.get(Xwalk.euclideanDistLabel)>distHash2.get(Xwalk.euclideanDistLabel)){
						distHash = distHash2;
					}
				}
			}
		}
		else if(pair1!=null){
			distHash = pair1.get(atom2);
		}
		else if(pair2!=null){
			distHash = pair2.get(atom1);
		}
		return distHash;
	}
	
	//-----------------------------------------------------------------------------------------------------------------//
	
	private String[] getSortedResidueIds(StructureAtom atom1, StructureAtom atom2){
		Vector<StructureAtom> v = new Vector();
		v.add(atom1);
		v.add(atom2);
		Collections.sort(v, new Comparator() {
		    public int compare(Object o1, Object o2) {
				StructureAtom atom1 = ((StructureAtom) o1);
				char c1 = atom1.getChainId();if(c1==' ') c1='_';
				String residueId1 = c1+"-"+atom1.getResNo()+"-"+atom1.getResName().trim();
				
				StructureAtom atom2 = ((StructureAtom) o2);
				char c2 = atom2.getChainId();if(c2==' ') c2='_';
				String residueId2 = c2+"-"+atom2.getResNo()+"-"+atom2.getResName().trim();
				
				return residueId1.compareTo(residueId2);
		    }
		});
		atom1 = v.get(0);
		atom2 = v.get(1);

		String[] residueIds = new String[2];
		residueIds[0] = atom1.getResName().trim()+"-"+atom1.getResNo();
		if(atom1.getChainId()==' ') atom1.setChainId('_');
		residueIds[0] += "-"+atom1.getChainId()+"-"+atom1.getName().trim();

		residueIds[1] = atom2.getResName()+"-"+atom2.getResNo();
		if(atom2.getChainId()==' ') atom2.setChainId('_');
		residueIds[1] += "-"+atom2.getChainId()+"-"+atom2.getName().trim();
	return residueIds;
	}

	//--------------------------------------------------------------------------------------------------------------------------------------//
	private ArrayList<Vector> sortElemHash(Hashtable<String,Vector> elemHash){
		ArrayList<Vector> dists = new ArrayList();
		for(Enumeration<Vector> e = elemHash.elements(); e.hasMoreElements();){
			dists.add(e.nextElement());
		}
		Collections.sort(dists,new Comparator(){
			public int compare(Object o1, Object o2){
				Vector v1 = ((Vector) o1);
				Vector v2 = ((Vector) o2);
				Double dist1 = (Double)v1.get(v1.size()-1);
				Double dist2 = (Double)v2.get(v2.size()-1);

				// sort such that smallest distance is at the beginning
				if (dist1 > dist2) {
					return 1;
				}				
				else if (dist1 < dist2) {
					return -1;
				}
				else {
					return 0;
				}
			}
		});
		return dists;
	}
	
	//--------------------------------------------------------------------------------------------------------------------------------------//
	public void sort(){
		Hashtable<String,Double> nonRedundant = new Hashtable();
		Hashtable<String,Vector> elemHash = new Hashtable();
		Hashtable<String,String> rel = new Hashtable();
		
		for(Enumeration<StructureAtom> e = this.xlinksHash.keys(); e.hasMoreElements();){
			StructureAtom atom1 = e.nextElement();
			for(Enumeration<StructureAtom> f=xlinksHash.get(atom1).keys(); f.hasMoreElements();){
				StructureAtom atom2 = f.nextElement(); 

				// before we go on, we first need to check whether the reverse distance is shorter, in which case we continue with the reverse distance
				Hashtable<String, Double> distHash = this.getDistHashOfShortestDist(atom1, atom2);
				
				String[] residueIds = this.getSortedResidueIds(atom1, atom2);
				Vector elem = new Vector();
				elem.add(residueIds[0]);
				elem.add(residueIds[1]);
				elem.add(distHash.get(Xwalk.sequenceDistLabel).intValue());
				elem.add(distHash.get(Xwalk.euclideanDistLabel));
				double dist = distHash.get(Xwalk.euclideanDistLabel);
				if(this.solvDist){
					elem.add(distHash.get(Xwalk.gridDistLabel));
					dist = distHash.get(Xwalk.gridDistLabel);
				}
				
				String[] residueIdArray1 = residueIds[0].split("-");
				String resName1 = residueIdArray1[0];
				String resNo1 = residueIdArray1[1];
				String chainId1 = residueIdArray1[2];
				String atomName1 = residueIdArray1[3];
				String[] residueIdArray2 = residueIds[1].split("-");
				String resName2 = residueIdArray2[0];
				String resNo2 = residueIdArray2[1];
				String chainId2 = residueIdArray2[2];
				String atomName2 = residueIdArray2[3];

				String residueId1 = resName1+"-"+resNo1;
				if(!this.atomType1.equals("")) residueId1 += "-"+atomName1;
				if(!this.onlyHomo) residueId1 += "-"+chainId1;
				String residueId2 = resName2+"-"+resNo2;
				if(!this.atomType2.equals("")) residueId2 += "-"+atomName2;
				if(!this.onlyHomo) residueId2 += "-"+chainId2;
				
				if(nonRedundant.get(residueId1+"-"+residueId2)==null){
					nonRedundant.put(residueId1+"-"+residueId2, dist);
					elemHash.put(residueId1+"-"+residueId2, elem);
				}
				else{
					if(dist < nonRedundant.get(residueId1+"-"+residueId2)){
						nonRedundant.put(residueId1+"-"+residueId2, dist);							
						elemHash.put(residueId1+"-"+residueId2, elem);
					}
				}
			}
		}
		
		xlinksArray = this.sortElemHash(elemHash);
	}
	
	//--------------------------------------------------------------------------------------------------------------------------------------//	
	public String toString(boolean outputWarning){
		String filename = infile.replaceAll(".*"+Xwalk.fileSeperator, "");

		StringBuffer output = new StringBuffer();
		                 
		if(this.xlinksArray.size()!=0){
			for(int i=0; i<xlinksArray.size(); i++){
				output.append((i+1)+"\t"+filename+"\t");
				Vector v = xlinksArray.get(i);
				// second is residue ID
				String[] residueId1 = ((String)v.get(0)).split("-");
				output.append(residueId1[0]+"-"+residueId1[1]);
				if(!this.onlyHomo) output.append("-"+residueId1[2]);
				if(!this.atomType1.equals("")) output.append("-"+residueId1[3]);

				output.append("\t");

				String[] residueId2 = ((String)v.get(1)).split("-");
				output.append(residueId2[0]+"-"+residueId2[1]);
				if(!this.onlyHomo) output.append("-"+residueId2[2]);
				if(!this.atomType1.equals("")) output.append("-"+residueId2[3]);

				output.append("\t"+v.get(2));
				output.append("\t"+v.get(3));
				if(this.solvDist){
					output.append("\t"+v.get(4));
				}
				output.append(nl);
			}
			
		}
		else{
			if(outputWarning)
				
				System.err.print("WARNING: "+infile.replaceAll(".*"+fileSeperator, "")+"\tThere is no pair of "+this.residueType1.replaceAll("#", "")+"\t"+this.residueType2.replaceAll("#", "")+" residues that has a distance smaller than "+maxDist+"."+nl);
		}
		return output.toString();
	}
	
	//--------------------------------------------------------------------------------------------------------------------------------------//

	public String outputPymolScript(){
		Hashtable<Double, Vector<String>> hash = new Hashtable();

		StringBuffer output = new StringBuffer(); 
		output.append("load "+this.infile+nl);
		output.append("disable all"+nl);
		output.append("hide"+nl);
		output.append("set dash_radius, 1"+nl);
		output.append("bg_color white"+nl);
		output.append("util.cbc"+nl);
		output.append("create het, hetatm"+nl);
		output.append("show sticks, het"+nl);
		output.append("color grey, het"+nl); 
		output.append("disable het"+nl);
		
		boolean emptyChainId = false;
		
		if(this.xlinksArray.size()!=0){
			for(int i=0; i<xlinksArray.size(); i++){
				Vector v = xlinksArray.get(i);
				// second is residue ID
				String residueId1 = (String)v.get(0);
				String residueId2 = (String)v.get(1);

				String[] residueIdArray1 = residueId1.split("-");
				String resName1 = residueIdArray1[0];
				String resNo1 = residueIdArray1[1];
				String chainId1 = residueIdArray1[2];
				String atomName1 = residueIdArray1[3];

				String[] residueIdArray2 = residueId2.split("-");
				String resName2 = residueIdArray2[0];
				String resNo2 = residueIdArray2[1];
				String chainId2 = residueIdArray2[2];
				String atomName2 = residueIdArray2[3];
				
				if(chainId1.equals("_") || chainId2.equals("_")){
					emptyChainId = true;					
				}

				if(!chainId1.equals("_")) 
					output.append("create chain"+chainId1+", chain "+chainId1+nl);
				String selection1 = "\"resn "+resName1+" and resi "+resNo1+" and chain "+chainId1+" and name "+atomName1+"\"";
				
				if(!chainId2.equals("_")) 
					output.append("create chain"+chainId2+", chain "+chainId2+nl);
				String selection2 = "\"resn "+resName2+" and resi "+resNo2+" and chain "+chainId2+" and name "+atomName2+"\"";

				String distName = (i+1)+"_";
				if(this.solvDist) distName += v.get(v.size()-1)+"_";
				distName += residueId1+"_"+residueId2;

				output.append("cmd.select(\"pk1\","+selection1+")"+nl);
				output.append("cmd.select(\"pk2\","+selection2+")"+nl);
				output.append("cmd.show(\"spheres\",\"pk1\")"+nl);			
				output.append("cmd.show(\"spheres\",\"pk2\")"+nl);
				output.append("cmd.distance(\""+distName+"\", \"(pk1)\", \"(pk2)\")"+nl);
				output.append("cmd.color(\"red\", \""+distName+"\")"+nl); 

			}
			
			output.append("delete pk1"+nl); 
			output.append("delete pk2"+nl); 
			
			output.append("show ribbon"+nl);
			output.append("show surface"+nl);
			if(emptyChainId || onlyHomo) output.append("set transparency, 0.5"+nl);
			else output.append("set transparency, 0.5, chain*"+nl);
			output.append("center");
		}
		return output.toString();
	}
	

	//--------------------------------------------------------------------------------------------------------------------------------------//

	public static void main(String args[]) {
		Commandline commandline = new Commandline();
		String infile			= "";
		String outfile			= "";
		String distFile			= ""; 
		String residueType1  	= "";
		String residueType2 	= "";
		String aa1resNo  		= "-999";
		String aa2resNo  		= "-999";
		String chainIds1 		= Xwalk.alphabet;
		String chainIds2 		= Xwalk.alphabet;
		String altLocs1  		= Xwalk.alphabet;
		String altLocs2  		= Xwalk.alphabet;
		String atomType1 		= "";
		String atomType2 		= "";
		double maxDist   		= 30;
		boolean pymolOutput 	= false;
		boolean solvDist 		= true;
		double templateRadius 	= 3.0;
		double gridCellLength 	= 1.0;
		double minSasRatio		= -1.0;
		boolean homo 			= false;
		boolean showIntra 		= false;
		boolean showInter 		= false;
		boolean verbose 		= false;
		boolean verboseGrid 	= false;
		if (args.length == 0){
			System.out.print(nl+"Xwalk -in 1brs.pdb -aa1 LYS -aa2 lys"+nl+
							 "More informations with \"Xwalk -help\""+nl+nl);
			System.exit(0);
			}
		// if "-help" is written as parameter. a describtion about the parameter
		// will display

		
		if (commandline.get(args,"-help", false).equals("EXISTS") || commandline.get(args,"-h", false).equals("EXISTS")){
			System.out.print(nl+"EXAMPLARY command for program execution:" +nl+
							"java Xwalk -in 1brs.pdb -aa1 LYS#ARG -aa2 lys#arg -max 21"+nl+nl+
							"ABOUT"+nl+
							"Version 2.2"+nl+
							"Xwalk calculates and outputs distances in Angstroem for potential cross-links between -aa1 type amino acids and -aa2 type amino acids in the PDB file -in."+nl+nl+
							"IMPORTANT"+nl+
							"If large protein complexes are processed, the Java heap size might need to be increased from the default 64MB to 256MB, with the Java parameter -Xmx256m "+nl+nl+
							"OUTPUT FORMAT:"+nl +
							"IndexNo\tInfileName\tResidue1info\tResidue2info\tDistanceInPDBsequence\tEuclideanDistance\tSolventPathDistance"+nl +
							nl+
							"COMMANDLINE PARAMETER:"+nl+
							"INPUT/OUTPUT:"+nl+
							"\t-in\t<path>\tAny PDB file [required]."+nl +
							"\t-dist\t<path>\tAny Xwalk distance file, from which all residue information will be extracted [optional]."+nl +
					 		"\t-out\t<path>\tWrites output to this file, otherwise output is directed to the STDOUT channel. If -p is set than filename must have .pml filename ending [optional]."+nl +
							"\t-p\t[switch]\tOutputs a PyMOL (http://www.pymol.org/) script highlighting the calculated distances of the potential cross-links [optional]."+nl +
							"\t-v\t[switch]\tOutputs various information other than distances [optional]."+nl +
							"\t-vg\t[switch]\tOutputs on STDOUT channel the Solvent-Path-Distance grids in PDB format with distances in the B-factor column [optional]."+nl +
							nl+
							"RESIDUE/ATOM SELECTION:"+nl +
							"\t-aa1\t[String]\tThree letter code of 1st amino acid. To specify more than one amino acid use '#' as a delimeter [required, if -r1 or -dist is not set]."+nl +
							"\t-aa2\t[String]\tThree letter code of 2nd amino acid. To specify more than one amino acid use '#' as a delimeter [required, if -r2 or -dist is not set]."+nl +
							"\t-r1\t[String]\tAmino acid residue number. To specify more than one residue number use '#' as a delimeter. [required, if -aa1 or -dist is not set]."+nl +
							"\t-r2\t[String]\tAmino acid residue number. To specify more than one residue number use '#' as a delimeter. [required, if -aa2 or -dist is not set]."+nl +
							"\t-c1\t[String]\tChain ids for -aa1 or -r1. For blank chain Id use '_'. To specify more than one chain Id, append chain ids to a single string, e.g. ABC [optional](default: all chain Ids)."+nl +
							"\t-c2\t[String]\tChain ids for -aa2 or -r2. For blank chain Id use '_'. To specify more than one chain Id, append chain ids to a single string, e.g. ABC [optional](default: all chain Ids)."+nl +
							"\t-a1\t[String]\tAtom type for -aa1 or -r1. To specify more than one atom type use '#' as a delimeter. [optional]."+nl +
							"\t-a2\t[String]\tAtom type for -aa2 or -r2. To specify more than one atom type use '#' as a delimeter. [optional]."+nl +
							"\t-l1\t[String]\tAlternative location id for -aa1 or -r1. To specify more than one alternative location, append alternative location ids to a single string, e.g. AB [optional]."+nl +
							"\t-l2\t[String]\tAlternative location id for -aa2 or -r1. To specify more than one alternative location, append alternative location ids to a single string, e.g. AB [optional]."+nl +
							"\t-sas\t[double]\tMinimum solvent accessibility of both amino acids in % (recommended value 0.25 )."+nl +
							"\t-intra\t[switch]\tOutputs only \"intra-molecular\" distances [optional]."+nl +
							"\t-inter\t[switch]\tOutputs only \"inter-molecular\" distances [optional]."+nl +
							"\t-homo\t[double]\tOutputs only shortest distance of potential cross-links between equally numbered residues. Reduces redundancy if PDB file is a homomeric protein complex. [optional]."+nl +
							nl +
							"DISTANCE RELATED:"+nl +
							"\t-max\t[double]\tCalculates distances in Angstroem only up-to this value (default: 30.00)."+nl +
							"\t-euc\t[switch]\tSkips Solvent-Path-Distance calculation and outputs only Euclidean distances [optional]. "+nl +
							nl +
							"SOLVENT-PATH-DISTANCE GRID RELATED:"+nl +
							"\t-radius\t[double]\tCross-section radius in Angstroem of cross-linker [optional](default 3.0)."+nl +
							"\t-space\t[double]\tSpacing in Angstroem between grid cells. [optional](default 1.0)."+nl +
							nl);
			System.exit(0);
		}

		//----------------------------
		// getting parameter "-in"
		if (commandline.get(args,"-in", true).equals("ERROR")){
			System.err.print(nl+"ERROR: Could NOT find value for parameter \"-in1\"." +
			                 "For more information, please type \"Xwalk -help\""+nl+nl);
			System.exit(1);
		}
		else {		
			infile = commandline.get(args,"-in", true);
			try{
				new FileReader(infile);
			}
			catch(FileNotFoundException e){
				System.err.print(nl+"ERROR: File \""+infile+"\" NOT found!!!"+nl+nl);
				System.exit(2);
			}
		}
		//----------------------------
		// getting parameter "-p"
		if (commandline.get(args,"-p", false).equals("EXISTS")){
			pymolOutput = true;
			if(!commandline.get(args,"-out", true).equals("ERROR")){
				if(!commandline.get(args,"-out", true).endsWith(".pml")){
					System.err.print(nl+"ERROR: Please use the file ending \".pml\" for the output file \""+commandline.get(args,"-out", true)+"\""+nl+nl);
					System.exit(-1);
				}
			}
			else{
				System.err.print(nl+"WARNING: Please make sure that your filename ending is \".pml\", if STDOUT is redirected into a file. "+nl+nl);
			}
		}
		//----------------------------
		if (!commandline.get(args,"-out", true).equals("ERROR")){
			outfile = commandline.get(args,"-out", true);
			if(ReadFile.exists(outfile)){
				System.out.print(nl+"File \""+outfile+"\" already exists. Overwrite? [y/n]: ");
				String respond = Commandline.get();
				if(!respond.equalsIgnoreCase("y") && !respond.equalsIgnoreCase("yes")){
					System.exit(0);
				}
			}
		}
		else if(commandline.get(args,"-out", false).equals("EXISTS")){
			System.err.print(nl+"ERROR: Please specify output filename."+nl+nl);
			System.exit(1);
		}
		//----------------------------
		// getting parameter "-aa1"
		if (commandline.get(args,"-aa1", true).equals("ERROR")){
			if(commandline.get(args,"-r1", true).equals("ERROR")){
				if(commandline.get(args,"-dist", true).equals("ERROR")){
					System.err.print(nl+"ERROR: Could NOT find value neither for parameter \"-aa1\" nor for \"-r1\"." +
							"For more information, please type \"Xwalk -help\""+nl+nl);			
					System.exit(3);
				}
			}
		}
		else{
			residueType1 = "#"+commandline.get(args,"-aa1", true).toUpperCase()+"#";
		}
		//----------------------------
		// getting parameter "-aa2"
		if (commandline.get(args,"-aa2", true).equals("ERROR")){
			if(commandline.get(args,"-r2", true).equals("ERROR")){
				if(commandline.get(args,"-dist", true).equals("ERROR")){
					System.err.print(nl+"ERROR: Could NOT find value neither for parameter \"-aa2\" nor for \"-r2\"." +
							"For more information, please type \"Xwalk -help\""+nl+nl);			
					System.exit(3);
				}
			}
		}
		else{		
			residueType2 = "#"+commandline.get(args,"-aa2", true).toUpperCase()+"#";
		}
		//----------------------------
		// getting parameter "-c1/c2"
		if (!commandline.get(args,"-c1", true).equals("ERROR")){
			chainIds1 = commandline.get(args,"-c1", true).toUpperCase();
		}
		if (!commandline.get(args,"-c2", true).equals("ERROR")){
			chainIds2 = commandline.get(args,"-c2", true).toUpperCase();
		}
		//----------------------------
		// getting parameter "-l1/l2"
		if (!commandline.get(args,"-l1", true).equals("ERROR")){
			altLocs1 = commandline.get(args,"-l1", true).toUpperCase();
		}
		if (!commandline.get(args,"-l2", true).equals("ERROR")){
			altLocs2 = commandline.get(args,"-l2", true).toUpperCase();
		}

		//----------------------------

		if(commandline.get(args,"-r1", true).equals("ERROR")){
			if(residueType1.equals("")){
				if(commandline.get(args,"-dist", true).equals("ERROR")){
					System.err.print(nl+"ERROR: Could NOT find value neither for parameter \"-r1\" nor for \"-aa1\"." +
							"For more information, please type \"Xwalk -help\""+nl+nl);			
					System.exit(4);
				}
			}
		}
		else{
				aa1resNo = "#"+commandline.get(args,"-r1", true).toUpperCase()+"#";
				if(!residueType1.equals("")){
					residueType1 = "";
					System.err.print(nl+"WARNING: -aa1 and -r1 are both set. Disregarding -aa1 and considering only \"-r1 "+aa1resNo.replaceAll("#", "")+"\""+nl+nl);
				}
		}
		
		//----------------------------
		if(commandline.get(args,"-r2", true).equals("ERROR")){
			if(residueType2.equals("")){
				if(commandline.get(args,"-dist", true).equals("ERROR")){
					System.err.print(nl+"ERROR: Could NOT find value neither for parameter \"-r2\" nor for \"-aa2\"." +
							"For more information, please type \"Xwalk -help\""+nl+nl);			
					System.exit(4);
				}
			}
		}
		else{
			aa2resNo = "#"+commandline.get(args,"-r2", true).toUpperCase()+"#";
			if(!residueType2.equals("")){
				residueType2 = "";
				System.err.print(nl+"WARNING: -aa2 and -r2 are both set. Disregarding -aa1 and considering only \"-r2 "+aa2resNo.replaceAll("#", "")+"\""+nl+nl);
			}
		}
		
		//----------------------------
		// getting parameter "-a1"
		if (!commandline.get(args,"-a1", true).equals("ERROR")){
			atomType1 = "#"+commandline.get(args,"-a1", true).toUpperCase()+"#";
		}
		//----------------------------
		// getting parameter "-a2"
		if (!commandline.get(args,"-a2", true).equals("ERROR")){
			atomType2 = "#"+commandline.get(args,"-a2", true).toUpperCase()+"#";
		}
		//----------------------------
		// getting parameter "-max"
		if (!commandline.get(args,"-max", true).equals("ERROR")){
			maxDist = Double.parseDouble(commandline.get(args,"-max", true));
		}
		//----------------------------
		// getting parameter "-euc"
		if (commandline.get(args,"-euc", false).equals("EXISTS")){
			solvDist = false;
		}
		//----------------------------
		// getting parameter "-radius"
		if (!commandline.get(args,"-radius", true).equals("ERROR")){
			templateRadius = Double.parseDouble(commandline.get(args,"-radius", true));
		}
		//----------------------------
		// getting parameter "-space"
		if (!commandline.get(args,"-space", true).equals("ERROR")){
			gridCellLength = Double.parseDouble(commandline.get(args,"-space", true));
		}
		//----------------------------
		// getting parameter "-sas"
		if (!commandline.get(args,"-sas", true).equals("ERROR")){
			if(!solvDist){
				System.err.print(nl+"ERROR: Please set \"-dist\" in order to use \"-sas\""+nl+nl);
				System.exit(8);
			}
			minSasRatio = Double.parseDouble(commandline.get(args,"-sas", true));;
		}
		
		//----------------------------
		// getting parameter "-intra"
		if (commandline.get(args,"-intra", false).equals("EXISTS")){
			showIntra = true;
		}
		//----------------------------
		// getting parameter "-inter"
		if (commandline.get(args,"-inter", false).equals("EXISTS")){
			showInter = true;
		}
		//----------------------------
		// getting parameter "-homo"
		if (commandline.get(args,"-homo", false).equals("EXISTS")){
			homo = true;
		}
		//----------------------------
		// getting parameter "-v"
		if (commandline.get(args,"-v", false).equals("EXISTS")){
			verbose = true;
		}
		//----------------------------
		// getting parameter "-vg"
		if (commandline.get(args,"-vg", false).equals("EXISTS")){
			verboseGrid = true;
		}
		//----------------------------
		// getting parameter "-dist"
		if (!commandline.get(args,"-dist", true).equals("ERROR")){
			distFile = commandline.get(args,"-dist", true);
			try{
				new FileReader(distFile);
			}
			catch(FileNotFoundException e){
				System.err.print(nl+"ERROR: File \""+distFile+"\" NOT found!!!"+nl+nl);
				System.exit(2);
			}
		}

		// ------------------- END OF COMMANDLINE PARAMETER ANALYSIS --------------------------------------//
		
		Xwalk xwalk = new Xwalk(infile, solvDist, gridCellLength, templateRadius, verbose, verboseGrid);		

		StringBuffer output 		= new StringBuffer();
		StringBuffer outputPymol 	= new StringBuffer();
		
		output.append("#-----\t--------\t---------\t---------\t---\t---");
		if(solvDist) output.append("\t--------");
		output.append(nl);			
		output.append("#Index\tFileName\tResi1info\tResi2info\tSeq\tEuc");
		if(solvDist) output.append("\tSolvPath");
		output.append(nl);
		output.append("#-----\t--------\t---------\t---------\t---\t---");
		if(solvDist) output.append("\t--------");
		output.append(nl);

		if(distFile.equals("")){
			xwalk.findXlinks(atomType1, atomType2, 
							 residueType1, residueType2, 
							 aa1resNo, aa2resNo, 
							 chainIds1, chainIds2, 
							 altLocs1, altLocs2, 
							 showIntra, showInter, homo, 
							 maxDist, minSasRatio);
			output.append(xwalk.toString(true));
			outputPymol.append(xwalk.outputPymolScript());
		}
		else{
			ReadFile file = new ReadFile(distFile);

			for(Iterator<String> i=file.iterator(); i.hasNext();){
				String line = i.next();
				if(!line.startsWith("#")){
					String[] array = line.split("\t");
					int index = Integer.parseInt(array[0]);
					String fileName = array[1];
					String res1 = array[2];
					String res2 = array[3];
					double seqDist = Double.parseDouble(array[4]);
					double eucDist = Double.parseDouble(array[5]);
					double solDist = -1; 
					if(array.length>6){
						solDist = Double.parseDouble(array[6]);
					}
					array = res1.split("-");
					String resName1 = array[0];
					String resNo1 = array[1];
					char chainId1 = array[2].charAt(0);
					String atomName1 = array[3];

					array = res2.split("-");
					String resName2 = array[0];
					String resNo2 = array[1];
					char chainId2 = array[2].charAt(0);
					String atomName2 = array[3];
										
					residueType1 = "";
					residueType2 = "";
					aa1resNo = "";
					aa2resNo = "";
					chainIds1 = "";
					chainIds2 = "";
					atomType1 = "";
					atomType2 = "";

					residueType1 += "#"+resName1+"#";
					residueType2 += "#"+resName2+"#";
					aa1resNo += "#"+resNo1+"#";
					aa2resNo += "#"+resNo2+"#";
					chainIds1 += "#"+chainId1+"#";
					chainIds2 += "#"+chainId2+"#";
					atomType1 += "#"+atomName1+"#";
					atomType2 += "#"+atomName2+"#";

					xwalk.findXlinks(atomType1, atomType2, 
									 residueType1, residueType2, 
									 aa1resNo, aa2resNo, 
									 chainIds1, chainIds2, 
									 altLocs1, altLocs2, 
									 showIntra, showInter, homo, 
									 maxDist, minSasRatio);

					output.append(xwalk.toString(false).replaceAll("^\\d+", Integer.toString(index)));
					outputPymol.append(xwalk.outputPymolScript());
					
					xwalk.resetXlinks();
				}
			}
		}
		if(outfile.equals("")){
			if(pymolOutput)	System.out.print(outputPymol.toString());
			else System.out.print(output.toString());
		}
		else{
			WriteFile write = new WriteFile(outfile);
			if(pymolOutput)	write.write(outputPymol.toString());
			else write.write(output.toString());
		}
	}
}
