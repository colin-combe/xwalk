package cX;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.*;

/*
 * Created on Jun 24, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/**
 * @author abdullah
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
public class Grid{	
	/*--------------------------------------------------------------------------*/
	double cellLength;	
	// cell diagonale is c=sqrt(a*a+b*b). for radius only half of cell diagonal needed
	double cellRadius;

	double maxX = -99999999;
	double maxY = -99999999;
	double maxZ = -99999999;

	double minX = 99999999;
	double minY = 99999999;
	double minZ = 99999999;

	double potInitiate = 0;
	
	int noOfxCells = -1;
	int noOfyCells = -1;
	int noOfzCells = -1;
	
	boolean[][][] isMolecule;
	boolean[][][] isMoleculeCopy;
	double[][][] gridPot;
	boolean[][][] alreadyPassed;
	int[][][] cluster;
	Vector<Vector<Integer>> occupied = new Vector();
	double[][][] gridXcoords;
	double[][][] gridYcoords;
	double[][][] gridZcoords;
	double[][][] distances;
	
	private int cellHits  = 0;
	public int clusterNo = 0;
	private double templateRadius;

	static String nl = System.getProperty("line.separator");

	/*--------------------------------------------------------------------------*/
	public Grid(){
		
	}
	/*--------------------------------------------------------------------------*/
	
	public Grid(StructureCoords coords, double cellLength, int offSetNoCells, double offSetMinMax){

		this.cellLength = cellLength;
		
		// cell diagonale is c=sqrt(a*a+b*b). for radius only half of cell diagonal needed
		cellRadius = Math.sqrt(2*(cellLength*cellLength))/2;	

		this.calcGrid(coords, offSetNoCells, offSetMinMax);
//		this.cluster();
	}
				
	/*--------------------------------------------------------------------------*/

	public double getCellLength(){
		return cellLength;
	}
			
	/*--------------------------------------------------------------------------*/
	private void calcMinMax(StructureCoords coords, double offset){
		for(Iterator<StructureAtom> i=coords.iterator();i.hasNext();){
			StructureAtom atom = i.next();
			double x = atom.getX();
			double y = atom.getY();
			double z = atom.getZ();
			double r = atom.getRadius();

			this.maxX = Math.max(this.maxX, x+r+offset);
			this.maxY = Math.max(this.maxY, y+r+offset);
			this.maxZ = Math.max(this.maxZ, z+r+offset);

			this.minX = Math.min(this.minX, x-r-offset);
			this.minY = Math.min(this.minY, y-r-offset);
			this.minZ = Math.min(this.minZ, z-r-offset);
		}
	}

	/*--------------------------------------------------------------------------*/
	private void calcNoOfCells(int offset){
		this.noOfxCells = (int)((this.maxX-this.minX)/this.cellLength)+offset;
		this.noOfyCells = (int)((this.maxY-this.minY)/this.cellLength)+offset;
		this.noOfzCells = (int)((this.maxZ-this.minZ)/this.cellLength)+offset;
	}
	/*--------------------------------------------------------------------------*/
	public void initiate(){
		this.isMolecule 	= new boolean[this.noOfxCells][this.noOfyCells][this.noOfzCells];
		this.isMoleculeCopy	= new boolean[this.noOfxCells][this.noOfyCells][this.noOfzCells];
		this.gridXcoords 	= new double[this.noOfxCells][this.noOfyCells][this.noOfzCells];
		this.gridYcoords 	= new double[this.noOfxCells][this.noOfyCells][this.noOfzCells];
		this.gridZcoords 	= new double[this.noOfxCells][this.noOfyCells][this.noOfzCells];

		this.distances		= new double[this.noOfxCells][this.noOfyCells][this.noOfzCells];

		this.gridPot 		= new double[this.noOfxCells][this.noOfyCells][this.noOfzCells];

		this.alreadyPassed = new boolean[this.noOfxCells][this.noOfyCells][this.noOfzCells];
		this.cluster = new int[this.noOfxCells][this.noOfyCells][this.noOfzCells];

		
		for(int i=0; i<noOfxCells; i++){
			for(int j=0; j<noOfyCells; j++){
				for(int k=0; k<noOfzCells; k++){
					double x = minX + i*cellLength + cellLength/2;
					double y = minY + j*cellLength + cellLength/2;
					double z = minZ + k*cellLength + cellLength/2;
					this.gridXcoords[i][j][k] 	= x; 
					this.gridYcoords[i][j][k] 	= y; 
					this.gridZcoords[i][j][k] 	= z;
					this.distances[i][j][k] 	= 999.99;
					this.gridPot[i][j][k] 	= this.potInitiate; 
				}
			}
		}
	}

	/*--------------------------------------------------------------------------*/
	public double[] calcGridNo(StructureAtom atom) {

		double[] values = new double[4];
		
		double radius = atom.getRadius();
		// if min values are negative, gridNox gets pos
		values[0] = (int)((atom.getX()-this.minX)/this.cellLength)+0.0001;
		values[1] = (int)((atom.getY()-this.minY)/this.cellLength)+0.0001;
		values[2] = (int)((atom.getZ()-this.minZ)/this.cellLength)+0.0001;

		if(radius == 0.0)
			values[3] = 0.01;
		else
			values[3] = radius/this.cellLength;
		
		return values;
	}
	/*--------------------------------------------------------------------------*/
	public boolean isInGrid(StructureAtom atom, double threshhold) {
		
		int hit=0; 
		int nonHit=0; 
		int count=0;
		double[] atomGridNo = this.calcGridNo(atom);
		int ratio = atomGridNo[3] - (int)atomGridNo[3] - 0.001 > 0 ? (int)atomGridNo[3]+1: (int)atomGridNo[3];

		
		for(int i=(int)atomGridNo[0]-ratio;i<=atomGridNo[0]+ratio;i++){
			for(int j=(int)atomGridNo[1]-ratio;j<=atomGridNo[1]+ratio;j++){
				for(int k=(int)atomGridNo[2]-ratio;k<=atomGridNo[2]+ratio;k++){
					count++;	
					if(k<0 || j<0 || i<0 || k>this.noOfzCells-1 || j>this.noOfyCells-1 || i>this.noOfxCells-1 || 
					   !this.isMolecule[i][j][k]){
						nonHit++;
					}
					else{
						hit++;
					}
				}
			}
		}
		if(hit+nonHit!=count)System.err.print("WARNING: hit ("+hit+") nonHit ("+nonHit+") != count("+count+")"+nl);
		
		return (hit/count)>=threshhold ? true : false;			
	}
	/*--------------------------------------------------------------------------*/
	public boolean isInGrid2(StructureAtom atom) {
		double[] atomGridNo = this.calcGridNo(atom);
		for(int i=(int)(atomGridNo[0]-atomGridNo[3]);i<=atomGridNo[0]+atomGridNo[3];i++){
			for(int j=(int)(atomGridNo[1]-atomGridNo[3]);j<=atomGridNo[1]+atomGridNo[3];j++){
				for(int k=(int)(atomGridNo[2]-atomGridNo[3]);k<=atomGridNo[2]+atomGridNo[3];k++){
					
					if(k>0 && j>0 && i>0 && k<this.noOfzCells && j<this.noOfyCells && i<this.noOfxCells && 
					   this.isMolecule[i][j][k]){
						return true;
					}
				}
			}
		}
		return false;
	}
	/*--------------------------------------------------------------------------*/
	private void calcGrid(StructureCoords coords, int offSetNoCells, double offSetMinMax) {
		
		this.calcMinMax(coords, offSetMinMax);
		this.calcNoOfCells(offSetNoCells);
		this.initiate();
						
		for(int i=0; i<coords.size(); i++){
			StructureAtom atom = coords.getAtom(i);
			double[] atomGridNo = this.calcGridNo(atom);
			int ratio = atomGridNo[3] - (int)atomGridNo[3] - 0.001 > 0 ? (int)atomGridNo[3]+1: (int)atomGridNo[3];

			// apply grid mask
			boolean hit = false;
			double maxDist = 0;
			
			// sometimes it happens that a point has no grid point assigned as the next grid point is too far away. To allow assignments for those points
			// increase maximal allowed distance.
			while(!hit){
				for(int j=(int)atomGridNo[0]-ratio;j<=atomGridNo[0]+ratio;j++){
					for(int k=(int)atomGridNo[1]-ratio;k<=atomGridNo[1]+ratio;k++){
						for(int l=(int)atomGridNo[2]-ratio;l<=atomGridNo[2]+ratio;l++){

							double x = this.gridXcoords[j][k][l];
							double y = this.gridYcoords[j][k][l];
							double z = this.gridZcoords[j][k][l];
							double dist = Math.sqrt(Math.pow(atom.getX()-x,2) + 
									                Math.pow(atom.getY()-y,2) + 
											        Math.pow(atom.getZ()-z,2) );
							dist = dist-atom.getRadius()-cellRadius;
								
							if(dist <= 0){
								if(!this.isMolecule[j][k][l]){
									this.isMolecule[j][k][l] = true;
									this.isMoleculeCopy[j][k][l] = true;
									this.gridPot[j][k][l] = atom.getPotential();
									cellHits++;
						 			Vector<Integer> tmp = new Vector();
									tmp.add(j);
									tmp.add(k);
									tmp.add(l);
									occupied.add(tmp);									
								}
							}
							hit=true;
						}
					}
				}
				maxDist += 0.05;
			}
		}
	}
	//--------------------------------------------------------------------------------
	private void cluster(){

		for(Iterator<Vector<Integer>> i=occupied.iterator(); i.hasNext();){
			Vector<Integer> v = i.next();
			int j = v.get(0);
			int k = v.get(1);
			int l = v.get(2);
			
			if(!this.alreadyPassed[j][k][l]){
				this.clusterNo++;
				this.cluster[j][k][l] = this.clusterNo;
				this.alreadyPassed[j][k][l] = true;
				this.search4neighbours(j, k, l);
			}
		}
	}
	//--------------------------------------------------------------------------------
	private void search4neighbours(int x, int y, int z){

		for(		int m=x-1;m>=0 && m<=x+1 && m<this.noOfxCells;m++){
			for(	int n=y-1;n>=0 && n<=y+1 && n<this.noOfyCells;n++){
				for(int o=z-1;o>=0 && o<=z+1 && o<this.noOfzCells;o++){

					if(this.isMolecule[m][n][o]){
						if(!this.alreadyPassed[m][n][o]){
							this.cluster[m][n][o] = this.clusterNo;
							this.alreadyPassed[m][n][o] = true;
							this.search4neighbours(m, n, o);
						}
					}
				}
			}
		}	
	}
	
	//--------------------------------------------------------------------------------
	public void reset4DistCalc(){
		
		for(int i=0; i<noOfxCells; i++){
			for(int j=0; j<noOfyCells; j++){
				for(int k=0; k<noOfzCells; k++){
					this.distances[i][j][k] = 999.99;						
					this.alreadyPassed[i][j][k] = false;
					this.isMolecule[i][j][k] = isMoleculeCopy[i][j][k];
				}
			}
		}
	}
	//--------------------------------------------------------------------------------
	public double getSolventAccessRatio(StructureCoords coords, double templateRadius){

		int countSAS = 0;
		int countTotalNNgridPoints = 0;
	
		double sas = 0;
		int count = 0;
		double tmp = templateRadius/this.cellLength;
		for(Iterator<StructureAtom> i=coords.iterator(); i.hasNext();){
			StructureAtom atom = i.next();
			double[] no = this.calcGridNo(atom);
			int j = (int)no[0];
			int k = (int)no[1];
			int l = (int)no[2];

			int increment = no[3] - (int)no[3] - 0.001 > 0 ? (int)no[3]+1: (int)no[3];			
			if(templateRadius >= 0) increment = tmp - (int)tmp - 0.001 > 0 ? (int)tmp+1: (int)tmp; 

			sas+=this.getSolventAccessRatio(j, k, l, increment);
			count++;
		}
		return sas/count;
	}
	//--------------------------------------------------------------------------------
	private double getSolventAccessRatio(int j, int k, int l, int increment){

		int countSAS = 0;
		int countTotalNNgridPoints = 0;

		for(int m=-increment; m<=2*increment && j+m>=0 && j+m<this.noOfxCells; m+=2*increment){
			for(int n=-increment; n<=2*increment && k+n>=0 && k+n<this.noOfyCells; n+=2*increment){
				for(int o=-increment; o<=2*increment && l+o>=0 && l+o<this.noOfzCells; o+=2*increment){
					// gridStatus is set true if grid cell is occupied by protein
					// gridStatus is set false if grid cell is occupied by solvent
					if(!this.isMolecule[j+m][k+n][l+o]){
						countSAS++;
					}
					countTotalNNgridPoints++;
				}
			}
		}
		return countSAS/(double)countTotalNNgridPoints;
	}
		
	//--------------------------------------------------------------------------------
	public void treatCoordAsSolvent(Vector<StructureResidue> coords){
		
		for(Iterator<StructureResidue> i=coords.iterator(); i.hasNext();) 
			this.treatCoordAsSolvent(i.next());
	}
	//--------------------------------------------------------------------------------
	public void treatCoordAsSolvent(StructureCoords coords){
		for(Iterator<StructureAtom> i=coords.iterator(); i.hasNext();){
			StructureAtom atom = i.next();
			double[] no = this.calcGridNo(atom);
			int j = (int)no[0];
			int k = (int)no[1];
			int l = (int)no[2];

			int ratio = no[3] - (int)no[3] - 0.001 > 0 ? (int)no[3]+1 : (int)no[3];

			for(int m=-ratio; m<=ratio && j+m>=0 && j+m<this.noOfxCells; m++){
				for(int n=-ratio; n<=ratio && k+n>=0 && k+n<this.noOfyCells; n++){
					for(int o=-ratio; o<=ratio && l+o>=0 && l+o<this.noOfzCells; o++){
						this.distances[j+m][k+n][l+o]=999.99;
						this.alreadyPassed[j+m][k+n][l+o] = false;
						this.isMolecule[j+m][k+n][l+o] = false;
					}
				}
			}
		}
	}

	//--------------------------------------------------------------------------------
	public void setSurfaceDistances(StructureCoords coords, double maxDist, double templateRadius, double minSasRatio, boolean verbose){

		for(Iterator<StructureAtom> i=coords.iterator(); i.hasNext();){
			StructureAtom atom = i.next();
			double xyz1[] = {atom.getX(), atom.getY(), atom.getZ()};
			
			double[] no = this.calcGridNo(atom);
			int j = (int)no[0];
			int k = (int)no[1];
			int l = (int)no[2];
			int ratio = no[3] - (int)no[3] - 0.001 > 0 ? (int)no[3]+1 : (int)no[3];
			
			// initiate distance calculation with distances assigned to the grid points closest to the respective atom
			if(ratio>=1){
				for(int m=-ratio; m<=ratio && j+m>=0 && j+m<this.noOfxCells; m++){
					for(int n=-ratio; n<=ratio && k+n>=0 && k+n<this.noOfyCells; n++){
						for(int o=-ratio; o<=ratio && l+o>=0 && l+o<this.noOfzCells; o++){
							if(!this.isMolecule[j+m][k+n][l+o]){
								double xyz2[] = {this.gridXcoords[j+m][k+n][l+o], this.gridYcoords[j+m][k+n][l+o], this.gridZcoords[j+m][k+n][l+o]};
								double dist = Mathematics.distance(xyz1, xyz2);
								this.distances[j+m][k+n][l+o]=dist;
							}
						}
					}
				}
			}
			else{
				if(!this.isMolecule[j][k][l]){
					double xyz2[] = {this.gridXcoords[j][k][l], this.gridYcoords[j][k][l], this.gridZcoords[j][k][l]};
					double dist = Mathematics.distance(xyz1, xyz2);
					this.distances[j][k][l] = dist;
				}
			}
		}

		// start distance calculation, after initiation. Thereby not important which atom to take as a starting point as all grid cells associated to coords are set to 0.
		double[] no = this.calcGridNo(coords.getAtom(0));
		int j = (int)no[0];
		int k = (int)no[1];
		int l = (int)no[2];
		this.assignDistance2neighbours(j, k, l, maxDist, templateRadius, minSasRatio, verbose);
	}
		
	//--------------------------------------------------------------------------------
	// calculate the distance between grid point jkl and each grid point within maxDist distance by using the onion shell approach
	private void assignDistance2neighbours(int j, int k, int l, double maxDist, double templateRadius, double minSasRatio, boolean verbose){

		// number of shell to the calculate distance
		double ratio = maxDist/this.cellLength;
		int gridSize = ratio - (int)ratio - 0.001 > 0 ? (int)ratio+1 : (int)ratio; 

		double tmp = templateRadius/this.cellLength;
		int increment = tmp - (int)tmp - 0.001 > 0 ? (int)tmp+1: (int)tmp;
		
		for(int i=0; i<=gridSize; i++){
			for(		int m=j-i;m>=0 && m<=j+i && m<this.noOfxCells;m++){
				for(	int n=k-i;n>=0 && n<=k+i && n<this.noOfyCells;n++){
					for(int o=l-i;o>=0 && o<=l+i && o<this.noOfzCells;o++){
						if(!this.alreadyPassed[m][n][o]){
							if(minSasRatio > 0) if(this.getSolventAccessRatio(m, n, o, increment) < minSasRatio) continue;
							this.alreadyPassed[m][n][o] = true;
							this.calculateDistances(m, n, o);
							if(verbose){
								StructureAtom atom = new StructureAtom();
								atom.setName("C");
								atom.setX(this.gridXcoords[m][n][o]);
								atom.setY(this.gridYcoords[m][n][o]);
								atom.setZ(this.gridZcoords[m][n][o]);
								double dist = this.distances[m][n][o];
								atom.setTemp(dist);
								System.out.print(atom.toString());
							}
						}
					}
				}
			}
		}
		if(verbose)
			System.out.print("END"+nl);
	}
	//--------------------------------------------------------------------------------
	private void calculateDistances(int j, int k, int l){
		double parentX = this.gridXcoords[j][k][l];
		double parentY = this.gridYcoords[j][k][l];
		double parentZ = this.gridZcoords[j][k][l];
		double[] parentXYZ = {parentX,parentY,parentZ};
		for(		int m=j-1;m>=0 && m<=j+1 && m<this.noOfxCells;m++){
			for(	int n=k-1;n>=0 && n<=k+1 && n<this.noOfyCells;n++){
				for(int o=l-1;o>=0 && o<=l+1 && o<this.noOfzCells;o++){

					double childX = this.gridXcoords[m][n][o];
					double childY = this.gridYcoords[m][n][o];
					double childZ = this.gridZcoords[m][n][o];
					double[] childXYZ = {childX,childY,childZ};
					
					double distParentChild = Mathematics.distance(parentXYZ, childXYZ);
					double distParent 	   = this.distances[j][k][l];
					double distChild = this.distances[m][n][o];
					double distCurrent 	   = distChild+distParentChild;
					// assign distances only to non-atom containing grid points
					if( distCurrent < distParent){
						this.distances[j][k][l] = distCurrent;
					}
				}
			}
		}
	}
	
	//--------------------------------------------------------------------------------
	public double getSurfaceDistances(StructureCoords coords2){

		double minDist = 999.99;
		for(Iterator<StructureAtom> i=coords2.iterator(); i.hasNext();){
			StructureAtom atom = i.next();
			double[] no = this.calcGridNo(atom);
			int j = (int)no[0];
			int k = (int)no[1];
			int l = (int)no[2];			

//			this.proceedWithDistanceAssignment(j, k, l, templateRadius);
			double xyz1[] = {atom.getX(),atom.getY(),atom.getZ()};
			double xyz2[] = {this.gridXcoords[j][k][l],this.gridYcoords[j][k][l],this.gridZcoords[j][k][l]};
			double dist1 = Mathematics.distance(xyz1, xyz2);
			double dist = this.distances[j][k][l]-dist1;
			if(minDist>dist){
				minDist = this.distances[j][k][l]-dist1;
			}
		}
		return minDist;
	}
	/*--------------------------------------------------------------------------*/
	public String toString(){
		StringBuffer output = new StringBuffer();
		for(int i=0; i<noOfxCells; i++){
			for(int j=0; j<noOfyCells; j++){
//				output += i+""+j+"\t";
				for(int k=0; k<noOfzCells; k++){
					if(isMolecule[i][j][k]){
						StructureAtom atom = new StructureAtom();
						atom.setX(gridXcoords[i][j][k]);
						atom.setY(gridYcoords[i][j][k]);
						atom.setZ(gridZcoords[i][j][k]);
						atom.setChainId('Y');
						output.append(atom.toString(atom.getOccupancy(), gridPot[i][j][k]));
//						output +="1\t";						
					}
					else{
						StructureAtom atom = new StructureAtom();
						atom.setX(gridXcoords[i][j][k]);
						atom.setY(gridYcoords[i][j][k]);
						atom.setZ(gridZcoords[i][j][k]);
						atom.setChainId('N');
						output.append(atom.toString(atom.getOccupancy(), gridPot[i][j][k]));
//						output +="0\t";						
					}
				}
//				output +="\n";
			}
//			output +="------------------------------------------------------------------------------------------------------------------------------------------------------\n";
		}
	return output.toString();
	}
	
	/*--------------------------------------------------------------------------*/

	public static void main(String[] args){
		String alphabet			= " ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
		Commandline commandline = new Commandline();
		String infile = "";
		//----------------------------
		// reading in from commandline "in="
		if (!commandline.get(args,"-in", true).equals("ERROR")){
			infile = commandline.get(args,"-in", true);
			try{
				FileReader fileReader = new FileReader(infile);
			}
			catch(FileNotFoundException e){
				System.err.print("File \""+infile+"\" not found!!!\n\n");
				System.exit(1);
			}
		}

		StructureCoords coords = new StructureCoords(infile);
		coords.extractAtoms("ATOM  #HETATM",alphabet, alphabet, false, false);

		Grid grid = new Grid(coords, 3, 0, 18);
//		grid.setSurfaceDistances(coords, 15, false);
		
//		System.out.print(grid.toString());
//		System.err.print(grid.minX+" "+grid.minY+" "+grid.minZ+"\n"+grid.maxX+" "+grid.maxY+" "+grid.maxZ+"\n"+grid.noOfxCells+" "+grid.noOfyCells+" "+grid.noOfzCells+"\n");
		
//		DxFile dx = new DxFile();
//		System.out.print(dx.toString(grid));
//		System.err.print("2\n");
		
//		Clefts clefts = new Clefts(coords, true, false);
//		System.out.print(clefts+"asdfadsf\n");
	}
	
} // end of class Grid
