package cX;

import java.util.Vector;
import java.util.Collections;

import java.text.DecimalFormat;
import java.util.Locale;



/**
 * Ueberschrift:  PropSearch
 * Beschreibung:  Statistics. Class with statistic operations and functions.
 * Copyright:     Copyright (c) 2003
 * Organisation:  FH-Giessen
 * @author Abdullah Kahraman
 * @version 1.0
 */

public class Mathematics {

	/*--------------------------------------------------------------------------*/
	// Classmethods
	/*--------------------------------------------------------------------------*/    

	static public double electronCharge 	= 1.6021765314E-19; 	// in C
	static public double electronMass 		= 9.109382616E-28;	 	// in g
	static public double diracsConstant 	= 1.05457159682E-34; // in J s//6.5821191556E-16; 	// eV s
	
	
	//----------------------------------------------------------------------------//
	static public double sigmoidPotential(double dist, double hydrophobEnvRadius){
		double relDist = dist/hydrophobEnvRadius;
		if(relDist<=1){
			return 1-0.5*(7*Math.pow(relDist,2)-9*Math.pow(relDist,4)+5*Math.pow(relDist,6)-Math.pow(relDist,8));
		}
		else{
			return 0;
		}
	}
	
	/*--------------------------------------------------------------------------*/        
	/** calculate avg value from several double values out of a vector object.
	 * @param Vector values: Double values from which the avg value has to be build.
	 * @return double: return the avg value from 'values'.
	 * @throws ClassCastException: if the values in the Vector are not Double object
	 * method will throw a ClassCastException.
	 */
	static public double avg(Vector values) throws ClassCastException{

		double avg = 0.0;

		for(int i=0;i<values.size();i++){
			avg += ((Double)values.get(i)).doubleValue();
		}
		avg = avg/values.size();
		
	return avg;
	} // End of method avg()

	/*--------------------------------------------------------------------------*/        
	static public double max(Vector values) throws ClassCastException{

		double max = ((Double)values.get(0)).doubleValue();
		for(int i=0;i<values.size()-1;i++){
			double value = ((Double)values.get(i)).doubleValue();
		 	max = Math.max(max, value);
		}
		
	return max;
	} // End of method max()

	/*--------------------------------------------------------------------------*/        
	static public double max(double[] values) throws ClassCastException{

		double max = values[0];
		for(int i=0;i<values.length-1;i++){
			double value = values[i];
			max = Math.max(max, value);
		}
		
	return max;
	} // End of method max()

	/*--------------------------------------------------------------------------*/        
	static public int max(int[] values) throws ClassCastException{

		int max = values[0];
		for(int i=0;i<values.length-1;i++){
			int value = values[i];
			max = Math.max(max, value);
		}
		
	return max;
	} // End of method max()

	/*--------------------------------------------------------------------------*/        

	/** calculate the spherical coordinates out of a point with x,y,z coordinates
	 * @param double x - x coordinate of a point
	 * @param double y - y coordinate of a point
	 * @param double z - z coordinate of a point
	 * @return double[]: return a array with double values:
	 * 					 first is the radius, 
	 * 					 second is degree phi,
	 * 					 third is degree theta
	 */
	static public double[] xyz2rtp(double x, double y, double z) {

		double r = Math.sqrt(Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2));
		double theta = 0;
		double phi = 0;
		
		double eps = 0.0000000001; // zero in double format is never exactly 0.0
		
		if (Math.abs(r) != 0.0) {
			theta = Math.acos(z/r);	

			if (    Math.abs(x) > eps && Math.abs(y) > eps) {
				phi = Math.atan(y/x);
			}
			else if(Math.abs(x) > eps && Math.abs(y) < eps){
				phi = 0;
			}
			else if(Math.abs(x) < eps && y < 0){
				phi = -Math.PI/2;
			}
			else if(Math.abs(x) < eps && y > 0){
				phi = Math.PI/2;
			}

		}
		else {
			// doesn't matter which value theta has, because there is no degree necessary
			// when there is no vector (caused of radius = 0)
			theta = 0;
			phi = 0;
		}		

		/* phi describes the angle in the xy-plane from the x-axis, whereas theta 
		 * describes the angle from the z-axis to the radius vector.
		 * The equation above for phi works, but has the problem that arctan can not 
		 * distinguish between y/x and -y/-x or -y/x and y/-x. In both cases the function
		 * will return the value for y/x and -y/x which conform to values between -PI/2 and 
		 * PI/2. With other words, we get only values in the first and fourth coordinate 
		 * area. To get values in the second and third coordinate, which would conform to
		 * a negative x value (x < 0) we had to sum a PI to the phi value.
		 * Furthermore the phi value should be used for the recalculation of the x and y,
		 * which caused that we had to sum 2*PI to every negative phi value, to get a 
		 * identical positiv phi value.
		 */

		if(x < 0.0 && Math.abs(x) > eps) {
			phi = phi + Math.PI; 
		}
		if (phi < 0.0 && Math.abs(phi) > eps) {
			phi = phi + 2*Math.PI;
		}
		
		// testing with recalculating x, y, z
		double[] xyz = new double[3];
		xyz[0] = r * Math.cos(phi) * Math.sin(theta);
		xyz[1] = r * Math.sin(phi) * Math.sin(theta);
		xyz[2] = r * Math.cos(theta);
		
		double x1 = xyz[0];
		double y1 = xyz[1];
		double z1 = xyz[2];
		
		if ( Math.abs(x1-x) > 0.001 ) {
			System.err.print("x:"+x+" "+x1+"\tinconsistence of value x\n"+
							 "y:"+y+" "+y1+"\n"+"z:"+z+" "+z1+"\n"+
							 "radius:"+r+"\t"+"theta:"+theta+"\t"+"phi:"+phi+"\n");
//			System.exit(1);
		}
		
		else if ( Math.abs(y1-y) > 0.001 ) {
			System.err.print("y:"+y+" "+y1+"\tinconsistence of value y\n"+
							 "x:"+x+" "+x1+"\n"+"z:"+z+" "+z1+"\n"+
							 "radius:"+r+"\t"+"theta:"+theta+"\t"+"phi:"+phi+"\n");
//			System.exit(1);
		}

		if ( Math.abs(z1-z) > 0.001 ) {
			System.err.print("z:"+z+":"+z1+"\tinconsistence of value z\n"+
							 "x:"+x+" "+x1+"\n"+"y:"+y+" "+y1+"\n"+
							 "radius:"+r+"\t"+"theta:"+theta+"\t"+"phi:"+phi+"\n");
//			System.exit(1);
		}

		double[] sphericalCoord = new double[3];
		sphericalCoord[0] = r;
		sphericalCoord[1] = theta;
		sphericalCoord[2] = phi;
			
	return sphericalCoord;
	} // End of method xyz2rpt()

	/*--------------------------------------------------------------------------*/        
	/** calculate the cartesian coordinates out of a spherical coordinates 
	 *   r, phi, theta
	 * @param double r - radius
	 * @param double phi - angle in the xy-plane
	 * @param double theta - angle in the z-axis
	 * @return double[]: return a array with double values:
	 * 					 first is the x, 
	 * 					 second is y,
	 * 					 third is z
	 */
	static public double[] rtp2xyz(double r, double theta, double phi) {

		double eps = 0.0000000001; // zero in double format is never exactly 0.0

		double x = r * Math.cos(phi) * Math.sin(theta);
		double y = r * Math.sin(phi) * Math.sin(theta);
		double z = r *                 Math.cos(theta);
		
//System.out.print(cartCoord[0]+","+cartCoord[1]+","+cartCoord[2]+"\n");
		// testing
		double[] rtp = xyz2rtp(x,y,z);
		double r1 = rtp[0];
		double t1 = rtp[1];
		double p1 = rtp[2];
				
		if ( Math.abs(r1-r) > 0.01) {
			System.err.print("r:"+r+" "+r1+"\tinconsistence of value r\n"+
							 "theta:"+theta+" "+t1+"\n"+"phi:"+phi+" "+p1+"\n"+
				 			 "x:"+x+"\n"+"y:"+y+"\n"+"z:"+z+"\n");
//			System.exit(1);
		}
		
		
		if ( Math.abs(t1-theta) > 0.01 && Math.abs(r) > eps) {
			System.err.print("theta:"+theta+" "+t1+"\tinconsistence of value theta\n"+
							 "phi:"+phi+" "+p1+"\n"+"r:"+r+" "+r1+"\n"+
							 "x:"+x+"\n"+"y:"+y+"\n"+"z:"+z+"\n");
//			System.exit(1);
		}

		if ( Math.abs(p1-phi) > 0.01 && Math.abs(theta) > eps && Math.abs(r) > eps) {
			System.err.print("phi:"+phi+" "+p1+"\tinconsistence of value phi\n"+
							 "theta:"+theta+" "+t1+"\n"+"r:"+r+" "+r1+"\n"+
							 "x:"+x+"\n"+"y:"+y+"\n"+"z:"+z+"\n");
//			System.exit(1);
		}

		double[] cartCoord = new double[3];
		// when calculating XYZ values should be three digits after komma. If no
		// rounded up, calculations can be unprecise.
		Locale.setDefault(Locale.US);
		DecimalFormat decFormat = new DecimalFormat("0.000");
		cartCoord[0] = Double.parseDouble(decFormat.format(x));
		cartCoord[1] = Double.parseDouble(decFormat.format(y));
		cartCoord[2] = Double.parseDouble(decFormat.format(z));

	return cartCoord;
	} // End of method xyz2rpt()

	/*--------------------------------------------------------------------------*/        
	/** calculate deviation from several double values out of a vector object.
	 * @param Vector values: Double values from which the deviation value has to build.
	 * @return double: return the deviation value from 'values'.
	 * @throws ClassCastException: if the values in the Vector are not Double object
	 * method will throw a ClassCastException.
mo	 */
	static public double deviation(Vector values) throws ClassCastException{
		
		double avg = 0.0;
		double div = 0.0;
		
		avg = avg(values);
		
		for(int i=0;i<values.size();i++){
			div += Math.pow(((Double)values.get(i)).doubleValue() - avg,2);
		}
		if(values.size() == 1){
			div = 0;
		}
		else {
			div = Math.sqrt(div/(values.size()-1));
		}
	return div;
	} // End of method deviation()

	/*--------------------------------------------------------------------------*/        
	/** calculate the median value from several double values out of a vector object.
	 * @param Vector values: Double values from which the median value has to build.
	 * @return double: return the median value from 'values'.
	 * @throws ClassCastException: if the values in the Vector are not Double object
	 * method will throw a ClassCastException.
	 */
	static public double median(Vector values) throws ClassCastException{
		
		double median = 0.0;
		
		Collections.sort(values);

		// the median is the value at the middle of a sorted list
		int arraySize = values.size();
		if(arraySize%2 == 0){
		   median = (((Double)values.get((arraySize/2)-1)).doubleValue() + 
		   ((Double)values.get(arraySize/2)).doubleValue() ) / 2;
		}
		else{
		   median = ((Double)values.get(arraySize/2)).doubleValue();
		}

	return median;
	} // End of method median()

	/*--------------------------------------------------------------------------*/        
	/** calculates the ttest, does not assume equal variance in both samples.
	 * @param double ave1: avg value of set1.
	 * @param double ave2: avg value of set2.
	 * @param double sd1: deviation value of set1.
	 * @param double sd2: deviation value of set2.
	 * @param int numberOfFiles1: count of values from which the avg/deviation value was build.
	 * @param int numberOfFiles2: count of values from which the avg/deviation value was build.
	 * @return double: t-test Value or pValue (input NULL if not interested in result).
	 * @return in case of problems (degrees of freedom = 0 or variance = 0, returns -1. 
	*/
	public double tTest(double avg1, double avg2, double sd1, double sd2, int numberOfFiles1, int numberOfFiles2) {
	  double df;
	  double v1,v2;
	  double t1,d;
	  double prob = -1.0;
	  double t = -1.0;
          
	  v1 = sd1*sd1;
	  v2 = sd2*sd2;
	  d = Math.sqrt (v1/numberOfFiles1 + v2/numberOfFiles2);
          
	  // d is only null if both std deviation are null. 
	  // This is possible only if a single value, namely avg is given.
	  // Without std dev no t-test can be performed, thus accept null-hypothesis, i.e. -1.
	  if (d == 0.0) {
					if ( t != 0.0 )
							t = -1.0;
					if ( prob != 0.0)
							prob = -1.0;
					return (prob);
	  }
          
	  t1 = (avg1-avg2)/d;
          
	  if (t != 0)
					t = t1;
	  if (prob != 0) {
					// for calc. standard deviation one requires at least two numbers. Thus if only one number given accept null hypothesis, i.e. -1
					if (numberOfFiles1 < 2 || numberOfFiles2 < 2) {
						prob = -1.0;
						return prob;
					}
		  
					/* formula used by Excel */
					df = Math.pow (v1/numberOfFiles1 + v2/numberOfFiles2, 2) /
					 ((v1*v1/(numberOfFiles1*numberOfFiles1))/(numberOfFiles1-1) + (v2*v2/(numberOfFiles2*numberOfFiles2))/(numberOfFiles2-1));

					/* formula from Lothar Sachs, Statistische Auswertungsmethoden
					df = pow (v1/numberOfFiles 1 + v2/numberOfFiles 2,2) /
					 ((v1*v1/(numberOfFiles 1*numberOfFiles 1))/(numberOfFiles 1+1) + (v2*v2/(numberOfFiles 2*numberOfFiles 2))/(numberOfFiles 2+1)) - 2;
					*/
			prob = rcp_betai (0.5*df,0.5,df/(df+t1*t1));
	  }
	  return prob;
	} // End of method ttest

	/*--------------------------------------------------------------------------*/        

	private double rcp_gammln (double xx)
	{
	  double x,y,tmp,ser;
	  int j;
	  double cof[] = {76.18009172947146,-86.50532032941677,24.01409824083091,
									  -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	  y = x = xx;
	  tmp = x+5.5;
	  tmp -= (x+0.5)*Math.log (tmp);
	  ser = 1.000000000190015;
	  for (j=0;j<6;j++) {
			ser += cof[j]/++y;
	  }
	  return -tmp + Math.log (2.5066282746310005*ser/x);
	} // End of method rcp_gammln()

	/*--------------------------------------------------------------------------*/        

	private double rcp_betai (double a,double b,double x) {
          
	  double bt;

	  if (x < 0.0 || x > 1.0)
			System.err.println("rcp_betai: x out of range");
	  if (x == 0.0 || x == 1.0)
			bt = 0.0;
	  else
			bt = Math.exp (rcp_gammln (a+b) - rcp_gammln (a) - rcp_gammln (b) +
							  a*Math.log (x) + b*Math.log (1.0-x));
	  if (x < (a+1.0)/(a+b+2.0))
			return bt*rcp_betacf (a,b,x)/a;
	  else
			return 1.0-bt*rcp_betacf (b,a,1.0-x)/b;
	}// End of method rcp_betai()

	/*--------------------------------------------------------------------------*/        

	private double rcp_betacf (double a,double b,double x)
	{
	  int itmax = 100;
	  double eps = 3.0e-7;
	  double fpmin = 1.0e-30;
	  double aa,c,d,del,h,qab,qam,qap;
	  int m,m2;

	  qab = a+b;
	  qap = a+1.0;
	  qam = a-1.0;
	  c = 1.0;
	  d = 1.0-qab*x/qap;
	  if (Math.abs (d) < fpmin)
			d = fpmin;
	  d = 1.0/d;
	  h = d;
	  for (m=1;m<=itmax;m++) {
			m2 = 2*m;
			aa = m*(b-m)*x/((qam+m2)*(a+m2));
			d = 1.0 + aa*d;
			if (Math.abs (d) < fpmin)
			  d = fpmin;
			c = 1.0 + aa/c;
			if (Math.abs (c) < fpmin)
			  c = fpmin;
			d = 1.0/d;
			h *= d*c;
			aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
			d = 1.0 + aa*d;
			if (Math.abs (d) < fpmin)
			  d = fpmin;
			c = 1.0 + aa/c;
			if (Math.abs (c) < fpmin)
			  c = fpmin;
			d = 1.0/d;
			del = d*c;
			h *= del;
			if (Math.abs (del-1.0) < eps)
			  break;
	  }
	  if (m > itmax)
			System.err.println("rcp_betacf: a or b too big, or itmax too small");

	  return h;
	}// End of method rcp_betacf()
	/*--------------------------------------------------------------------------*/        
	static public double distance(StructureAtom atom1, StructureAtom atom2) {
		
		double x1 = atom1.getX();
		double y1 = atom1.getY();
		double z1 = atom1.getZ();

		double x2 = atom2.getX();
		double y2 = atom2.getY();
		double z2 = atom2.getZ();
		
		double dist = Math.sqrt( Math.pow(x2-x1,2) + Math.pow(y2-y1,2) + Math.pow(z2-z1,2) );
	
//System.err.print("dist: "+dist+"\nradius1: "+atom1.getRadius()+"\nradius2:"+atom2.getRadius()+"\n");
		
		return dist;
	}
	
	/*--------------------------------------------------------------------------*/        
	static public double distance(double[] xyz1, double[] xyz2) {				
		return Math.sqrt( Math.pow(xyz2[0]-xyz1[0],2) + Math.pow(xyz2[1]-xyz1[1],2) + Math.pow(xyz2[2]-xyz1[2],2) );
	}
	
	/*--------------------------------------------------------------------------*/        
	static public double rad2dec (double radAngle){
		return (radAngle * 180/Math.PI);
	}
	/*--------------------------------------------------------------------------*/        
	static public double dec2rad (double decAngle){
		return (decAngle * Math.PI / 180);
	}
	
	/*--------------------------------------------------------------------------*/        
	// matrix taken from
	//    z
	//    |
	//    |
	//    /------y
	//   /
	//  x
	// phi: rotation around z axis
	// theta: rotation around x axis
	// teta: rotation around z axis (again)
	// http://mathworld.wolfram.com/EulerAngles.html
	static public double[][] getEulerRotationMatrix (double phi, double theta, double teta){
		double[][] a = new double[3][3];
		a[0][0] = Math.cos(teta)*Math.cos(phi)-Math.cos(theta)*Math.sin(phi)*Math.sin(teta);
		a[0][1] = Math.cos(teta)*Math.sin(phi)+Math.cos(theta)*Math.cos(phi)*Math.sin(teta);
		a[0][2] = Math.sin(teta)*Math.sin(theta);

		a[1][0] = -Math.sin(teta)*Math.cos(phi)-Math.cos(theta)*Math.sin(phi)*Math.cos(teta);
		a[1][1] = -Math.sin(teta)*Math.sin(phi)+Math.cos(theta)*Math.cos(phi)*Math.cos(teta);
		a[1][2] =  Math.cos(teta)*Math.sin(theta);

		a[2][0] = Math.sin(theta)*Math.sin(phi);
		a[2][1] = -Math.sin(theta)*Math.cos(phi);
		a[2][2] = Math.cos(theta);
		
		return a;
	}
	
	/*--------------------------------------------------------------------------*/        
	static public double length (double[] v){
		
	return Math.sqrt(
					 Math.pow(v[0],2)+
					 Math.pow(v[1],2)+
					 Math.pow(v[2],2)
					 );
	}

	/*--------------------------------------------------------------------------*/        
	static public double length (int[] i){
		
	return Math.sqrt(
					 Math.pow(i[0],2)+
					 Math.pow(i[1],2)+
					 Math.pow(i[2],2)
					 );
	}
	/*--------------------------------------------------------------------------*/        
	static public int[] subtraction (int[] vector1, int[] vector2){

		int[] sub = {vector2[0]-vector1[0],
						vector2[1]-vector1[1],
						vector2[2]-vector1[2]
						};

	return sub;
	}
	/*--------------------------------------------------------------------------*/        
	static public double[] subtraction (double[] vector1, double[] vector2){
		double[] sub = {	vector2[0]-vector1[0],
					 		vector2[1]-vector1[1],
					 		vector2[2]-vector1[2]};
	return sub;
	}
	/*--------------------------------------------------------------------------*/        
	static public double[] addition (double[] vector1, double[] vector2){

		double[] sub = {vector2[0]+vector1[0],
					 vector2[1]+vector1[1],
					 vector2[2]+vector1[2]
					};

	return sub;
	}
	/*--------------------------------------------------------------------------*/        
	static public double[] multFactor (double[] vector, double factor){

		double[] sub = {factor * vector[0],
						factor * vector[1],
						factor * vector[2]
						};

	return sub;
	}
	
	/*--------------------------------------------------------------------------*/        

	static public double[] negative(double[] vector){
		double[] neg = vector.clone();
		neg[0] = -vector[0];
		neg[1] = -vector[1];
		neg[2] = -vector[2];
		return neg;
	}

	//--------------------------------------------------------------------------//
	
	static public double[] mult(double[][] matrix, double[] vector){
		double[] result = new double[vector.length];

		for (int row=0; row<matrix.length; row++) {			
			if (matrix[row].length != vector.length)	
				throw new IllegalArgumentException("Matrix column numbers are not equal to vector row numbers\n");

			double s = 0;
			for (int column=0; column<matrix[row].length; column++ ) {
				s += matrix[row][column] * vector[column];
			}
			result[row] = s;
		}
		return result;

	}
	/*--------------------------------------------------------------------------*/        
/*
	static public void cartesian2internal(IMolecule mol){
		double[][] zMatrix = new double[mol.getAtomCount()][3];

		zMatrix[0][0] = 0;
		zMatrix[0][1] = 0;
		zMatrix[0][2] = 0;
		
		zMatrix[1][0] = Mathematics.distance(mol.getAtom(0), mol.getAtom(1));
		zMatrix[1][1] = 0;
		zMatrix[1][2] = 0;
		
		zMatrix[2][0] = Mathematics.distance(mol.getAtom(0), mol.getAtom(1));
		zMatrix[2][1] = 0;
		zMatrix[2][2] = 0;
		
		
		for(int i=0; i<mol.getAtomCount(); i++){
			IAtom atom = mol.getAtom(i);
			
		}
		
		
	}
	
*/
	/*--------------------------------------------------------------------------*/        

	/** main-method. output the deviation value of '1.0, 2.0, 3.0, 4.0'.
	 */ 
	public static void main(String[] args) {
		
/*		Vector v = new Vector(); 
		v.add(new Double(1.0));
		v.add(new Double(2.0));
		v.add(new Double(3.0));
		v.add(new Double(4.0));
*/
		double[] a = {13.782,54.683,118.968};
		double[] b = {13.438,54.083,119.693};
//		double[] c = {14.424,53.785,118.733};
		double[] c = {14.636,54.203,118.119};
		
		double[] v1 = new double[3];
		v1[0] = b[0]-a[0];
		v1[1] = b[1]-a[1];
		v1[2] = b[2]-a[2];

		double[] v2 = new double[3];
		v2[0] = c[0]-a[0];
		v2[1] = c[1]-a[1];
		v2[2] = c[2]-a[2];
		
		
	} // End of method main()

//	-----------------------------------------------------------------
} // End of class Statistics
