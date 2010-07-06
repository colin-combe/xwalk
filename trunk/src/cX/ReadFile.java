package cX;

import java.util.Vector;
import java.util.ArrayList;
import java.util.Iterator;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.File;

/**
 * Ueberschrift:  Bioinformatik
 * Beschreibung:  Bioinformatic programs. Class to read a file in a ArrayList object.
 * Copyright:     Copyright (c) 2003
 * Organisation:  FH-Giessen
 * @author Abdullah Kahraman
 * @version 1.0
 */
public class ReadFile extends ArrayList{

	/*--------------------------------------------------------------------------*/
	// Classobjects
	/*--------------------------------------------------------------------------*/
	
	private String fileName;
		  						
	/*--------------------------------------------------------------------------*/
	// Classmethods
	/*--------------------------------------------------------------------------*/    

	/**Konstruktor
	 * @param fileName: name of the file to read. If an exception produced during 
	 * creating a output channel from fileName, it will be catched from this method.
	 */
/*	public ReadFile(String fileName) {

		String lineBuffer;
			
		try{			
			BufferedReader br = new BufferedReader( 
				 	new InputStreamReader( 
			 			new FileInputStream(fileName))); // Dateistrom von der PDBFIND.txt Datei
			while ((lineBuffer = br.readLine()) != null ) {
					this.add(lineBuffer+"\n");
			}
			br.close();
		}
		catch(Exception e){
			System.err.print("###Error while opening file \"" + fileName + "\".###\n" +
							 "###"+e + "###\n");
		}
		this.fileName = fileName;
	} // End of constructor
*/
	//--------------------------------------------------------------------------------------------------
	/**Konstruktor
	 * @param fileName: name of the file to read. If an exception produced during 
	 * creating a output channel from fileName, it will be catched from this method.
	 */
	public ReadFile(BufferedReader br) {

		String lineBuffer;
			
		try{			
			while ((lineBuffer = br.readLine()) != null ) {
					this.add(lineBuffer+"\n");
			}
			br.close();
		}
		catch(Exception e){
			System.err.print("ERROR: while opening BufferedReader\n" +
							 e+"\n");
		}
	} // End of constructor

	//--------------------------------------------------------------------------------------------------
	
	public ReadFile(String fileName){
		try {
			
			String lineBuffer;

			InputStream ins = this.getClass().getClassLoader().getResourceAsStream(fileName);
   			BufferedReader br;
   			if(ins==null){
   		   		br = new BufferedReader (new FileReader (fileName));
   			}
   			else{
   				br = new BufferedReader(new InputStreamReader(ins));
   			}

			while ((lineBuffer = br.readLine()) != null ) {
				this.add(lineBuffer+"\n");
		}
		br.close();
		
	   }
		catch(Exception e){
			System.err.print("###Error while opening file \"" + fileName + "\" within CLASSPATH directory.###\n" +
							 "###"+e + "###\n");
		}
		this.fileName = fileName;
	}
	//--------------------------------------------------
	public static boolean exists(String fileName) {

		try{
			BufferedReader br = new BufferedReader( 
					new FileReader (fileName));
			String lb;
			while ((lb = br.readLine()) != null );
			br.close();

			return true;
		}
		catch(Exception e){

			return false;
		}
	} 
	//--------------------------------------------------
	public static boolean dirExists(String dirPath) {
	    boolean exists = (new File(dirPath)).exists();
	    if (exists) {
	    	return true;
	    } else {
	    	return false;
	    }
	} 


	//--------------------------------------------------
	/** Reads infile in a String object.
	 * @return String: String with the lines of infile.
	 */
	public String toString() {
		
		StringBuffer fileContent = new StringBuffer();

		for(Iterator i=this.iterator(); i.hasNext();){
			fileContent.append((String)i.next());
		}

	return fileContent.toString();
	} // End of method read()
	
	//--------------------------------------------------
	/** Read only a specific column out of infile.
	 * @param int columnNumber: number of column which should be read out [count beginning from 0].
	 * @return Vector: in this vector each entry is a value of readed column.   
	 */
	public Vector getColumn(int columnNumber, String delim){
		
		Vector column = new Vector();
	        	        						
		for(Iterator i=this.iterator(); i.hasNext();){
			String[] lineArray = ((String)(i.next())).split(delim);
			if (lineArray.length <= columnNumber) {
				System.err.print("###Can't find column no.\""+columnNumber+"\". "+
								 "Only \""+(lineArray.length-1)+"\" columns exists "+
								 "in file \""+fileName+"\"###\n");
 			}
			else{
				column.add(lineArray[columnNumber]);
			}
		}
	
	return column;
	} // End of method readColumn()
	
	//--------------------------------------------------
	public void rename(String newName){
		WriteFile file = new WriteFile(newName);
		file.write(this.toString());
		
		File file2delete = new File(fileName);
		file2delete.delete();
		
	}
	
	
	//--------------------------------------------------
	public void copy(String newName){
		WriteFile file = new WriteFile(newName);
		file.write(this.toString());
			
	}
	
	//-----------------------------------------------------------------
	/** main-method. Reads the parameters infile, and optional outfile out 
	 * of commandline and put the content of infile out.
	 */
	public static void main(String[] args) {

		String infile = new String();
		String outfile = new String();
		
		/*--------------------------------------------------------------------------*/
		// to read all needed parameters from commandline
		// @param : String[] args, list all of parameters
	
		Commandline commandline = new Commandline();
	
		/*-----------user information---------------*/
		if (args.length == 0){
			System.out.print("\njava ReadFile infile=ReadFile.java\n\n"+
								"Mehr Informationen gibt es unter \"java ReadFile -help\"\n\n");
			System.exit(0);
		}
		// if "-help" is written as parameter. a describtion about the parameter
		// will display
		if (args.length==1){
			if (args[0].equalsIgnoreCase("-help")){
				System.out.print("\njava ReadFile infile=ReadFile.java outfile=test.mrg\n"+
								 "\t-  infile=\t file to read\n"+
								 "\t-  outfile=\t if output to a file is desired[optional]\n");
				System.exit(0);
			 }
		}
	
		//----------------------------
		// bearbeiten von Parameter "infile="
		if (commandline.get(args,"infile=", true).equals("ERROR")){
			System.out.print("\nFehler beim Einlesen vom Parameter \"infile=\".\n");
			System.exit(1);
		}
		else{
			infile = commandline.get(args,"infile=",true);
			try{
				FileReader fileReader = new FileReader(infile);
			}
			catch(FileNotFoundException e){
				System.out.print("\nDie Datei \""+infile+"\" konnte nicht gefunden werden\n\n");
				System.exit(1);
			}
		}

		//----------------------------
		// search for parameter "outfile="
		outfile = commandline.get(args,"outfile=", true);

		//-----------End of analyzing the commandline parameter-------------------------------

		ReadFile readFile = new ReadFile(infile);
		
		if(outfile.equalsIgnoreCase("ERROR")){
			System.out.print(readFile);
		}
		else{
			readFile.rename(outfile);
		}
	} // End of main method

	//-----------------------------------------------------------------
	//-----------------------------------------------------------------
    
} // End of class ReadFile
