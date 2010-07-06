package cX;

import java.io.FileWriter;
import java.io.IOException;
import java.io.File;


/**
 * Ueberschrift:  Bioinformatik
 * Beschreibung:  Bioinformatic programs. Class to write a file.
 * Copyright:     Copyright (c) 2003
 * Organisation:  FH-Giessen
 * @author Abdullah Kahraman
 * @version 1.0
 */
public class WriteFile {

	/*--------------------------------------------------------------------------*/
	// Classobjects
	/*--------------------------------------------------------------------------*/
	
	/**stream object to write a file.
	 */
	private FileWriter fileWriter;	
	
	/**name of file to write.
	 */	
	private String fileName;  
	  						
	/*--------------------------------------------------------------------------*/
	// Classmethods
	/*--------------------------------------------------------------------------*/    

	/** Constructor I
	 * @param String fileName: (path and) name of file to write.
	 */
	public WriteFile(String fileName){
		
		this.fileName = fileName;
		
		try{
			fileWriter = new FileWriter(fileName);
		}
		catch(IOException e){
			System.err.print("\n\n###Error while trying to create file \"" + fileName + "\".###" );
			System.err.print("\n###" + e + "###\n\n");
		}
	} // End of constructor I
	
	//----------------------------------------------------------
	/** Constructor II
	 * @param String fileName: (path and) name of file to write.
	 * @param boolean append: determin if fileName should be opend appendly. 
	 */
	public WriteFile(String fileName, boolean append){
		
		this.fileName = fileName;
		
		try{
			fileWriter = new FileWriter(fileName, append);
		}
		catch(IOException e){
			System.err.print("\n\n###ERROR while creating the file \"" + fileName + "\".###" );
			System.err.print("\n###" + e + "###\n\n");
		}
	} // End of constructor II

	//-----------------------------------------------------------------------------------
	/** writes the data in file.
	 * @param String data: information which has to be write in file.
	 */
	public void write(String data){
    	
    	try{
			fileWriter.write(data); 
			fileWriter.close();
			fileWriter = null;
    	}
    	catch(IOException e){
			System.err.print("\n\n###ERROR while creating the file \"" + fileName + "\".###" );
			System.err.print("\n###" + e + "###\n\n");
    	}
	} // End of method write()

	//-----------------------------------------------------------------------------------
	/** creates a directory.
	 * @param String dir: directory to create.
	 */
	public static void createDir(String dir){
    	
		try{
			new File(dir).mkdir();
				
		}
		catch(SecurityException e){
			System.err.print("\n\n###ERROR while creating directory \"" + dir + "\".###" );
			System.err.print("\n###" + e + "###\n\n");
		}
	} // End of method write()


//	 Deletes all files and subdirectories under dir.
    // Returns true if all deletions were successful.
    // If a deletion fails, the method stops attempting to delete and returns false.
    public static boolean deleteDir(String dirPath) {
    	File dir = new File(dirPath);
    	
    	if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
            	boolean success = new File(dirPath+File.separator+children[i]).delete();
                if (!success) {
                    return false;
                }
            }
        }
    
        // The directory is now empty so delete it
        return dir.delete();
    }

    // Delete a file.
    // Returns true if deletion was successful.
    // If deletion fails, the method stops attempting to delete and returns false.
    public static boolean deleteFile(String filePath) {
    	File file = new File(filePath);
    	
    	if (file.isFile()) {
    		return file.delete();
        }
    	else{
    		return false;
    	}
    	
    }

//	-----------------------------------------------------------------
//	-----------------------------------------------------------------
} // End of class WriteFile
