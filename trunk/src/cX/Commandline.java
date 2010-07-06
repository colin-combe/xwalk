package cX;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Ueberschrift:  PropSearch
 * Beschreibung:  Commandline. Klasse zum Erhalten von Parameter aus der Commandshell. 
 * Copyright:     Copyright (c) 2003
 * Organisation:  FH-Giessen
 * @author Abdullah Kahraman
 * @version 1.0
 */


public class Commandline {

  /*--------------------------------------------------------------------------*/
  // Classmethods
  /*--------------------------------------------------------------------------*/

  public Commandline() {
  }
  /*--------------------------------------------------------------------------*/

  /**@param args - commandline parameter. Normally args from the main-method
   *        param - desired parameter for look for
   * @return - if the desired parameter is not found, "ERROR" is returned
   *
   * read a desired parameter from the commandline
   */
  	public String get(String[] args, String param, boolean hasParameter){
	// loop over all parameters in the commandline
    	for(int i=0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase(param) && !hasParameter){
				return "EXISTS";
			}
			if (args[i].equalsIgnoreCase(param) && hasParameter && args.length > i+1){
				if(args[i+1].charAt(0) != '-'){
					return args[i+1];
				}
				
				else{
//					System.err.print("\nCouldn't find value for parameter \""+param+"\"\n\n");
//					System.exit(2);				
					// if required the next parameter or negative number can be extracted from here, by deleting EXISTS from return string
					return "EXISTS"+args[i+1];
				}
			}
    	}
	return "ERROR";

	} // End of method get()

    /*--------------------------------------------------------------------------*/
    /**@param null
     * @return - string which was typed on the commandline by the user.
     */
  	static public String get(){
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
	    try {
	    	return br.readLine();
	    } catch (IOException ioe) {
	    	return "";
	    }
  	}
  /*--------------------------------------------------------------------------*/
} // End of class Commandline