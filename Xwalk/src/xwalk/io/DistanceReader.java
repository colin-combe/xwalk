/*
 * (C) 2010 Abdullah Kahraman
 *
 * This software is part of the open-source project "Xwalk". You can use this
 * software under the terms of the
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 * (http://creativecommons.org/licenses/by-nc-sa/3.0/).
 * This means that you
 * 1.) can copy, modify, distribute the software
 * 2.) must give credit to the author
 * 3.) must not use this work for commercial purposes
 * 4.) must license derivative works under the same or a similar license.
 *
 */

package xwalk.io;

import java.io.IOException;

import structure.constants.Constants;
import structure.grid.GridCell.Value;
import structure.io.ReadFile;
import structure.matter.Atom;
import xwalk.crosslink.CrossLink;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.CrossLinkParameter.Parameter;

/**
 * This class converts distance files into CrossLink objects.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class DistanceReader {
    /**
     * Constructor.
     */
    protected DistanceReader() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Reads in fileName and converts these into CrossLink objects.
     * @param fileName
     *        - String object holding the path to a distance file.
     * @param onlyIntraLinks
     *        - boolean object determining whether only intra links should be
     *          read out.
     * @param onlyInterLinks
     *        - boolean object determining whether only inter links should be
     *          read out.
     * @return List of CrossLink objects extracted from the distance file.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     */
    public static CrossLinkList getCrossLinks(final String fileName,
                                              final boolean onlyIntraLinks,
                                              final boolean onlyInterLinks)
                                                            throws IOException {
        CrossLinkList set = new CrossLinkList();

        ReadFile read = new ReadFile(fileName);
        for (String line : read) {
            if (!line.startsWith("#") && line.trim().length() >= 1) {
                Atom atom1 = new Atom();
                Atom atom2 = new Atom();
                String file = "";
                int index = 0;
                int seqDist = -1;
                double eucDist = -1;
                double solvDist = 
                                Double.parseDouble(Value.DISTANCE.getDefault());
                try {
                    String[] array = line.trim().split("\t");
                    index = Integer.parseInt(array[0]);
                    file = array[1];
                    String atom1info = array[2];
                    String atom2info = array[3];
// no distance information should be read in from a distance file
// which will leave the distance information at -1, which can be used
// to recognize whether the cross-link has been found in the current structure.                    
/*                    if (array.length > 4) {
                        seqDist = Integer.parseInt(array[4]);
                    }
                    if (array.length > 5) {
                        eucDist = Double.parseDouble(array[5]);
                    } 
                    if (array.length > 6) {
                        solvDist = Double.parseDouble(
                                                     Value.DISTANCE.getDefault()
                                                     );
                        if (!array[6].equals("-")) {
                            solvDist = Double.parseDouble(array[6]);
                        }
                    }
                    if (array.length > 7) {
                        String peptideSequence = array[7];
                    }
*/
                    array = atom1info.split("-");
                    if (array.length < 2) {
                        System.err.println("WARNING: First atom of cross-link "
                                         + "number " + index + " must list a "
                                         + "residue name and residue "
                                         + "number.");
                    } else {
                        atom1.setResidueName(array[0].trim());
                        atom1.setResidueNumber(
                                               Integer.parseInt(array[1].trim())
                                              );
                    }
                    if (array.length >= 3) {
                        if(array[2].trim().length() != 0) {
                            atom1.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                        }
                    }
                    if (array.length >= 4) {
                            atom1.setName(array[3].trim());
                    }

                    array = atom2info.split("-");
                    if (array.length < 2) {
                        System.err.println("WARNING: Second atom of cross-link "
                                         + "number " + index + " must list a "
                                         + "residue name and residue "
                                         + "number.");
                    } else {
                        atom2.setResidueName(array[0].trim().toUpperCase());
                        atom2.setResidueNumber(
                                               Integer.parseInt(array[1].trim())
                                              );
                    }
                    if (array.length >= 3) {
                        if(array[2].trim().length() != 0) {
                            atom2.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                            if (onlyInterLinks) {
                               if (atom1.getChainId() == atom2.getChainId()) {
                                   continue;
                               }
                            }
                            if (onlyIntraLinks) {
                                if (atom1.getChainId() != atom2.getChainId()) {
                                    continue;
                                }
                             }
                        }
                    }
                    if (array.length >= 4) {
                        atom2.setName(array[3].trim());
                    }
                } catch (Exception e) {
                    System.err.println("WARNING: Distance file \"" + fileName
                                     + "\" does not conform to distance file "
                                     + "format" + Constants.LINE_SEPERATOR + e);
                }

                CrossLink crossLink = new CrossLink(atom1, atom2, seqDist,
                                                    eucDist);
                crossLink.setFileName(file);
                crossLink.setSolventPathDistance(solvDist);
                crossLink.setIndex(index);

                if (set.get(atom1, atom2) == null
                    &&
                    set.get(atom2, atom1) == null) {
                        set.add(crossLink);
                }
            }
        }
    return set;
    }
}
