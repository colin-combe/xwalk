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

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import structure.constants.Constants;
import structure.io.WriteFile;
import structure.matter.Atom;
import xwalk.crosslink.CrossLink;
import xwalk.crosslink.CrossLinkParameter;
import xwalk.crosslink.CrossLinkList;
import xwalk.crosslink.CrossLinkParameter.Parameter;

/**
 * This class holds the distance file format and can be used to output or write
 * CrossLink object in distance file format.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class DistanceWriter extends WriteFile {

    /**
     * Distances should be reported up to one digit after the comma.
     */
    private NumberFormat decFormat = new DecimalFormat("0.0");

    //--------------------------------------------------------------------------
    /**
     * Returns the header for the distance file.
     * @return String object holding the distance file header.
     */
    private static String getDistanceFileHeader() {
        StringBuffer output = new StringBuffer();
        output.append("#-----\t--------\t---------\t---------\t---\t---");
        output.append("\t---\t----------------");
        output.append(Constants.LINE_SEPERATOR);
        output.append("#Index\tFileName\tResi1info\tResi2info\tSeq\tEuc");
        output.append("\tSpd\tPepSeq");
        output.append(Constants.LINE_SEPERATOR);
        output.append("#-----\t--------\t---------\t---------\t---\t---");
        output.append("\t---\t----------------");
        output.append(Constants.LINE_SEPERATOR);
    return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Writes all cross-link objects in a distance file format into a file.
     * @param crossLinkList
     *        - CrossLinkList object holding all cross-links found for a protein
     *          complex.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return {@code TRUE} if creating/writing file is successful,
     *         {@code FALSE} if IOException is thrown and caught.
     */
    public final boolean writeFile(final CrossLinkList crossLinkList,
                                   final CrossLinkParameter parameter) {
        return super.write(DistanceWriter.toString(crossLinkList, parameter));
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the set of CrossLink objects in distance file format.
     * @param crossLinkList
     *        - CrossLinkList object holding all cross-links found for a protein
     *          complex.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return String object holding all cross-links in distance file format.
     */
    public static String toString(final CrossLinkList crossLinkList,
                                  final CrossLinkParameter parameter) {
        StringBuffer output = new StringBuffer();

        // get necessary values from CrossLinkParameter object.
        if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                       )
                                )
           ) {
            output.append(DistanceWriter.getDistanceFileHeader());
        }

        output.append(crossLinkList.toString());
        return output.toString();
    }
    //--------------------------------------------------------------------------
    /**
     * Writes a the shortest paths of a list of virtual cross-links into a
     * PDB formated file.
     * @param pathFileName
     *        - String object holding the path to the PDB file to be written.
     * @param crossLinkList
     *        - List of cross-link objects holding each the the shortest path
     */
    public final void writeSASDpdbFile(final String pathFileName,
                                       final CrossLinkList crossLinkList) {
        StringBuffer content = new StringBuffer();
        for (CrossLink crossLink : crossLinkList) {
            Atom atom1 = crossLink.getPreAtom();
            Atom atom2 = crossLink.getPostAtom();
            String distName = crossLink.getIndex() + "_"
                            + this.decFormat.format(
                                          crossLink.getSolventPathDistance()
                                              ) + "_"
                            + atom1.getResidueName().trim() + ""
                            + atom1.getResidueNumber() + ""
                            + atom1.getChainId()
//                            + "_" + atom1.getName().trim()
                            + "-"
                            + atom2.getResidueName().trim() + ""
                            + atom2.getResidueNumber() + ""
                            + atom2.getChainId()
//                            + "_" + atom2.getName().trim()
                            + "_solventPath";

            content.append("HEADER " + distName + Constants.LINE_SEPERATOR
                         + crossLink.getPath().toString(crossLink.getIndex())
                         + "END" + Constants.LINE_SEPERATOR);
        }
        WriteFile file = new WriteFile();
        file.setFile(pathFileName);
        file.write(content.toString());

    }

    //--------------------------------------------------------------------------
    /**
     * Writes a PyMol script, which loads the complex, highlights all virtual
     * cross-linked atoms and draws dashed lines in-between them as indications
     * for their cross-link.
     * @param crossLinkList
     *        - CrossLinkList object holding all cross-links found for a protein
     *          complex.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return {@code TRUE} if creating/writing file is successful,
     *         {@code FALSE} if IOException is thrown and caught.
     */
    public final boolean writePymolScript(
                                          final CrossLinkList crossLinkList,
                                          final CrossLinkParameter parameter
                                         ) {
        StringBuffer output = new StringBuffer();

        // get necessary values from CrossLinkParameter object.
        String infile = parameter.getParameter(Parameter.INFILE_PATH);
        String nl = Constants.LINE_SEPERATOR;

        String infileWithoutExtension = new File(infile).getName().replaceAll(
                                                                     "\\..*", ""
                                                                             );

        output.append("load " + infile + nl);
        output.append("set transparency, 0.5, " + infileWithoutExtension + nl);
        //        output.append("disable " + infileWithoutExtension + nl);
        output.append("hide everything, " + infileWithoutExtension + nl);
        output.append("set dash_radius, 1, " + infileWithoutExtension + nl);
        output.append("set dash_width, 15, " + infileWithoutExtension + nl);
        output.append("set dash_gap, 0, " + infileWithoutExtension + nl);
        output.append("bg_color white" + nl);
        output.append("util.cbc" + nl);
        output.append("create het, hetatm and " + infileWithoutExtension + nl);
        output.append("show sticks, het" + nl);
        output.append("color grey, het" + nl);
        output.append("disable het" + nl);

        // maximum index number of cross-links, which will be used to set the
        // final iteration step in the for loops applied to color the
        // cross-links in PyMOL.
        int maxIndex = -1;

        for (CrossLink crossLink : crossLinkList) {

            maxIndex = Math.max(maxIndex, crossLink.getIndex());

            Atom atom1 = crossLink.getPreAtom();
            Atom atom2 = crossLink.getPostAtom();


            if (atom1.getChainId() != '_') {
                output.append("create chain" + atom1.getChainId() + ", chain "
                             + atom1.getChainId() + " and "
                             + infileWithoutExtension + nl);
            }
            if (atom2.getChainId() != '_') {
                output.append("create chain" + atom2.getChainId()
                            + ", chain " + atom2.getChainId() + " and "
                            + infileWithoutExtension + nl);
            }
            String selection1 = "resn " + atom1.getResidueName().trim()
                              + " and resi " + atom1.getResidueNumber()
                              + " and chain " + atom1.getChainId()
                              + " and name " + atom1.getName().trim() + " and "
                              + infileWithoutExtension;
            String selection2 = "resn " + atom2.getResidueName().trim()
                              + " and resi " + atom2.getResidueNumber()
                              + " and chain " + atom2.getChainId()
                              + " and name " + atom1.getName().trim() + " and "
                              + infileWithoutExtension;

            String distName = crossLink.getIndex() + "_";
            if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                           )
                                    )
               ) {
                distName += this.decFormat.format(
                                              crossLink.getSolventPathDistance()
                                                 ) + "_";
            } else {
                distName += this.decFormat.format(
                                                crossLink.getEuclideanDistance()
                                                 ) + "_";
            }

            distName += atom1.getResidueName().trim() + ""
                      + atom1.getResidueNumber() + ""
                      + atom1.getChainId()
//                      + "_" + atom1.getName().trim()
                      + "-"
                      + atom2.getResidueName().trim() + ""
                      + atom2.getResidueNumber() + ""
                      + atom2.getChainId();
//                      + "_" + atom2.getName().trim();

            output.append("select pk1, " + selection1 + nl);
            output.append("select pk2, " + selection2 + nl);
            output.append("show spheres, pk1" + nl);
            output.append("show spheres, pk2" + nl);
            output.append("distance " + distName + nl);
//            output.append("cmd.color(\"red\", \"" + distName + "\")" + nl);
        }

        output.append("delete pk1" + nl);
        output.append("delete pk2" + nl);

        // Write solvent path distances into a file to be loaded
        if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                      )
                               )
          ) {

            String outputDir = new File(parameter.getParameter(
                                                Parameter.OUTFILE_PATH
                                                              )).getParent();
            if (outputDir == null) {
                outputDir = "";
            }

            String sasdFileName = infileWithoutExtension
                                + "_solventPathDistances.pdb";
            String pathFileName = outputDir
                                + File.separatorChar
                                + sasdFileName;

            WriteFile.deleteFile(pathFileName);

            this.writeSASDpdbFile(pathFileName, crossLinkList);

            output.append("load " + sasdFileName + ", solventPaths" + nl);
            output.append("hide everything, *solvent*" + nl);
            output.append("show spheres, *solvent*" + nl);
//            output.append("set sphere_scale, 0.7, *solvent*" + nl);
//            output.append("color red, *solvent*" + nl);

        }

        output.append("show ribbon, chain*" + nl);
        output.append("show surface, chain*" + nl);

        int n = crossLinkList.size();
        output.append("for i in range(1," + n + 1 + "): "
                    + "cmd.set_color(\"col\"+str(i), "
                    + "[1-float((i*20)%" + n + "/" + n + "), float((i*30)%" + n
                    + ")/" + n + ",0])" + nl);
        for (int i = 0; i < n; i++) {
            output.append("set sphere_color, col" + (i + 1) + ", "
                        + crossLinkList.get(i).getIndex() + "_*-*"
                        + nl
                        + "set dash_color, col" + (i + 1) + ", "
                        + crossLinkList.get(i).getIndex() + "_*-*"
                        + nl
                        + "set label_color, col" + (i + 1) + ", "
                        + crossLinkList.get(i).getIndex() + "_*-*"
                        + nl);
        }
        output.append("reset" + nl);

    return super.write(output.toString());
    }
    //--------------------------------------------------------------------------
}
