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

        NumberFormat decFormat = new DecimalFormat("0.0");

        // get necessary values from CrossLinkParameter object.
        String infile = parameter.getParameter(Parameter.INFILE_PATH);
        String nl = Constants.LINE_SEPERATOR;

        String infileWithoutExtension = new File(infile).getName().replaceAll(
                                                                     "\\..*", ""
                                                                             );

        output.append("load " + infile + nl);
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

        boolean emptyChainId = false;
        // maximum index number of cross-links, which will be used to set the
        // final iteration step in the for loops applied to color the
        // cross-links in PyMOL.
        int maxIndex = -1;

        for (CrossLink crossLink : crossLinkList) {

            maxIndex = Math.max(maxIndex, crossLink.getIndex());

            Atom atom1 = crossLink.getPreAtom();
            Atom atom2 = crossLink.getPostAtom();

            if (atom1.getChainId() == '_' || atom2.getChainId() == '_') {
                emptyChainId = true;
            }

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
            String selection1 = "\"resn " + atom1.getResidueName().trim()
                              + " and resi " + atom1.getResidueNumber()
                              + " and chain " + atom1.getChainId()
                              + " and name " + atom1.getName().trim() + " and "
                              + infileWithoutExtension + "\"";
            String selection2 = "\"resn " + atom2.getResidueName().trim()
                              + " and resi " + atom2.getResidueNumber()
                              + " and chain " + atom2.getChainId()
                              + " and name " + atom1.getName().trim() + " and "
                              + infileWithoutExtension + "\"";

            String distName = crossLink.getIndex() + "_";
            if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                           )
                                    )
               ) {
                distName += decFormat.format(crossLink.getSolventPathDistance())
                            + "_";
            } else {
                distName += decFormat.format(crossLink.getEuclideanDistance())
                            + "_";
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

            output.append("cmd.select(\"pk1\"," + selection1 + ")" + nl);
            output.append("cmd.select(\"pk2\"," + selection2 + ")" + nl);
            output.append("cmd.show(\"spheres\",\"pk1\")" + nl);
            output.append("cmd.show(\"spheres\",\"pk2\")" + nl);
            output.append("cmd.distance(\"" + distName
                        + "\", \"(pk1)\", \"(pk2)\")" + nl
                          );
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
            String pathFileName = infileWithoutExtension
                                + "_solventPathDistances.pdb";

            WriteFile.deleteFile(pathFileName);

            for (CrossLink crossLink : crossLinkList) {
                Atom atom1 = crossLink.getPreAtom();
                Atom atom2 = crossLink.getPostAtom();
                String distName = crossLink.getIndex() + "_"
                                + decFormat.format(
                                              crossLink.getSolventPathDistance()
                                                  ) + "_"
                                + atom1.getResidueName().trim() + ""
                                + atom1.getResidueNumber() + ""
                                + atom1.getChainId()
//                                + "_" + atom1.getName().trim()
                                + "-"
                                + atom2.getResidueName().trim() + ""
                                + atom2.getResidueNumber() + ""
                                + atom2.getChainId()
//                                + "_" + atom2.getName().trim()
                                + "_solventPath";

                WriteFile file = new WriteFile();
                file.setFile(pathFileName, true);
                file.write("HEADER " + distName + nl
                           + crossLink.getPath().toString(crossLink.getIndex())
                           + "END" + nl);
            }
            output.append("load " + pathFileName + ", solventPaths" + nl);
            output.append("hide everything, *solvent*" + nl);
            output.append("show spheres, *solvent*" + nl);
            output.append("set sphere_scale, 0.7, *solvent*" + nl);
//            output.append("color red, *solvent*" + nl);

        }

        output.append("show ribbon, chain*" + nl);
        output.append("show surface, chain*" + nl);
        if (emptyChainId || Boolean.parseBoolean(parameter.getParameter(
                                                          Parameter.IS_HOMOMERIC
                                                                       )
                                                )
           ) {
            output.append("set transparency, 0.5, " + infileWithoutExtension
                        + nl);
        } else {
            output.append("set transparency, 0.5, chain*" + nl);
        }

        output.append("for i in range(1,100): "
                    + "cmd.set_color(\"col\"+str(i), "
                    + "[1-float((i*20)%100/100), float((i*30)%100)/100,"
                    + "0])" + nl);
        output.append("for i in range(1," + (maxIndex + 1) + "): "
                    + "cmd.set(\"sphere_color\",\"col\"+str(i),str(i)+\"_*-*\")"
                    + nl
                    + "for i in range(1," + (maxIndex + 1) + "): "
                    + "cmd.set(\"dash_color\",\"col\"+str(i),str(i)+\"_*-*\")"
                    + nl
                    + "for i in range(1," + (maxIndex + 1) + "): "
                    + "cmd.set(\"label_color\",\"col\"+str(i),str(i)+\"_*-*\")"
                    + nl);
        output.append("reset" + nl);

    return super.write(output.toString());
    }
    //--------------------------------------------------------------------------
}
