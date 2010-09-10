package xwalk.io;

import java.io.IOException;

import structure.constants.Constants;
import structure.grid.GridCell.Value;
import structure.io.ReadFile;
import structure.matter.Atom;
import xwalk.crosslink.CrossLink;
import xwalk.crosslink.CrossLinkList;

/**
 * This class converts distance files into CrossLink objects.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
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
     * @return List of CrossLink objects extracted from the distance file.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     */
    public static CrossLinkList getCrossLinks(final String fileName)
                                                            throws IOException {
        CrossLinkList set = new CrossLinkList();

        ReadFile read = new ReadFile(fileName);

        for (String line : read) {
            if (!line.startsWith("#") && line.trim().length() >= 1) {
                Atom atom1 = new Atom();
                Atom atom2 = new Atom();
                int index = 0;
                int seqDist = Integer.parseInt(Value.DISTANCE.getDefault());
                double eucDist = Double.parseDouble(
                                                    Value.DISTANCE.getDefault()
                                                   );
                double solvDist = Double.parseDouble(
                                                     Value.DISTANCE.getDefault()
                                                    );
                try {
                    String[] array = line.trim().split("\t");
                    index = Integer.parseInt(array[0]);
                    String file = array[1];
                    String atom1info = array[2];
                    String atom2info = array[3];
                    seqDist = Integer.parseInt(array[4]);
                    eucDist = Double.parseDouble(array[5]);
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

                    array = atom1info.split("-");
                    atom1.setResidueName(array[0].trim());
                    atom1.setResidueNumber(Integer.parseInt(array[1].trim()));
                    atom1.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                    atom1.setName(array[3].trim());

                    array = atom2info.split("-");
                    atom2.setResidueName(array[0].trim());
                    atom2.setResidueNumber(Integer.parseInt(array[1].trim()));
                    atom2.setChainId(array[2].trim().charAt(0) == '_'
                                                    ? ' ' : array[2].charAt(0));
                    atom2.setName(array[3].trim());
                } catch (Exception e) {
                    System.err.println("WARNING: Distance file \"" + fileName
                                     + "\" does not conform to distance file "
                                     + "format" + Constants.LINE_SEPERATOR + e);
                }
                CrossLink crossLink = new CrossLink(atom1, atom2, seqDist,
                                                    eucDist);
                crossLink.setSolventPathDistance(solvDist);
                crossLink.setIndex(index);
                set.add(crossLink);
            }
        }
    return set;
    }
}
