package xwalk.crosslink;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.TreeSet;
import java.util.zip.DataFormatException;

import structure.constants.Constants;
import structure.constants.Constants.ParameterSets;
import structure.grid.AtomGrid;
import structure.grid.GridUtilities;
import structure.grid.Path;
import structure.grid.GridCell.Value;
import structure.io.pdb.PDBreader;
import structure.math.Mathematics;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.parameter.AtomType;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.Digestion;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;

import xwalk.io.DistanceReader;
import xwalk.math.SolventPathDistance;
import xwalk.crosslink.CrossLinkParameter.Parameter;

/**
 * Class that holds all relevant methods to calculate virtual cross-links on PDB
 * protein structures.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public final class CrossLinkUtilities {

    /**
     * Constructor.
     */
    protected CrossLinkUtilities() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------
    /**
     * Returns a list of potential virtual cross-links.
     * @param complexes
     *        - List of protein complex objects holding all protein complex
     *          molecules for which virtual cross-links should be calculated.
     * @param parameter -
     *        CrossLinkParameter object, holding all parameter that are
     *        necessary for the virtual cross-link calculation.
     * @throws IOException if an error occurred while reading the infile.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     * @return CrossLinkList object that holds all virtual cross-links on a
     *         protein complex.
     */
    public static CrossLinkList getVirtualCrossLinks(
                                  final ArrayList < PolyPeptideList > complexes,
                                  final CrossLinkParameter parameter
                                                    )
                                                  throws IOException,
                                                         DataFormatException {
        //--------------------------------
        // read in cross-links from distance file
        CrossLinkList distXl = null;
        if (!parameter.getParameter(Parameter.DISTANCE_FILE_PATH).equals("")) {
            distXl = DistanceReader.getCrossLinks(parameter.getParameter(
                                                    Parameter.DISTANCE_FILE_PATH
                                                                        )
                                                 );
        }
        //--------------------------------
        // digest protein
        ArrayList < PolyPeptideList > allDigest = null;
        if (Boolean.parseBoolean(parameter.getParameter(
                                                 Parameter.DO_TRYPSIN_DIGEST
                                                       )
                                )
        ) {
            allDigest = CrossLinkUtilities.digestProteinComplex(
                                                             parameter,
                                                             complexes
                                                            );
        }

        //--------------------------------
        // calculate distances

        CrossLinkList allCrossLinkList = new CrossLinkList();
        // find and create virtual cross-links on the protein complexes.
        for (int i = 0; i < complexes.size(); i++) {
            PolyPeptideList complex = complexes.get(i);
            PolyPeptideList digest = null;

            if (allDigest != null) {
                digest = allDigest.get(i);
            }

            //---------------------------------
            // First find cross-links based on Euclidean distance.
            CrossLinkList crossLinkList =
                            CrossLinkUtilities.crossLinkByEuclideanDistance(
                                                                      parameter,
                                                                      complex,
                                                                      distXl,
                                                                      digest
                                                                           );
            //---------------------------------
            // remove redundant cross-links if the complex should be labeled by
            // the user as homomeric.
            if (Boolean.parseBoolean(parameter.getParameter(
                                                          Parameter.IS_HOMOMERIC
                                                           ))) {
             CrossLinkUtilities.removeRedundanciesInHomomers(crossLinkList);
            }

            //---------------------------------
            // If requested by the user check further for Solvent Path distance.
            if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                           ))) {
                if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                               ))) {
                    System.err.print("Checking \"" + crossLinkList.size() + "\""
                                   + " Euclidean distances for becoming solvent"
                                   + " accessible surface distances.\n");
                }
                crossLinkList =
                           CrossLinkUtilities.calculatesSolventPathDistance(
                                                                   parameter,
                                                                   complex,
                                                                   crossLinkList
                                                                           );
            }
            //---------------------------------
            // set file name of each grid to the file name of its protein
            // complex and assign peptide sequences to cross-linked atoms.
            for (CrossLink xlink : crossLinkList) {
                xlink.setFileName(complex.getName());
            }

            //---------------------------------
            // sort list of cross-links by distance.
            crossLinkList.sort();

            //---------------------------------
            // set indices of cross-links
            CrossLinkUtilities.setCrossLinkIndices(crossLinkList, parameter);

            allCrossLinkList.addAll(crossLinkList);
        }

        return allCrossLinkList;
    }
    //--------------------------------------------------------------------------
    /**
     * Digest a list of protein complex objects.
     * @param complexes
     *        List of Protein complex object to be all digested.
     * @param parameter
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return List of PolyPeptides that are formed by digestion.
     */
    public static ArrayList <PolyPeptideList> digestProteinComplex(
            final CrossLinkParameter parameter,
            final ArrayList < PolyPeptideList > complexes
        ) {
        ArrayList < PolyPeptideList > allDigest =
                                            new ArrayList < PolyPeptideList >();
        for (PolyPeptideList complex : complexes) {

            PolyPeptideList digest = CrossLinkUtilities.trypsinate(complex,
                                                   Boolean.parseBoolean(
                                                      parameter.getParameter(
                                                        Parameter.DO_EXPASY_RULE
                                                                            )
                                                                       )
                                                  );
            allDigest.add(digest);
            if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                           ))) {
                // output digested peptides
                for (int i = 0; i < digest.size(); i++) {
                    System.err.println(i + 1 + ". "
                                     + digest.get(i).toStringOneLetterCode());
                }
            }
        }
        return allDigest;
    }

    //--------------------------------------------------------------------------
    /**
     * CrossLinks all atoms in a protein complex that have an Euclidean distance
     * smaller than a user set maxDist value and if set by the user are found
     * in a distance file.
     * @param complex
     *        Protein complex object.
     * @param parameter
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @param distanceFileCrossLinks
     *        List of CrossLink objects extracted from a distance file. If no
     *        distance file has been read, just submit a {@code NULL}.
     * @param digest
     *        - List of PolyPeptides that are formed by digestion. If no
     *          digested peptides exist, than submit {@code NULL}.
     * @return List of CrossLink object that all have a Euclidean distance
     *         smaller then the user set maxDist.
     *         and have a Euclidean distance smaller then the user set maxDist.
     * @throws IOException if an error occurred while reading the distance file.
     */
    public static CrossLinkList crossLinkByEuclideanDistance(
                                     final CrossLinkParameter parameter,
                                     final PolyPeptideList complex,
                                     final CrossLinkList distanceFileCrossLinks,
                                     final PolyPeptideList digest)
                                            throws IOException {
        Hashtable < Atom, AtomList > relevantAtomPairs =
            new Hashtable < Atom, AtomList >();

        if (distanceFileCrossLinks == null) {
            relevantAtomPairs = CrossLinkUtilities.findRelevantPairs(
                                                                  parameter,
                                                                  complex
                                                                    );
        } else {
            relevantAtomPairs = CrossLinkUtilities.extractRelevantPairs(
                                                     complex,
                                                     distanceFileCrossLinks,
                                                     parameter
                                                                        );
        }


        // create CrossLinks object from all relevant atom pairs.
        CrossLinkList crossLinks = new CrossLinkList();
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            PolyPeptide atom1TrypticPeptide = null;
            for (Atom atom2 : relevantAtomPairs.get(atom1)) {
                PolyPeptide atom2TrypticPeptide = null;
                if (digest == null) {
                    CrossLink xl = new CrossLink(atom1, atom2);
                    crossLinks.add(xl);
                } else {
                    for (PolyPeptide peptide : digest) {
                        for (AminoAcid aa : peptide) {
                            if (aa.getAllAtoms().contains(atom1)) {
                                atom1TrypticPeptide = peptide;
                                // Added these lines to allow for
                                // self cross-links within a peptide.
                                if (aa.getAllAtoms().contains(atom2)) {
                                    atom2TrypticPeptide = peptide;
                                }
                                break;
                            }
                            if (aa.getAllAtoms().contains(atom2)) {
                                atom2TrypticPeptide = peptide;
                                break;
                            }
                        }
                    }
                    if (atom1TrypticPeptide != null
                        &&
                        atom2TrypticPeptide != null) {
                        CrossLink xl = new CrossLink(atom1, atom2);
                        xl.setPeptides(atom1TrypticPeptide,
                                       atom2TrypticPeptide);
                        crossLinks.add(xl);
                    }
                }
            }
        }

        return crossLinks;
    }

    //--------------------------------------------------------------------------
    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance using a single grid
     * object.
     * @param complex
     *      - PolyPeptideList object holding all atoms of the protein.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkList object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkList calculateCrossLinksOnGlobalGrid(
                              final PolyPeptideList complex,
                              final CrossLinkList crossLinksByEuclideanDistance,
                              final CrossLinkParameter parameter
                                                              ) {

        boolean doSAS = Boolean.parseBoolean(parameter.getParameter(
                                                                Parameter.DO_SAS
                                                                   )
                                            );
        // set solvent radius
        double solventRadius = 0;
        if (doSAS) {
            solventRadius = Double.parseDouble(parameter.getParameter(
                                                        Parameter.SOLVENT_RADIUS
                                                                     )
                                              );
        }
        // create a single global atom grid for entire protein complex.
        AtomGrid grid = new AtomGrid(
                                     complex.getAllAtoms(),
                                     Double.parseDouble(parameter.getParameter(
                                                        Parameter.GRID_CELL_SIZE
                                                                              )
                                                       ),
                                     solventRadius + 1,
                                     doSAS
                                    );

        // measure distances in the atom grid.
        PolyPeptideList dummy = null;
        CrossLinkList xlinkSet = CrossLinkUtilities.evaluateSolventPathDistance(
                                                  dummy,
                                                  grid,
                                                  crossLinksByEuclideanDistance,
                                                  parameter
                                                                              );

        return xlinkSet;
    }
    //--------------------------------------------------------------------------
    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance using a grid object
     * for each atom1 object.
     * @param complex
     *      - PolyPeptideList object holding all atoms of the protein.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkList object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkList calculateCrossLinksOnLocalGrid(
                              final PolyPeptideList complex,
                              final CrossLinkList crossLinksByEuclideanDistance,
                              final CrossLinkParameter parameter
                                                               ) {
        AtomGrid dummy = null;
        CrossLinkList xlinkSet =
            CrossLinkUtilities.evaluateSolventPathDistance(
                                                  complex,
                                                  dummy,
                                                  crossLinksByEuclideanDistance,
                                                  parameter
                                                         );
        return xlinkSet;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts all necessary atom coordinates from the PDB file as defined in
     * the CrossLinkParameters.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return List of PolyPeptideList objects, each holding all protein
     *         coordinates labeled as ATOM up to an END flag or end of file.
     * @throws IOException if input file could not be read.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     */
    public static ArrayList < PolyPeptideList > getComplexesCoordinates(
                                             final CrossLinkParameter parameter
                                                                       )
                                                   throws IOException,
                                                           DataFormatException {

        ArrayList < PDBreader > pdbReaders =
            PDBreader.createPDBreaders(parameter.getParameter(
                                                           Parameter.INFILE_PATH
                                                             ));

        ArrayList < PolyPeptideList > proteinComplexes =
              CrossLinkUtilities.extractProteinComplexes(pdbReaders, parameter);

        // assign vdW radius to protein atoms
        for (PolyPeptideList polyPeptideList : proteinComplexes) {
             CrossLinkUtilities.setRadius(polyPeptideList, parameter);
        }

    return proteinComplexes;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts user set chain and alternative location based PDBcomplex objects
     * from PDBreader objects.
     * @param pdbReaders
     *        - List of PDBreader objects each holding the content of a single
     *          PDB file.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return List of PolyPeptideList objects holding only coordinates
     *         of atoms with user defined chain and alternative location
     *         information.
     */
    public static ArrayList < PolyPeptideList > extractProteinComplexes(
                                       final ArrayList < PDBreader > pdbReaders,
                                       final CrossLinkParameter parameter
                                                                      ) {
        ArrayList < PolyPeptideList > proteinComplexes =
            new ArrayList < PolyPeptideList >();

        if (parameter.getParameter(
                Parameter.CHAIN_ID1).equals(
                        parameter.getParameter(Parameter.CHAIN_ID2)
                                           )
                                   &&
                 parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1).equals(
                 parameter.getParameter(Parameter.ALTERNATIVE_LOCATION2)
                                   )
            ) {
            // as the restrictions are equal for both ends of the virtual
            // cross-links, it is unimportant which Id informations are taken.
            for (PDBreader reader : pdbReaders) {
                proteinComplexes.addAll(reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                                                         )
                                );
            }
        } else {
            for (PDBreader reader : pdbReaders) {
                proteinComplexes.addAll(reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1)
                       + parameter.getParameter(Parameter.CHAIN_ID2),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                       + parameter.getParameter(Parameter.ALTERNATIVE_LOCATION2)
                                                         )
                                );
            }
        }

        if (Boolean.parseBoolean(parameter.getParameter(
                                                 Parameter.DO_BACKBONE_READ))) {
            // replace full atom coordinates of protein with backbone and
            // beta-carbon only coordinates
            for (PolyPeptideList proteinComplex : proteinComplexes) {
                for (PolyPeptide protein : proteinComplex) {
                    for (int i = 0; i < protein.size(); i++) {
                        AtomList backbone = new AtomList();
                        for (Atom atom : protein.get(i).getAllAtoms()) {
                            if (atom.getType() == AtomType.CARBON_ALPHA
                               ||
                               atom.getType() == AtomType.CARBON_BETA
                               ||
                               atom.getType() == AtomType.NITROGEN
                               ||
                               atom.getType() == AtomType.OXYGEN) {
                                backbone.add(atom);
                            }
                        }
                        AminoAcid bb = new AminoAcid(backbone);
                        protein.set(i, bb);
                    }
                }
            }
        }

        return proteinComplexes;
    }

    //--------------------------------------------------------------------------
    /**
     * Digest all protein components of a protein complex.
     * @param proteinComplex
     *        - PolyPeptideList objects to be digested.
     * @param useExpasyRules
     *        - boolean value indicating to use
     *          <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">ExPASy Exception rules</a>
                for digestion.
     * @return new PolyPeptideList object with PolyPeptide objects holding the
     *         digested peptide segments.
     */
    public static PolyPeptideList trypsinate(
                           final PolyPeptideList proteinComplex,
                           final boolean useExpasyRules
                                                     ) {
        PolyPeptideList digestedComplex = new PolyPeptideList();
        for (PolyPeptide protein : proteinComplex) {
            ArrayList < PolyPeptide > digest =
                                     Digestion.partialTrypticDigest(
                                                                  protein,
                                                                  useExpasyRules
                                                                   );
            digestedComplex.addAll(digest);
            digestedComplex.setName(proteinComplex.getName());
        }
        return digestedComplex;
    }

    //--------------------------------------------------------------------------
    /**
     * Sets the van der Waals radius of the protein complex atoms appropriately,
     * either only to SURFNET radii or if SAS calculation should be carried out
     * to SURFNET radius + Solvent molecule radius.
     * @param complex -
     *        Protein complex object.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @throws IOException if an error occurs while reading the parameter file.
     */
    private static void setRadius(final PolyPeptideList complex,
                                  final CrossLinkParameter parameter)
                                                            throws IOException {

        // Finally set atom radii to SURFNET ones.
        complex.setAtomRadii(ParameterSets.SURFNET);

        // If cross-links should be excluded by SAS, then van der Waals radii
        // must be increased by the radius of the solvent molecule.
        if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                       )
                                )
           &&
           Boolean.parseBoolean(parameter.getParameter(Parameter.DO_SAS))) {

            double solventRadius = Double.parseDouble(parameter.getParameter(
                                                        Parameter.SOLVENT_RADIUS
                                                                            )
                                                     );
            for (PolyPeptide protein : complex) {
                for (AminoAcid aa : protein) {
                    for (Atom atom : aa.getAllAtoms()) {
                        atom.setVanDerWaalsRadius(atom.getVanDerWaalsRadius()
                                                + solventRadius);
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Returns pairs of atoms that conform to the atom and amino acid
     * identifiers set by the user and which have a Euclidean distance smaller
     * then the user set maxDist.
     * @param complex -
     *        Protein complex object.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return Hashtable of atom pairs that conform to the user set identifiers.
     */
    private static Hashtable < Atom, AtomList > findRelevantPairs(
                                             final CrossLinkParameter parameter,
                                             final PolyPeptideList complex
                                                                 ) {
        ArrayList <TreeSet < AtomList >> relevantAtoms =
                                       new ArrayList < TreeSet <AtomList >>();
        relevantAtoms.add(CrossLinkUtilities.findAllRelevantAtoms1(
                                                                   parameter,
                                                                   complex)
                                                                  );
        relevantAtoms.add(CrossLinkUtilities.findAllRelevantAtoms2(
                                                                   parameter,
                                                                   complex)
                                                                  );

        Hashtable < Atom, AtomList > pairs =
                    CrossLinkUtilities.createPairsBetweenRelevantAtoms(
                                                           relevantAtoms.get(0),
                                                           relevantAtoms.get(1),
                                                           parameter
                                                                      );

        pairs = CrossLinkUtilities.fixIntraInterSelection(pairs, parameter);

        return pairs;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts all relevant pairs of atoms from a user given distance file.
     * @param complex -
     *        Protein complex object.
     * @param crossLinks
     *        - List of CrossLink objects extracted from the distance file.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return Hashtable of atom pairs that conform to the user set identifiers.
     * @throws IOException if input file could not be read.
     */
    private static Hashtable < Atom, AtomList > extractRelevantPairs(
                                             final PolyPeptideList complex,
                                             final CrossLinkList crossLinks,
                                             final CrossLinkParameter parameter
                                                                    )
                                                            throws IOException {

        Hashtable < Atom, AtomList > pairs = new Hashtable < Atom, AtomList >();

        // get all cross-link atoms
        AtomList uniqueCrossLinkAtoms = new AtomList();
        for (CrossLink xl : crossLinks) {
             Atom preAtom = xl.getPreAtom();
             Atom postAtom = xl.getPostAtom();
             if (!uniqueCrossLinkAtoms.contains(preAtom)) {
                 uniqueCrossLinkAtoms.add(preAtom);
             }
             if (!uniqueCrossLinkAtoms.contains(postAtom)) {
                 uniqueCrossLinkAtoms.add(postAtom);
             }
        }

        // get all cross-linked atoms from the complex
        AtomList uniqueCrossLinkAtomsInComplex = new AtomList();
        for (Atom atom1 : complex.getAllAtoms()) {
             for (Atom atom2 : uniqueCrossLinkAtoms) {
                  if (MatterUtilities.equalsResidue(atom1, atom2)
                      &&
                      atom1.getName().trim().equals(atom2.getName().trim())) {
                          uniqueCrossLinkAtomsInComplex.add(atom1);
                          break;
                  }
             }
        }

        for (CrossLink xl : crossLinks) {
             Atom preAtom = xl.getPreAtom();
             Atom postAtom = xl.getPostAtom();

             for (Atom atom : uniqueCrossLinkAtomsInComplex) {
                 if (MatterUtilities.equalsResidue(atom, preAtom)
                     &&
                     atom.getName().trim().equals(preAtom.getName().trim())) {
                     preAtom = atom;
                 }
                 if (MatterUtilities.equalsResidue(atom, postAtom)
                     &&
                     atom.getName().trim().equals(postAtom.getName().trim())) {
                     postAtom = atom;
                 }
             }

             if (preAtom == xl.getPreAtom() || postAtom == xl.getPostAtom()) {
                 // Output a WARNING message that atoms in distance files could
                 // not be found in the input file
                 xl.setFileName(parameter.getParameter(Parameter.INFILE_PATH));
                 System.err.print("WARNING: ");
                 if (preAtom == xl.getPreAtom()
                     &&
                     postAtom != xl.getPostAtom()) {
                     System.err.print("1st atom ");
                 } else if (preAtom != xl.getPreAtom()
                            &&
                            postAtom == xl.getPostAtom()) {
                     System.err.print("2nd atom ");
                 } else {
                     System.err.print("Both atoms ");
                 }
                 System.err.print("could not be found : " + xl.getIndex() + "\t"
                                + xl.toString());
                 continue;
             }

             double dist = Double.parseDouble(
                     xwalk.constants.Constants.DISTANCE_DEC_FORMAT.format(
                                                     Mathematics.distance(
                                                           preAtom.getPoint3d(),
                                                           postAtom.getPoint3d()
                                                                         )
                                                                 ));
             double errorRange = Constants.getCoordinateUncertainty(preAtom)
                                 +
                                 Constants.getCoordinateUncertainty(postAtom);

             // check whether both cross-linked atoms in the distance file
             // have in this complex a distance smaller than maxDist.
             if (dist > Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                 )
                                          ) + errorRange
                ) {
                 continue;
             }
             AtomList pair = pairs.get(preAtom);
             if (pair == null) {
                 pair = new AtomList();
                 pair.add(postAtom);
                 pairs.put(preAtom, pair);
             } else {
                 pair.add(postAtom);
                 pairs.put(preAtom, pair);
             }
        }
        return pairs;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates cross-link objects using the Solvent-Path distance algorithm
     * either on a global or local grid, depending on the user setting and the
     * dimension of the protein complex.
     * @param parameter -
     *        CrossLinkParameter object, holding all parameter that are
     *        necessary for the virtual cross-link calculation.
     * @param complex -
     *        Protein complex object.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     * @return List of CrossLink objects.
     */
    private static CrossLinkList calculatesSolventPathDistance(
                            final CrossLinkParameter parameter,
                            final PolyPeptideList complex,
                            final CrossLinkList crossLinksByEuclideanDistance
                                                              ) {

        boolean verbose = Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                                     )
                                              );

        if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_GLOBAL_GRID
                                                       )
                                )
           ) {
            if (verbose) {
                System.err.print("Solvent Path Distances are calculated on "
                               + "global grid." + Constants.LINE_SEPERATOR);
            }
            return CrossLinkUtilities.calculateCrossLinksOnGlobalGrid(
                                                  complex,
                                                  crossLinksByEuclideanDistance,
                                                  parameter
                                                                     );
        } else {
            if (verbose) {
                System.err.print("Solvent Path Distances are calculated on "
                               + "local grids." + Constants.LINE_SEPERATOR);
            }
            return CrossLinkUtilities.calculateCrossLinksOnLocalGrid(
                                                  complex,
                                                  crossLinksByEuclideanDistance,
                                                  parameter
                                                                    );
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the indices of the cross-linked objects either iteratively or if
     * a distance file is set by the user according to the indices in the
     * distance file.
     * @param crossLinkList
     *        - List of cross-link object.
     * @param parameter -
     *        CrossLinkParameter object, holding all parameter that are
     *        necessary for the virtual cross-link calculation.
     * @throws IOException if an error occurs while reading a distance file.
     */
    private static void setCrossLinkIndices(final CrossLinkList crossLinkList,
                                            final CrossLinkParameter parameter)
                                                            throws IOException {
        if (parameter.getParameter(Parameter.DISTANCE_FILE_PATH).equals("")) {
            for (int i = 0; i < crossLinkList.size(); i++) {
                crossLinkList.get(i).setIndex(i + 1);
            }
        } else {
            CrossLinkList distanceFileCrossLinks = DistanceReader.getCrossLinks(
                            parameter.getParameter(Parameter.DISTANCE_FILE_PATH)
                                                                              );
            // assign indices of distance file to newly found cross-links.
            if (distanceFileCrossLinks.size() != 0) {
                for (CrossLink dxl : distanceFileCrossLinks) {
                    for (CrossLink xl : crossLinkList) {
                        if (dxl.equals(xl)) {
                            xl.setIndex(dxl.getIndex());
                        }
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Returns all atoms that conform to the identifier of the first atom as set
     * by the user.
     * @param complex -
     *        Protein complex object.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return Set of AtomList objects that hold all atoms of amino acids
     *         that conform to the user set identifiers.
     */
    private static TreeSet < AtomList > findAllRelevantAtoms1(
                                             final CrossLinkParameter parameter,
                                             final PolyPeptideList complex) {

        TreeSet < AtomList > candidates1 = new TreeSet < AtomList >(
                new Comparator<AtomList>() {
                    public int compare(final AtomList a1, final AtomList a2) {
                        for (Atom atom1 : a1) {
                            for (Atom atom2 : a2) {
                                if (atom1.equals(atom2)) {
                                    return 0;
                                }
                            }
                        }
                        return 1;
                    }
                });

        for (PolyPeptide protein : complex) {
            for (AminoAcid residue : protein) {
                if (CrossLinkUtilities.isAminoAcid1Relevant(residue,
                                                            parameter)
                                                           ) {
                    if (parameter.getParameter(
                                                Parameter.ATOM_TYPE1).equals("")
                                              ) {
                        candidates1.add(residue.getAllAtoms());
                    } else {
                        AtomList list = new AtomList();
                        for (Atom atom : residue.getAllAtoms()) {
                            if (CrossLinkUtilities.isAtomRelevant1(atom,
                                                                   parameter)) {
                                list.add(atom);
                            }
                        }
                        if (list.size() != 0) {
                            candidates1.add(list);
                        } else {
                            System.err.println("WARNING: "
                                 + parameter.getParameter(Parameter.INFILE_PATH)
                                 + "\tAtom \"-aa1 "
                                 + parameter.getParameter(Parameter.ATOM_TYPE1)
                                 + "\" not found in residue "
                                 + residue.getAtom(0).getResidueName()
                                 + residue.getAtom(0).getResidueNumber()
                                 + residue.getAtom(0).getChainId());
                        }
                    }
                }
            }
        }
        return candidates1;
    }
    //--------------------------------------------------------------------------

    /**
     * Returns all atoms that conform to the identifier of the second atom as
     * set by the user.
     * @param complex
     *        - Protein complex object
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return TreeSet of AtomList objects that hold all atoms of amino acids
     *         that conform to the user set identifiers.
     */
    private static TreeSet < AtomList > findAllRelevantAtoms2(
                                             final CrossLinkParameter parameter,
                                             final PolyPeptideList complex
                                                               ) {

        TreeSet < AtomList > candidates2 = new TreeSet < AtomList >(
                new Comparator<AtomList>() {
                    public int compare(final AtomList a1, final AtomList a2) {
                        for (Atom atom1 : a1) {
                            for (Atom atom2 : a2) {
                                if (atom1.equals(atom2)) {
                                    return 0;
                                }
                            }
                        }
                        return 1;
                    }
                });

        StringBuffer dataNotFoundMessage = new StringBuffer();

        for (PolyPeptide protein : complex) {
            for (AminoAcid residue : protein) {
                if (CrossLinkUtilities.isAminoAcid2Relevant(residue,
                                                            parameter)) {
                    if (parameter.getParameter(
                                               Parameter.ATOM_TYPE2).equals(""))
                    {
                        candidates2.add(residue.getAllAtoms());
                    } else {
                        AtomList list = new AtomList();
                        for (Atom atom : residue.getAllAtoms()) {
                            if (CrossLinkUtilities.isAtomRelevant2(atom,
                                                                   parameter)) {
                                list.add(atom);
                            }
                        }
                        if (list.size() != 0) {
                            candidates2.add(list);
                        } else {
                            dataNotFoundMessage.append("WARNING: "
                                 + parameter.getParameter(Parameter.INFILE_PATH)
                                 + "\tAtom \"-aa2 "
                                 + parameter.getParameter(Parameter.ATOM_TYPE2)
                                 + "\" not found in residue "
                                 + residue.getAtom(0).getResidueName()
                                 + residue.getAtom(0).getResidueNumber()
                                 + residue.getAtom(0).getChainId()
                                 + Constants.LINE_SEPERATOR);
                        }
                    }
                }
            }
        }
    return candidates2;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an AminoAcid object conforms to the name, number and chain
     * Id of the first amino acid as set by the user.
     * @param acid
     *        - Amino acid object
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return {@code TRUE} if atom corresponds to the first amino acid
     *         identifier as set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAminoAcid1Relevant(
                                              final AminoAcid acid,
                                              final CrossLinkParameter parameter
                                               ) {
        boolean residueNumberFound = parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME1
                                                           ).equals("")
                                     &&
                                     parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                           ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                    ) != -1;
        boolean residueNameAndNumberFound = !parameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NAME1
                                                                   ).equals("")
                                            &&
                                         !parameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                                ).equals("-999")
                                            &&
                                          parameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NAME1
                                                                ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                         ) != -1
                                            &&
                                           parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                                 ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                        ) != -1;

        boolean residueNameFound = !parameter.getParameter(
                                           Parameter.AMINO_ACID_RESIDUE_NAME1
                                                          ).equals("")
                                    &&
                                    parameter.getParameter(
                                           Parameter.AMINO_ACID_RESIDUE_NUMBER1
                                                          ).equals("-999")
                                    &&
                                    parameter.getParameter(
                                           Parameter.AMINO_ACID_RESIDUE_NAME1
                                                          ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                   ) != -1;

        boolean chainFound = parameter.getParameter(
                                                    Parameter.CHAIN_ID1
                                                    ).indexOf(
                           Character.toString(acid.getAtom(0).getChainId())
                                                             ) != -1;
        if (residueNumberFound || residueNameAndNumberFound
                || residueNameFound) {
            if (chainFound) {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an AminoAcid object conforms to the name, number and chain
     * Id of the second amino acid as set by the user.
     * @param acid
     *        - Amino acid object
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return {@code TRUE} if atom corresponds to the second amino acid
     *         identifier as set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAminoAcid2Relevant(
                                              final AminoAcid acid,
                                              final CrossLinkParameter parameter
                                               ) {
        boolean residueNumberFound = parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME2
                                                           ).equals("")
                                     &&
                                     parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                           ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                     ) != -1;

        boolean residueNameAndNumberFound = !parameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NAME2
                                                                   ).equals("")
                                            &&
                                         !parameter.getParameter(
                                          Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                                ).equals("-999")
                                            &&
                                          parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME2
                                                                ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                         ) != -1
                                            &&
                                          parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                                ).indexOf(
                                              "#"
                                            + acid.getAtom(0).getResidueNumber()
                                            + "#"
                                                                        ) != -1;

        boolean residueNameFound = !parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME2
                                                          ).equals("")
                                    &&
                                    parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NUMBER2
                                                          ).equals("-999")
                                    &&
                                    parameter.getParameter(
                                            Parameter.AMINO_ACID_RESIDUE_NAME2
                                                          ).indexOf(
                                                "#"
                                              + acid.getAtom(0).getResidueName()
                                              + "#"
                                                                   ) != -1;
        boolean chainFound = parameter.getParameter(
                                                    Parameter.CHAIN_ID2
                                                   ).indexOf(
                                Character.toString(acid.getAtom(0).getChainId())
                                                            ) != -1;

        if (residueNumberFound || residueNameAndNumberFound
            || residueNameFound) {
            if (chainFound) {
                return true;
            }
        }
        return false;
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an Atom object conforms to the name and alternative
     * location of the first atom as set by the user.
     * @param atom
     *        - Atom object
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return {@code TRUE} if atom corresponds to the first atom identifier as
     *         set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAtomRelevant1(final Atom atom,
                                           final CrossLinkParameter parameter) {
        if (!parameter.getParameter(Parameter.ATOM_TYPE1).equals("")) {
            if (parameter.getParameter(Parameter.ATOM_TYPE1).indexOf(
                                                           "#"
                                                         + atom.getName().trim()
                                                         + "#"
                                                                    ) != -1
                &&
                parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1).indexOf(
                             Character.toString(atom.getAlternativeLocation())
                                                                      ) != -1) {
                return true;
            }
            return false;
        } else {
            return true;
        }
    }
    //--------------------------------------------------------------------------

    /**
     * Checks whether an Atom object conforms to the name and alternative
     * location of the second atom as set by the user.
     * @param atom
     *        - Atom object
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return {@code TRUE} if atom corresponds to the second atom identifier as
     *         set by the user, {@code FALSE} otherwise.
     */
    private static boolean isAtomRelevant2(
                                           final Atom atom,
                                           final CrossLinkParameter parameter
                                          ) {
        if (!parameter.getParameter(Parameter.ATOM_TYPE2).equals("")) {
            if (parameter.getParameter(Parameter.ATOM_TYPE2).indexOf(
                                                           "#"
                                                         + atom.getName().trim()
                                                         + "#") != -1
                &&
                parameter.getParameter(Parameter.ALTERNATIVE_LOCATION2).indexOf(
                          Character.toString(atom.getAlternativeLocation())
                                                                      ) != -1) {
                return true;
            }
            return false;
        } else {
            return true;
        }
    }

    //--------------------------------------------------------------------------

    /**
     * Returns those pairs of potential cross-linkable atoms that have an
     * Euclidean distance smaller than the user set maximum distance.
     * @param candidates1
     *        - Set of AtomList objects that fulfill all criteria set by the
     *          user on the commandline for the first cross-linked atoms.
     * @param candidates2
     *        - Set of AtomList objects that fulfill all criteria set by the
     *          user on the commandline for the second cross-linked atoms.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return Hashtable of potential cross-linkable atoms.
     */
    private static Hashtable < Atom, AtomList > createPairsBetweenRelevantAtoms(
                                       final TreeSet < AtomList > candidates1,
                                       final TreeSet < AtomList > candidates2,
                                       final CrossLinkParameter parameter) {
        Hashtable < Atom, AtomList > pairs = new Hashtable < Atom, AtomList >();
        for (AtomList list1 : candidates1) {
            for (AtomList list2 : candidates2) {
                // An amino acid can not be self-cross-linked.
                if (MatterUtilities.equalsResidue(list1.get(0), list2.get(0))) {
                    continue;
                } else {
                    AtomList minimumDistanceAtomPair =
                            MatterUtilities.getClosestAtomPair(list1, list2);

                    double dist =
                                  Mathematics.distance(
                                    minimumDistanceAtomPair.get(0).getPoint3d(),
                                    minimumDistanceAtomPair.get(1).getPoint3d()
                                                      );

                    double errorRange = Constants.getCoordinateUncertainty(
                                                  minimumDistanceAtomPair.get(0)
                                                                          )
                                        +
                                        Constants.getCoordinateUncertainty(
                                                minimumDistanceAtomPair.get(1)
                                                                          );

                    if (dist > Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                           )
                                                 ) + errorRange
                       ) {
                        continue;
                    } else {
                        // The distance between these two amino acids would
                        // indicate that both could be cross-linked, at least in
                        // terms of their Euclidean distance.
                        // At a later stage, these potential candidates can be
                        // checked further to have a Solvent-Path distance that
                        // conforms to the length of the cross-linker.
                        Atom atom0 = minimumDistanceAtomPair.get(0);
                        Atom atom1 = minimumDistanceAtomPair.get(1);
                        AtomList associate0 = pairs.get(atom0);
                        AtomList associate1 = pairs.get(atom1);
                        if (associate0 == null && associate1 == null) {
                            associate0 = new AtomList();
                            associate0.add(atom1);
                            pairs.put(atom0, associate0);
                        } else {
                            if (associate0 == null) {
                                if (!associate1.contains(atom0)) {
                                    associate1.add(atom0);
                                    pairs.put(atom1, associate1);
                                }
                            } else if (associate1 == null) {
                                if (!associate0.contains(atom1)) {
                                    associate0.add(atom1);
                                    pairs.put(atom0, associate0);
                                }
                            }
                        }
                    }
                }
            }
        }

        return pairs;
    }

    //--------------------------------------------------------------------------

    /**
     * Returns pairs of atoms that depending on the user set parameter, contain
     * only intra, inter or all potential cross-links.
     * @param pairs
     *        - Hashtable of all pairs of atom that conform to the atom and
     *          amino acid identifiers as set by the user.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return Hashtable of intra-, inter or intra/inter atom pairs.
     */
    private static Hashtable < Atom, AtomList > fixIntraInterSelection(
                                       final Hashtable < Atom, AtomList > pairs,
                                       final CrossLinkParameter parameter) {
        Hashtable < Atom, AtomList > newPairs =
                                            new Hashtable < Atom, AtomList >();

        for (Atom atom1 : pairs.keySet()) {
            for (Atom atom2 : pairs.get(atom1)) {
                if (Boolean.parseBoolean(parameter.getParameter(
                                           Parameter.DO_INTRAMOLECULAR_DISTANCE)
                                                               )
                                         &&
                   !Boolean.parseBoolean(parameter.getParameter(
                                           Parameter.DO_INTERMOLECULAR_DISTANCE)
                                                               )
                                        ) {
                    if (atom1.getChainId() != atom2.getChainId()) {
                        continue;
                    }
                } else if (!Boolean.parseBoolean(parameter.getParameter(
                                            Parameter.DO_INTRAMOLECULAR_DISTANCE
                                                                       )
                                                )
                           &&
                           Boolean.parseBoolean(parameter.getParameter(
                                            Parameter.DO_INTERMOLECULAR_DISTANCE
                                                                      )
                                               )
                          ) {
                    if (atom1.getChainId() == atom2.getChainId()) {
                        continue;
                    }
                }
                AtomList atomList = newPairs.get(atom1);
                if (atomList == null) {
                    atomList = new AtomList();
                    atomList.add(atom2);
                } else {
                    atomList.add(atom2);
                }
                newPairs.put(atom1, atomList);
            }
        }
    return newPairs;
    }

    //--------------------------------------------------------------------------

    /**
     * Removes redundant atom pairs within homologous structures.
     * @param crossLinkList
     *        - CrossLinkList object holding all user set parameters for
     *          calculating cross-links
     */
    private static void removeRedundanciesInHomomers(
                                               final CrossLinkList crossLinkList
                                                    ) {

        // get all redundant cross links
        Hashtable < CrossLink, ArrayList < CrossLink > >
                                                   redundantCrossLinksCandidates
                      = new Hashtable < CrossLink, ArrayList < CrossLink > >();
        for (CrossLink crossLink1 : crossLinkList) {
            for (CrossLink crossLink2 : crossLinkList) {
                if (crossLink1 != crossLink2) {
                    if (crossLink1.equalsInHomolog(crossLink2)) {
                        ArrayList < CrossLink > list =
                                  redundantCrossLinksCandidates.get(crossLink1);
                        if (list == null) {
                            list = new ArrayList < CrossLink >();
                        }
                        list.add(crossLink2);
                        redundantCrossLinksCandidates.put(crossLink1, list);
                    }
                }
            }
        }

        // remove all redundant cross links except of the one with the lowest
        // Euclidean/SolventPath distance.
        for (CrossLink crossLink1 : redundantCrossLinksCandidates.keySet()) {
            CrossLink minXL = crossLink1;
            double minDist =
                (minXL.getSolventPathDistance()
                 ==
                 Double.parseDouble(Value.DISTANCE.getDefault())
                 ?
                 minXL.getEuclideanDistance() : minXL.getSolventPathDistance());
            for (CrossLink crossLink2 : redundantCrossLinksCandidates.get(
                                                                      crossLink1
                                                                         )
                ) {
                double dist2 =
                    (crossLink2.getSolventPathDistance()
                     ==
                     Double.parseDouble(Value.DISTANCE.getDefault())
                     ?
                     crossLink2.getEuclideanDistance()
                                         : crossLink2.getSolventPathDistance());

                if (minDist < dist2) {
                    crossLinkList.remove(crossLink2);
                } else {
                    crossLinkList.remove(minXL);
                    minDist = dist2;
                    minXL = crossLink2;
                }
            }
        }
    }
   //--------------------------------------------------------------------------

    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance. The cross-links
     * objects can either be generated on a local or global grid, where the
     * former is used when complex != {@code NULL} while the latter is used
     * when complex == {@code NULL} and grid != {@code NULL}.
     * @param complex
     *      - PolyPeptideList object holding all atoms of the protein.
     * @param grid
     *      - Local AtomGrid object build around one potential cross-linkable
     *        atom.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkList object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkList evaluateSolventPathDistance(
                           final PolyPeptideList complex,
                           final AtomGrid grid,
                           final CrossLinkList crossLinksByEuclideanDistance,
                           final CrossLinkParameter parameter
                                                          ) {

        CrossLinkList crossLinksBySolventPathDist = new CrossLinkList();
        double maxDist = Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                  )
                                           );

        Hashtable <Atom, AtomList> pairs =
                                         crossLinksByEuclideanDistance.toHash();

        for (Atom atom : pairs.keySet()) {
            AtomList pairedAtoms = pairs.get(atom);
            ArrayList < Path > paths = null;
            // local grid
            if (complex != null) {
                paths = CrossLinkUtilities.calculateShortestPathThroughSolvent(
                                                                    complex,
                                                                    atom,
                                                                    pairedAtoms,
                                                                    parameter
                                                                          );
            // global grid
            } else if (grid != null) {
                SolventPathDistance solvDist  = new SolventPathDistance(
                                                                    atom,
                                                                    pairedAtoms,
                                                                    grid
                                                                       );
                paths = solvDist.getShortestPath(
                                        maxDist
                                        +
                                        Constants.getCoordinateUncertainty(atom)
                                        +
                                        Constants.getCoordinateUncertainty(
                                                                     pairedAtoms
                                                                          ));

                if (CrossLinkUtilities.cleanAtomPairs(paths,
                                                      pairedAtoms,
                                                      maxDist)
                    &&
                    Boolean.parseBoolean(parameter.getParameter(
                                                            Parameter.FIND_ALL))
                                                               ) {
                    paths = new ArrayList <Path>();
                }

                if (Boolean.parseBoolean(parameter.getParameter(
                                                        Parameter.DO_GRID_OUTPUT
                                                               )
                                        )) {
                    System.out.print(grid.toString());
                }

            }
            // in case no shortest path could be found. Can happen e.g. when
            // distance file is supplied and all distances are to be found but
            // not all distances could be found.
            if (paths.size() == 0
                &&
                Boolean.parseBoolean(
                                     parameter.getParameter(Parameter.FIND_ALL))
                                    ) {
                return new CrossLinkList();
            }

            for (int i = 0; i < pairedAtoms.size(); i++) {

                CrossLink crossLink = crossLinksByEuclideanDistance.get(
                                                              atom,
                                                              pairedAtoms.get(i)
                                                                       );

                double dist = SolventPathDistance.extractTargetDistances(
                                                                    paths.get(i)
                                                                        );
                double errorRange =
                         Constants.getCoordinateUncertainty(atom)
                         +
                         Constants.getCoordinateUncertainty(pairedAtoms.get(i));

                if (dist <= maxDist + errorRange) {
                    crossLink.setSolventPathDistance(dist);
                    crossLink.setPath(paths.get(i));

                    crossLinksBySolventPathDist.add(crossLink);
                }
            }
        }

    return crossLinksBySolventPathDist;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates solvent path distances using a local grid.
     * Note, that all atom objects in atoms2 which have a distance larger than
     * maxDist to atom1 will removed from atoms2.
     * @param complex
     *      - PolyPeptideList object holding all atoms of the protein.
     * @param atom1
     *        - First protein atom to be connected by the virtual cross-linker.
     * @param atoms2
     *        - List of atoms to be cross-linked to atom1, if distance is
     *          shorter than maxDist.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return List of Path objects that form the shortest path between atom1
     *         and those atom2 objects that have a distance shorter than
     *         maxDist.
     */
    private static ArrayList < Path > calculateShortestPathThroughSolvent(
                                              final PolyPeptideList complex,
                                              final Atom atom1,
                                              final AtomList atoms2,
                                              final CrossLinkParameter parameter
                                                            ) {

        boolean doSAS = Boolean.parseBoolean(parameter.getParameter(
                                                                Parameter.DO_SAS
                                                                   )
                                            );

        double gridCellSize = Double.parseDouble(parameter.getParameter(
                                                Parameter.GRID_CELL_SIZE
                                                                       )
                                                );
        double maxDist = Double.parseDouble(parameter.getParameter(
                                              Parameter.MAXIMUM_DISTANCE
                                                                  )
                                           )
                         +
                         Constants.getCoordinateUncertainty(atom1)
                         +
                         Constants.getCoordinateUncertainty(atoms2);

        AtomGrid grid = new AtomGrid(complex.getAllAtoms(),
                                     atom1,
                                     maxDist,
                                     gridCellSize
                                    );

        // Check for solvent accessibility, if necessary.
        if (doSAS) {
            AtomList toBremoved = new AtomList();
            for (Atom atom2 : atoms2) {
                if (!CrossLinkUtilities.isAccessible(atom1, atom2, grid)) {
                    toBremoved.add(atom2);
                }
            }
            atoms2.removeAll(toBremoved);
            if (atoms2.size() == 0) {
                return new ArrayList < Path >();
            }
        }

        SolventPathDistance solvDist  = new SolventPathDistance(
                                                                atom1,
                                                                atoms2,
                                                                grid
                                                               );
        ArrayList < Path > paths = solvDist.getShortestPath(maxDist);

        if (CrossLinkUtilities.cleanAtomPairs(paths, atoms2, maxDist)
            &&
            Boolean.parseBoolean(parameter.getParameter(Parameter.FIND_ALL))) {
                paths = new ArrayList <Path>();
        }

    return paths;
    }
    //--------------------------------------------------------------------------
    /**
     * Checks whether both atoms in a cross link are within the grid and
     * solvent accessible.
     * @param atom1
     *        - First atom to be checked for accessibility
     * @param atom2
     *        - Second atom to be checked for accessibility
     * @param grid
     *        - AtomGrid object used to derive the solvent accessibility
     * @return {@code TRUE} if both atoms are accessible, {@code FALSE}
     *         otherwise.
     */
    private static boolean isAccessible(final Atom atom1,
                                        final Atom atom2,
                                        final AtomGrid grid) {

        if (grid.get(atom1) == null || grid.get(atom2) == null) {
            return false;
        }

        if (!GridUtilities.isAccessible(atom1, grid)
            ||
            !GridUtilities.isAccessible(atom2, grid)) {
            return false;
        }
        return true;
    }
    //--------------------------------------------------------------------------
    /**
     * Cleans the list of paths and atom pairs from atom pairs that have a
     * distance larger than maxDist.
     * @param paths
     *        - List of Path objects that form the shortest path between atom1
     *          and those atom2 objects that have a distance shorter than
     *          maxDist.
     * @param atoms2
     *        - List of atoms to be cross-linked to atom1, if distance is
     *          shorter than maxDist.
     * @param maxDist
     *        - double value representing the maximum distance beyond which
     *          atom pairs and path objects will be removed.
     * @return {@code TRUE} if elements in paths and atoms2 object had to be
     *         removed, {@code FALSE} otherwise.
     */
    private static boolean cleanAtomPairs(final ArrayList < Path > paths,
                                       final AtomList atoms2,
                                       final double maxDist) {
        // remove all atoms that have a distance larger than maxDist.
        ArrayList < Path > paths2beRemoved = new ArrayList < Path >();
        AtomList toBremoved = new AtomList();
        for (int i = 0; i < paths.size(); i++) {
             if (paths.get(i).size() == 1) {
                toBremoved.add(atoms2.get(i));
                paths2beRemoved.add(paths.get(i));
             } else if (SolventPathDistance.extractTargetDistances(
                                                            paths.get(i)
                                                                  )
                                                       > maxDist
                       ) {
                 toBremoved.add(atoms2.get(i));
                 paths2beRemoved.add(paths.get(i));
             }
        }
        atoms2.removeAll(toBremoved);
        paths.removeAll(paths2beRemoved);

        if (toBremoved.size() > 0 || paths2beRemoved.size() > 0) {
            return true;
        }
        return false;
    }
    //--------------------------------------------------------------------------
}
