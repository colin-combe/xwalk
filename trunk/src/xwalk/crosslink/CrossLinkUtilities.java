package xwalk.crosslink;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.zip.DataFormatException;

import structure.constants.Constants;
import structure.constants.Constants.ParameterSets;
import structure.grid.AtomGrid;
import structure.grid.GridUtilities;
import structure.grid.Path;
import structure.io.GzipFileReader;
import structure.io.pdb.GzipPDBreader;
import structure.io.pdb.PDBreader;
import structure.io.pdb.TarPDBreader;
import structure.math.Mathematics;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.pdb.AminoAcid;
import structure.matter.pdb.Digestion;
import structure.matter.pdb.PolyPeptide;
import structure.matter.pdb.ProteinComplex;

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
     * @param parameter -
     *        CrossLinkParameter object, holding all parameter that are
     *        necessary for the virtual cross-link calculation.
     * @throws IOException if an error occurred while reading the infile.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     * @return CrossLinkList object that holds all virtual cross-links on a
     *         protein complex.
     */
    public static CrossLinkList getVirtualCrossLinks(
                                              final CrossLinkParameter parameter
                                                   )
                                                  throws IOException,
                                                         DataFormatException {
        // get all protein complex atom coordinates of the user given inputFile.
        ArrayList < ProteinComplex > complexes = 
                                   CrossLinkUtilities.getComplexesCoordinates(
                                                                       parameter
                                                                             );
        if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                       )
                                )) {
            // output digested peptides
            for (ProteinComplex complex : complexes) {
                for (int i = 0; i < complex.size(); i++) {
                    System.err.print((i + 1) + ". "
                                    + complex.get(i).toStringOneLetterCode());
                }
            }
        }

        CrossLinkList crossLinkList = new CrossLinkList();
        // find and create virtual cross-links on the protein complexes.
        for (ProteinComplex complex : complexes) {
            // First find cross-links based on Euclidean distance.
            crossLinkList = CrossLinkUtilities.crossLinkByEuclideanDistance(
                                                                    parameter,
                                                                    complex
                                                                           );

            // If requested by the user check further for Solvent Path distance.
            if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                           ))) {
                crossLinkList = 
                           CrossLinkUtilities.calculatesSolventPathDistance(
                                                                   parameter,
                                                                   complex,
                                                                   crossLinkList
                                                                           );
            }
            // set file name of each grid to the file name of its protein
            // complex.
            for (CrossLink xlink : crossLinkList) {
                xlink.setFileName(complex.getName());
            }

        }
        // remove redundant cross-links if the complex should be labeled by the
        // user as homomeric.
        if (Boolean.parseBoolean(parameter.getParameter(Parameter.IS_HOMOMERIC))
        ) {
         CrossLinkUtilities.removeRedundanciesInHomomers(crossLinkList);
        }
        // sort list of cross-links by distance.
        crossLinkList.sort();
        
        // set indices of cross-links
        CrossLinkUtilities.setCrossLinkIndices(crossLinkList, parameter);

        return crossLinkList;
    }
    //--------------------------------------------------------------------------
    /**
     * CrossLinks all atoms in a protein complex that have an Euclidean distance
     * smaller than a user set maxDist value and if set by the user are found
     * in a distance file.
     * @param complex -
     *        Protein complex object.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return List of CrossLink object that all have a Euclidean distance
     *         smaller then the user set maxDist.
     *         and have a Euclidean distance smaller then the user set maxDist.
     * @throws IOException if an error occurred while reading the distance file.
     */
    public static CrossLinkList crossLinkByEuclideanDistance(
                                             final CrossLinkParameter parameter,
                                             final ProteinComplex complex) 
                                            throws IOException {
        Hashtable < Atom, AtomList > relevantAtomPairs =
            new Hashtable < Atom, AtomList > ();
        CrossLinkList distanceFileCrossLinks = new CrossLinkList();


        if (parameter.getParameter(Parameter.DISTANCE_FILE_PATH).equals("")) {
            relevantAtomPairs = CrossLinkUtilities.findRelevantPairs(
                                                                  parameter,
                                                                  complex
                                                                    );
        } else {
            // if 
            distanceFileCrossLinks = DistanceReader.getCrossLinks(
                        parameter.getParameter(Parameter.DISTANCE_FILE_PATH)
                                                                 );
            relevantAtomPairs = CrossLinkUtilities.extractRelevantPairs(
                                                     complex,
                                                     distanceFileCrossLinks,
                                                     parameter
                                                                        );
        }
        
        // create CrossLinks object from all relevant atom pairs.
        CrossLinkList crossLinks = new CrossLinkList();
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            for (Atom atom2 : relevantAtomPairs.get(atom1)) {
                crossLinks.add(new CrossLink(atom1, atom2));
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
     *      - ProteinComplex object holding all atoms of the protein.
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
                              final ProteinComplex complex,
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
        CrossLinkList xlinkSet = CrossLinkUtilities.evaluateSolventPathDistance(
                                                  crossLinksByEuclideanDistance,
                                                  grid,
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
     *      - ProteinComplex object holding all atoms of the protein.
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
                              final ProteinComplex complex,
                              final CrossLinkList crossLinksByEuclideanDistance,
                              final CrossLinkParameter parameter
                                                               ) {
        CrossLinkList xlinkSet = 
            CrossLinkUtilities.evaluateSolventPathDistance(
                                                  complex,
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
     * @return List of ProteinComplex objects, each holding all protein 
     *         coordinates labeled as ATOM up to an END flag or end of file.
     * @throws IOException if input file could not be read.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    private static ArrayList < ProteinComplex > getComplexesCoordinates(
                                             final CrossLinkParameter parameter
                                                                       ) 
                                                   throws IOException,
                                                           DataFormatException {

        ArrayList < PDBreader > pdbReaders =
            CrossLinkUtilities.createPDBreaders(parameter.getParameter(
                                                           Parameter.INFILE_PATH
                                                                      )
                                               );
        
        
        ArrayList < ProteinComplex > proteinComplexes =
              CrossLinkUtilities.extractProteinComplexes(pdbReaders, parameter);
        
        // assign vdW radius to protein atoms 
        for (ProteinComplex proteinComplex : proteinComplexes) {
             CrossLinkUtilities.setRadius(proteinComplex, parameter);
        }
        
        // digest protein
        if (Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_TRYPSIN_DIGEST
                                                       )
                                )
           ) {
            return CrossLinkUtilities.trypsinate(proteinComplexes, 
                                                 Boolean.parseBoolean(
                                                  parameter.getParameter(
                                                        Parameter.DO_EXPASY_RULE
                                                                        )
                                                                     )
                                                 );
        }
    return proteinComplexes;
    }
    //--------------------------------------------------------------------------
    /**
     * Creates list of PDBreader objects from the user given input file, where
     * the input file can be a list of PDB files in a tar archive (.tar),
     * compressed by GNU zip (.gz, .tar.gz, .tgz) or simply a PDB file.
     * @param infile
     *        - String object holding the path to the input file.
     * @return List of PDBreader objects each holding the content of a single
     *         PDB file.
     * @throws IOException if input file could not be read.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    public static ArrayList < PDBreader > createPDBreaders(final String infile) 
                                       throws IOException, DataFormatException {

        ArrayList < PDBreader > pdbReaders = new ArrayList < PDBreader > ();
        if (infile.endsWith(".tar.gz") || infile.endsWith(".tgz")) {
            GzipFileReader gzip = new GzipFileReader(infile);
            TarPDBreader tarPdb = new TarPDBreader(gzip.getGZIPInputStream());
            pdbReaders.addAll(tarPdb.getPDBreaders());
            
        } else if (infile.endsWith(".gz")) {
            GzipPDBreader gzipReader = new GzipPDBreader(infile);
            pdbReaders.add(gzipReader.getPDBreader());
        } else if (infile.endsWith(".tar")) {
            TarPDBreader tarPdb = new TarPDBreader(infile);
            pdbReaders.addAll(tarPdb.getPDBreaders());
        } else {
            pdbReaders.add(new PDBreader(infile));
        }
    return pdbReaders;
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
     * @return List of ProteinComplex objects holding only coordinates
     *         of atoms with user defined chain and alternative location
     *         information.
     */
    public static ArrayList < ProteinComplex > extractProteinComplexes(
                                       final ArrayList < PDBreader > pdbReaders,
                                       final CrossLinkParameter parameter
                                                                      ) {
        ArrayList < ProteinComplex > proteinComplexes = 
            new ArrayList < ProteinComplex > ();

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
                         parameter.getParameter(Parameter.CHAIN_ID1),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                                                         )
                                );
                proteinComplexes.addAll(reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID2),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION2)
                                                         )
                                );
            }
        }
        return proteinComplexes;
    }
    
    //--------------------------------------------------------------------------
    /**
     * Digest all protein components of a protein complex.
     * @param proteinComplexes
     *        - List of ProteinComplex objects to be digested.
     * @param useExpasyRules
     *        - boolean value indicating to use 
     *          <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">ExPASy Exception rules</a>
                for digestion. 
     * @return new ProteinComplex object with PolyPeptide objects holding the 
     *         digested peptide segments.
     */
    public static ArrayList < ProteinComplex > trypsinate(
                            final ArrayList < ProteinComplex > proteinComplexes,
                            final boolean useExpasyRules
                                                     ) {

        ArrayList < ProteinComplex > digestedComplexes = 
                                            new ArrayList < ProteinComplex > ();

        for (ProteinComplex proteinComplex : proteinComplexes) {
            ProteinComplex digestedComplex = new ProteinComplex();
            for (PolyPeptide protein : proteinComplex) {
                ArrayList < PolyPeptide > digest = Digestion.trypsinate(
                                                                  protein,
                                                                  useExpasyRules
                                                                   );
                digestedComplex.addAll(digest);
            }
            digestedComplex.setName(proteinComplex.getName());
            digestedComplexes.add(digestedComplex);                 
        }
        return digestedComplexes;
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
    private static void setRadius(final ProteinComplex complex,
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
                for (AminoAcid aa : protein.getAminoAcids()) {
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
                                             final ProteinComplex complex
                                                                 ) {
        ArrayList <ArrayList < AtomList >> relevantAtoms =
                                       new ArrayList < ArrayList <AtomList >>();
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
                                             final ProteinComplex complex,
                                             final CrossLinkList crossLinks,
                                             final CrossLinkParameter parameter
                                                                    ) 
                                                            throws IOException {
        
        Hashtable < Atom, AtomList > pairs = new Hashtable < Atom, AtomList >();
        
        
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
        
        AtomList uniqueCrossLinkAtomsInComplex = new AtomList();
        for (Atom atom1 : complex.getAllAtoms()) {
             for (Atom atom2 : uniqueCrossLinkAtoms) {
                  if (MatterUtilities.equalsResidue(atom1, atom2) 
                      &&
                      atom1.getName().trim().equals(atom2.getName().trim())) {
                          uniqueCrossLinkAtomsInComplex.add(atom1);
                  }
             }
        }

        if (uniqueCrossLinkAtomsInComplex.size()
            !=
            uniqueCrossLinkAtoms.size()) {
            return pairs;
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
                 System.err.println("Cross-linked atom in distance file could "
                                  + "not be found : " + xl.toString());
             }
                          
             double dist = Double.parseDouble(
                             Constants.DISTANCE_DEC_FORMAT.format(
                                                     Mathematics.distance(
                                                           preAtom.getPoint3d(),
                                                           postAtom.getPoint3d()
                                                                         )
                                                                 ));
             // check whether both cross-linked atoms in the distance file
             // have in this complex a distance smaller than maxDist.
             if (dist > Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                 )
                                          )
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
     *        
     */
    private static CrossLinkList calculatesSolventPathDistance(
                            final CrossLinkParameter parameter,
                            final ProteinComplex complex,
                            final CrossLinkList crossLinksByEuclideanDistance
                                                              ) {
        
        boolean verbose = Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                                     )
                                              );

        // if the dimension of the protein complex is larger than 150
        // Angstroem, the cross-links will be automatically calculated
        // on a local grid.
        double dim = MatterUtilities.getDimension(complex);
        if (verbose) {
            System.err.print("Protein Complex's dimension is "
            		       + "\"" + (int) dim + "\" Angstroem."
                           + Constants.LINE_SEPERATOR);
        }        
        if (dim > Constants.MAX_PROTEIN_DIMENSION 
            ||
            Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_LOCAL_GRID
                                                       )
                                )
           ) {
            if (verbose) {
                System.err.print("Solvent Path Distances are calculated on "
                               + "local grids." + Constants.LINE_SEPERATOR);
            }
            return CrossLinkUtilities.calculateCrossLinksOnLocalGrid(
                                                  complex,
                                                  crossLinksByEuclideanDistance,
                                                  parameter
                                                                    );
        } else {
            if (verbose) {
                System.err.print("Solvent Path Distances are calculated on "
                               + "global grid." + Constants.LINE_SEPERATOR);
            }
            return CrossLinkUtilities.calculateCrossLinksOnGlobalGrid(
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
     * @return ArrayList of AtomList objects that hold all atoms of amino acids
     *         that conform to the user set identifiers.
     */
    private static ArrayList < AtomList > findAllRelevantAtoms1(
                                             final CrossLinkParameter parameter,
                                             final ProteinComplex complex) {

        ArrayList < AtomList > candidates1 = new ArrayList < AtomList > ();

        for (PolyPeptide protein : complex) {
            for (AminoAcid residue : protein.getAminoAcids()) {
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
     * @return ArrayList of AtomList objects that hold all atoms of amino acids
     *         that conform to the user set identifiers.
     */
    private static ArrayList < AtomList > findAllRelevantAtoms2(
                                             final CrossLinkParameter parameter,
                                             final ProteinComplex complex
                                                               ) {

        ArrayList < AtomList > candidates2 = new ArrayList < AtomList > ();

        StringBuffer dataNotFoundMessage = new StringBuffer();

        for (PolyPeptide protein : complex) {
            for (AminoAcid residue : protein.getAminoAcids()) {
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
     *        - List of AtomList objects that fulfill all criteria set by the
     *          user on the commandline for the first cross-linked atoms.
     * @param candidates2
     *        - List of AtomList objects that fulfill all criteria set by the
     *          user on the commandline for the second cross-linked atoms.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links.
     * @return Hashtable of potential cross-linkable atoms.
     */
    private static Hashtable < Atom, AtomList > createPairsBetweenRelevantAtoms(
                                       final ArrayList < AtomList > candidates1,
                                       final ArrayList < AtomList > candidates2,
                                       final CrossLinkParameter parameter) {
        Hashtable < Atom, AtomList > pairs = new Hashtable < Atom, AtomList >();
        for (AtomList list1 : candidates1) {
            for (AtomList list2 : candidates2) {
                // An amino acid can not be self-cross-linked.
                if (MatterUtilities.equalsResidue(list1.get(0), list2.get(0))) {
                    continue;
                } else {
                    AtomList minimumDistanceAtomPair =
                            CrossLinkUtilities.getClosestAtomPair(list1, list2);
                    if (Mathematics.distance(
                                    minimumDistanceAtomPair.get(0).getPoint3d(),
                                    minimumDistanceAtomPair.get(1).getPoint3d()
                                            ) > Double.parseDouble(
                                                    parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                          )
                                                                  )
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
                                            new Hashtable < Atom, AtomList > ();

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
                      = new Hashtable < CrossLink, ArrayList < CrossLink > > ();
        for (CrossLink crossLink1 : crossLinkList) {
            for (CrossLink crossLink2 : crossLinkList) {
                if (crossLink1.equalsInHomolog(crossLink2)) {
                    ArrayList < CrossLink > list =
                                  redundantCrossLinksCandidates.get(crossLink1);
                    if (list == null) {
                        list = new ArrayList < CrossLink > ();
                    }
                    list.add(crossLink2);
                    redundantCrossLinksCandidates.put(crossLink1, list);
                }
            }
        }

        // remove all redundant cross links except of the one with the lowest
        // Euclidean/SolventPath distance.
        ArrayList < CrossLink > redundantCrossLinks =
                                                 new ArrayList < CrossLink > ();
        for (CrossLink crossLink1 : redundantCrossLinksCandidates.keySet()) {
            CrossLink minXL = crossLink1;
            double minDist = minXL.getSolventPathDistance() == -1 ?
                  minXL.getEuclideanDistance() : minXL.getSolventPathDistance();
            for (CrossLink crossLink2 : redundantCrossLinksCandidates.get(
                                                                      crossLink1
                                                                         )
                ) {
                double dist2 = crossLink2.getSolventPathDistance() == -1 ?
                                            crossLink2.getEuclideanDistance()
                                          : crossLink2.getSolventPathDistance();
                if (minDist < dist2) {
                    redundantCrossLinks.add(crossLink2);
                } else {
                    redundantCrossLinks.add(minXL);
                    minDist = dist2;
                    minXL = crossLink2;
                }
            }
        }
        crossLinkList.removeAll(redundantCrossLinks);
    }
    //--------------------------------------------------------------------------

    /**
     * Returns those two atoms that are closest in two AtomList objects.
     * @param list1
     *        - AtomList object holding a first list of atom coordinates.
     * @param list2
     *        - AtomList object holding a second list of atom coordinates.
     * @return AtomList object holding the coordinates of the two closest atoms
     *         in both atom lists.
     */
    private static AtomList getClosestAtomPair(final AtomList list1,
                                               final AtomList list2) {
        double minDist = Integer.MAX_VALUE;
        AtomList minList = new AtomList();
        minList.add(list1.get(0));
        minList.add(list2.get(0));

        for (Atom atom1 : list1) {
            for (Atom atom2 : list2) {
                double dist = Mathematics.distance(atom1.getPoint3d(),
                                                   atom2.getPoint3d()
                                                  );
                if (dist < minDist) {
                    minDist = dist;
                    minList.set(0, atom1);
                    minList.set(1, atom2);
                }
            }
        }
    return minList;
    }

    //--------------------------------------------------------------------------

    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance within a global
     * atom grid.
     * @param crossLinksByEuclideanDistance
     *      - List of CrossLinks object that have all a Euclidean distance
     *        smaller then a user set maxDist value between cross-linked atoms.
     * @param grid
     *      - AtomGrid object of the protein complex molecule
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkList object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkList evaluateSolventPathDistance(
                              final CrossLinkList crossLinksByEuclideanDistance,
                              final AtomGrid grid,
                              final CrossLinkParameter parameter
                                                          ) {
        double maximumDistance = Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                          )
                                                   );
        
        Hashtable <Atom, AtomList> pairs = 
                                         crossLinksByEuclideanDistance.toHash();

        CrossLinkList xlinkList = new CrossLinkList();
        for (Atom atom : pairs.keySet()) {
            AtomList pairedAtoms = pairs.get(atom);
            
            
            SolventPathDistance solvDist  = new SolventPathDistance(
                                                                    atom,
                                                                    pairedAtoms,
                                                                    grid
                                                                   );
            ArrayList < Path > paths = solvDist.getShortestPath(
                                                                maximumDistance
                                                               );
            for (int i = 0; i < pairedAtoms.size(); i++) {
                if (paths.get(i).size() == 1) {
                    continue;
                } else if (SolventPathDistance.extractTargetDistances(
                                                 paths.get(i)) > maximumDistance
                                                                     ) {
                    continue;
                } else {
                    CrossLink crossLink = new CrossLink(atom,
                                                        pairedAtoms.get(i)
                                                       );
                    crossLink.setSolventPathDistance(
                                 SolventPathDistance.extractTargetDistances(
                                                                    paths.get(i)
                                                                           ));
                    crossLink.setPath(paths.get(i));
                    xlinkList.add(crossLink);
                }
            }
        }
    return xlinkList;
    }
   //--------------------------------------------------------------------------

    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance within local grids.
     * @param complex
     *      - ProteinComplex object holding all atoms of the protein.
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
                           final ProteinComplex complex,
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
            ArrayList < Path > paths = 
                CrossLinkUtilities.calculateShortestPathThroughSolvent(
                                                                complex,
                                                                atom,
                                                                pairedAtoms,
                                                                parameter
                                                                  );
            for (int i = 0; i < pairedAtoms.size(); i++) {
                
                CrossLink crossLink = crossLinksByEuclideanDistance.get(
                                                              atom,
                                                              pairedAtoms.get(i)
                                                                       );
                    
                double dist = SolventPathDistance.extractTargetDistances(
                                                                    paths.get(i)
                                                                        );
                if (dist < maxDist) {
                    crossLink.setSolventPathDistance(dist);
                    crossLink.setPath(paths.get(i));
                    crossLink.setFileName(complex.getName());
                    
                    crossLinksBySolventPathDist.add(crossLink);
                }
            }
        }

    return crossLinksBySolventPathDist;
    }
    //--------------------------------------------------------------------------
    /**
     * Calculates solvent path distances. Note, that all atom objects in atoms2
     * which have a distance larger than maxDist to atom1 will removed from
     * atoms2.
     * @param complex
     *      - ProteinComplex object holding all atoms of the protein.
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
                                              final ProteinComplex complex,
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
                                           );
    
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
                return new ArrayList < Path > ();
            }
        }
        
        SolventPathDistance solvDist  = new SolventPathDistance(
                                                                atom1,
                                                                atoms2,
                                                                grid
                                                               );
        ArrayList < Path > paths = solvDist.getShortestPath(maxDist);
              
        CrossLinkUtilities.cleanAtomPairs(paths, atoms2, maxDist);
        
        if (atoms2.size() != 0
            &&
            Boolean.parseBoolean(parameter.getParameter(
                                                        Parameter.DO_GRID_OUTPUT
                                                       )
                                )
           ) {
            System.out.print(grid.toString(atom1));
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
     */
    private static void cleanAtomPairs(final ArrayList < Path > paths,
                                       final AtomList atoms2,
                                       final double maxDist) {
        // remove all atoms that have a distance larger than maxDist.
        ArrayList < Path > paths2beRemoved = new ArrayList < Path > ();
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
    }
    //--------------------------------------------------------------------------
}
