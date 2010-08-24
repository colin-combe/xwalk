package xwalk.crosslink;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Hashtable;

import xwalk.constants.Constants;
import xwalk.constants.Constants.ParameterSets;
import xwalk.grid.AtomGrid;
import xwalk.grid.GridUtilities;
import xwalk.grid.Path;
import xwalk.io.WriteFile;
import xwalk.io.pdb.PDBreader;
import xwalk.math.Mathematics;
import xwalk.math.Point3d;
import xwalk.math.SolventPathDistance;
import xwalk.matter.Atom;
import xwalk.matter.AtomList;
import xwalk.matter.MatterUtilities;
import xwalk.matter.pdb.AminoAcid;
import xwalk.matter.pdb.PolyPeptide;
import xwalk.matter.pdb.ProteinComplex;
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
     * @throws FileNotFoundException if PDB file could not be found.
     * @return CrossLinkSet object that holds all virtual cross-links on a
     *         protein complex.
     */
    public static CrossLinkSet getVirtualCrossLinks(
                                              final CrossLinkParameter parameter
                                                   )
                                                  throws FileNotFoundException {
        
        ProteinComplex complex = CrossLinkUtilities.getComplexCoordinates(
                                                                       parameter
                                                                         );

        AtomList allAtoms = complex.getAllAtoms();
        Point3d max = MatterUtilities.getMaximumCooridnate(allAtoms);
        Point3d min = MatterUtilities.getMinimumCooridnate(allAtoms);
        Point3d diff = new Point3d(max.getX() - min.getX(),
                                   max.getY() - min.getY(),
                                   max.getZ() - min.getZ());
        double sum = Math.pow(diff.getX(), 2)
                   + Math.pow(diff.getY(), 2)
                   + Math.pow(diff.getZ(), 2);
        
        double dim = Math.sqrt(sum);

        // if dimension is larger than maxDist
        Hashtable < Atom, AtomList > relevantAtomPairs =
                   CrossLinkUtilities.findRelevantPairs(parameter, complex);

        if (dim > Constants.MAX_PROTEIN_DIMENSION 
            ||
            Boolean.parseBoolean(parameter.getParameter(
                                                         Parameter.DO_LOCAL_GRID
                                                       )
                                )
           ) {
            return CrossLinkUtilities.calculateCrossLinksOnLocalGrid(
                    complex,
                    relevantAtomPairs,
                    parameter
                          );
        } else {
            return CrossLinkUtilities.calculateCrossLinksOnGlobalGrid(
                    complex,
                    relevantAtomPairs,
                    parameter
                           );
        }
    }
    //--------------------------------------------------------------------------
    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance using a single grid
     * object.
     * @param complex
     *      - ProteinComplex object holding all atoms of the protein.
     * @param relevantAtomPairs
     *      - Hashtable of potential cross-linkable atoms.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkSet object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkSet calculateCrossLinksOnGlobalGrid(
                           final ProteinComplex complex,
                           final Hashtable < Atom, AtomList > relevantAtomPairs,
                           final CrossLinkParameter parameter
                                                              ) {
     
        boolean doSAS = Boolean.parseBoolean(parameter.getParameter(
                                                                Parameter.DO_SAS
                                                                   )
                                            );
        CrossLinkSet xlinkSet = null;
        if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                       )
                                )
           ) {
            double solventRadius = 0;
            if (doSAS) {
                solventRadius = Double.parseDouble(parameter.getParameter(
                                                        Parameter.SOLVENT_RADIUS
                                                                         )
                                                  );
            }
            AtomGrid grid = new AtomGrid(
                    complex.getAllAtoms(),
                    Double.parseDouble(parameter.getParameter(
                                                        Parameter.GRID_CELL_SIZE
                                                             )
                                      ), 
                    solventRadius + 1, 
                    doSAS
            );

            xlinkSet = CrossLinkUtilities.crossLinkRelevantAtomPairs(
                                                              relevantAtomPairs,
                                                              grid,
                                                              parameter);

            } else {
                xlinkSet = CrossLinkUtilities.crossLinkRelevantAtomPairs(
                                                              relevantAtomPairs,
                                                              null,
                                                              parameter);
            }
        return xlinkSet;
    }
    //--------------------------------------------------------------------------
    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance using a grid object
     * for each atom1 object.
     * @param complex
     *      - ProteinComplex object holding all atoms of the protein.
     * @param relevantAtomPairs
     *      - Hashtable of potential cross-linkable atoms.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkSet object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkSet calculateCrossLinksOnLocalGrid(
                           final ProteinComplex complex,
                           final Hashtable < Atom, AtomList > relevantAtomPairs,
                           final CrossLinkParameter parameter
                                                               ) {
        CrossLinkSet xlinkSet = 
            CrossLinkUtilities.crossLinkRelevantAtomPairs(
                                                          complex,
                                                          relevantAtomPairs,
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
     * @throws FileNotFoundException if PDB file could not be found.
     * @return ProteinComplex object that holds all atom coordinates of the
     *          input file as set by the user.
     */
    private static ProteinComplex getComplexCoordinates(
                                             final CrossLinkParameter parameter)
                                                  throws FileNotFoundException {
        ProteinComplex complex = null;
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
            PDBreader reader = new PDBreader(
                                   parameter.getParameter(Parameter.INFILE_PATH)
                                            );
            complex = reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                                              );
        } else {
            PDBreader reader = new PDBreader(
                                 parameter.getParameter(Parameter.INFILE_PATH));
            ProteinComplex complex1  = reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                                                               );
            ProteinComplex complex2  = reader.getProteinComplex(
                         parameter.getParameter(Parameter.CHAIN_ID1),
                         parameter.getParameter(Parameter.ALTERNATIVE_LOCATION1)
                                                               );
            
            complex = complex1;
            complex.addAll(complex2);
        }

        CrossLinkUtilities.setRadius(complex, parameter);

    return complex;
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
     */
    private static void setRadius(final ProteinComplex complex,
                                  final CrossLinkParameter parameter) {

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
     * identifiers set by the user.
     * @param complex -
     *        Protein complex object.
     * @param parameter -
     *        CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links.
     * @return Hashtable of atom pairs that conform to the user set identifiers.
     */
    private static Hashtable < Atom, AtomList > findRelevantPairs(
                                             final CrossLinkParameter parameter,
                                             final ProteinComplex complex) {
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
        StringBuffer dataNotFoundMessage = new StringBuffer();

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
                            dataNotFoundMessage.append("WARNING: "
                                 + parameter.getParameter(Parameter.INFILE_PATH)
                                 + "\tAtom \"-aa1 "
                                 + parameter.getParameter(Parameter.ATOM_TYPE1)
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
     * @param crossLinkSet
     *        - CrossLinkSet object holding all user set parameters for
     *          calculating cross-links
     */
    private static void removeRedundanciesInHomomers(
                                                 final CrossLinkSet crossLinkSet
                                                    ) {

        // get all redundant cross links
        Hashtable < CrossLink, ArrayList < CrossLink > >
                                                   redundantCrossLinksCandidates
                      = new Hashtable < CrossLink, ArrayList < CrossLink > > ();
        for (CrossLink crossLink1 : crossLinkSet) {
            for (CrossLink crossLink2 : crossLinkSet) {
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
        crossLinkSet.removeAll(redundantCrossLinks);
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
     * that are closer than the user set maximum distance.
     * @param relevantAtomPairs
     *      - Hashtable of potential cross-linkable atoms.
     * @param grid
     *      - AtomGrid object of the protein complex molecule
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkSet object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkSet crossLinkRelevantAtomPairs(
                           final Hashtable < Atom, AtomList > relevantAtomPairs,
                           final AtomGrid grid,
                           final CrossLinkParameter parameter
                                                          ) {
        boolean doSAS = Boolean.parseBoolean(parameter.getParameter(
                                                                Parameter.DO_SAS
                                                                   )
                                            );
        
        CrossLinkSet xlinkSet = new CrossLinkSet();
        double maximumDistance = Double.parseDouble(parameter.getParameter(
                                                     Parameter.MAXIMUM_DISTANCE
                                                                          )
                                                   );
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            if (doSAS && !GridUtilities.isAccessible(atom1, grid)) {
                continue;
            }

            AtomList atoms2 = new AtomList();

            for (Atom atom2 : relevantAtomPairs.get(atom1)) {
                if (doSAS && !GridUtilities.isAccessible(atom2, grid)) {
                    continue;
                }
                
                if (Mathematics.distance(atom1.getPoint3d(), 
                                         atom2.getPoint3d()
                                        ) <= maximumDistance
                    ) {
                    atoms2.add(atom2);
                }
            }

            ArrayList < Path > paths = null;
            if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                           )
                                    )
               ) {
                SolventPathDistance solvDist  = new SolventPathDistance(
                                                                        atom1,
                                                                        atoms2,
                                                                        grid
                                                                       );
                paths = solvDist.getShortestPath(maximumDistance);
            }
            for (int i = 0; i < atoms2.size(); i++) {
                if (Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                               )
                                        )
                   ) {
                    if (paths.get(i).size() == 1) {
                        continue;
                    } else if (SolventPathDistance.extractTargetDistances(
                                                 paths.get(i)) > maximumDistance
                                                                   ) {
                        continue;
                    } else {
                        CrossLink crossLink = new CrossLink(atom1,
                                                            atoms2.get(i)
                                                           );
                        crossLink.setSolventPathDistance(
                                 SolventPathDistance.extractTargetDistances(
                                                                    paths.get(i)
                                                                           )
                                                        );
                        crossLink.setPath(paths.get(i));
                        xlinkSet.add(crossLink);
                    }
                } else {
                    CrossLink crossLink = new CrossLink(atom1, atoms2.get(i));
                    xlinkSet.add(crossLink);
                }
            }
        }

        if (Boolean.parseBoolean(parameter.getParameter(
                                                          Parameter.IS_HOMOMERIC
                                                       )
                                )
            ) {
            CrossLinkUtilities.removeRedundanciesInHomomers(xlinkSet);
        }
    return xlinkSet;
    }
   //--------------------------------------------------------------------------

    /**
     * Creates CrossLink objects between potential cross-linkable atom pairs
     * that are closer than the user set maximum distance.
     * @param complex
     *      - ProteinComplex object holding all atoms of the protein.
     * @param relevantAtomPairs
     *      - Hashtable of potential cross-linkable atoms.
     * @param parameter
     *      - CrossLinkParameter object holding all user set parameters for
     *        calculating cross-links
     * @return CrossLinkSet object that holds all virtual cross-links found
     *         between the potential cross-linkable atom pairs.
     */
    private static CrossLinkSet crossLinkRelevantAtomPairs(
                           final ProteinComplex complex,
                           final Hashtable < Atom, AtomList > relevantAtomPairs,
                           final CrossLinkParameter parameter
                                                          ) {

        CrossLinkSet xlinkSet = new CrossLinkSet();

        boolean doSpd = Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                                   )
                                            );
        
        for (Atom atom1 : relevantAtomPairs.keySet()) {
            AtomList atoms2 = relevantAtomPairs.get(atom1);
            if(atom1.getSerialNumber()==14602) {
                int h=0;
            }
            
            ArrayList < Path > paths = new ArrayList < Path > (); 
            if (doSpd) {
                paths = CrossLinkUtilities.calcShortestSolventPathDistance(
                                                                       complex,
                                                                       atom1,
                                                                       atoms2,
                                                                       parameter
                                                                          );
            }
            for (int i = 0; i < atoms2.size(); i++) {
                CrossLink crossLink = new CrossLink(atom1, atoms2.get(i));
                
                if (doSpd) {
                    double dist = SolventPathDistance.extractTargetDistances(
                                                                    paths.get(i)
                                                                            );
                    crossLink.setSolventPathDistance(dist);
                    crossLink.setPath(paths.get(i));
                }
                xlinkSet.add(crossLink);
            }
        }

        if (Boolean.parseBoolean(parameter.getParameter(
                                                          Parameter.IS_HOMOMERIC
                                                       )
                                )
            ) {
            CrossLinkUtilities.removeRedundanciesInHomomers(xlinkSet);
        }
    return xlinkSet;
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
    private static ArrayList < Path > calcShortestSolventPathDistance(
                                              final ProteinComplex complex,
                                              final Atom atom1, 
                                              AtomList atoms2,
                                              final CrossLinkParameter parameter
                                                            ) {
        
        ArrayList < Path > paths = new ArrayList < Path > (); 
        
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
        if (doSAS && !GridUtilities.isAccessible(atom1, grid)) {
            atoms2.removeAll(atoms2);
            return paths;
        }

        AtomList atoms2beRemoved = new AtomList();
        for (Atom atom2 : atoms2) {
            if (grid.get(atom2) == null) {
                atoms2beRemoved.add(atom2);
            }
            if (doSAS && !GridUtilities.isAccessible(atom2, grid)) {
                atoms2beRemoved.add(atom2);
            }
        }

        atoms2.removeAll(atoms2beRemoved);
        
        SolventPathDistance solvDist  = new SolventPathDistance(
                                                                atom1,
                                                                atoms2,
                                                                grid
                                                               );
        paths = solvDist.getShortestPath(maxDist);
        
        // remove all atoms that have a distance larger than maxDist.
        ArrayList < Path > paths2beRemoved = new ArrayList < Path > ();
        for (int i = 0; i < paths.size(); i++) {
             if (paths.get(i).size() == 1) {
                atoms2beRemoved.add(atoms2.get(i));
                paths2beRemoved.add(paths.get(i));
             } else if (SolventPathDistance.extractTargetDistances(
                                                            paths.get(i)
                                                                  )
                                                       > maxDist
                       ) {
                 atoms2beRemoved.add(atoms2.get(i));
                 paths2beRemoved.add(paths.get(i));
             }
        }
        atoms2.removeAll(atoms2beRemoved);
        paths.removeAll(paths2beRemoved);
        
        if (atoms2.size() != 0
            &&
            Boolean.parseBoolean(parameter.getParameter(
                                                        Parameter.DO_GRID_OUTPUT
                                                       )
                                )
           ) {
            System.out.print("HEADER "
                             + atom1.getSerialNumber()
                             + atom1.getResidueName()
                             + atom1.getResidueNumber()
                             + atom1.getChainId()
                             + Constants.LINE_SEPERATOR
                             + grid.toString()
                             + "END"
                             + Constants.LINE_SEPERATOR);
        }
    return paths;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the header for the distance file.
     * @param outputSolventPathHeader
     *      - boolean object indicating whether to print out the header
     *        information for solvent path distances.
     * @return String object holding the distance file header.
     */
    private static String getDistanceFileHeader(
                                           final boolean outputSolventPathHeader
                                               ) {
        StringBuffer output = new StringBuffer();
        output.append("#-----\t--------\t---------\t---------\t---\t---");
        if (outputSolventPathHeader) {
            output.append("\t--------");
        }
        output.append(Constants.LINE_SEPERATOR);
        output.append("#Index\tFileName\tResi1info\tResi2info\tSeq\tEuc");
        if (outputSolventPathHeader) {
            output.append("\tSolvPath");
        }
        output.append(Constants.LINE_SEPERATOR);
        output.append("#-----\t--------\t---------\t---------\t---\t---");
        if (outputSolventPathHeader) {
            output.append("\t--------");
        }
        output.append(Constants.LINE_SEPERATOR);
    return output.toString();
    }
    //--------------------------------------------------------------------------

    /**
     * Returns cross-link objects in the distance file format.
     * @param crossLinkSet
     *        - CrossLinkSet object holding all cross-links found for a protein
     *          complex.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return String object holding the distance file formatted output of the
     *         set of CrossLink objects.
     */
    public static String toString(final CrossLinkSet crossLinkSet,
                                  final CrossLinkParameter parameter) {
        StringBuffer output = new StringBuffer();

        // get necessary values from CrossLinkParameter object.
        boolean doVerbose = Boolean.parseBoolean(parameter.getParameter(
                                                     Parameter.DO_VERBOSE_OUTPUT
                                                                       )
                                                );
        boolean doSolventPathDist = Boolean.parseBoolean(parameter.getParameter(
                                              Parameter.DO_SOLVENT_PATH_DISTANCE
                                                                               )
                                                        );
        String infile = parameter.getParameter(Parameter.INFILE_PATH).trim();
        String fileName = new File(infile).getName().trim();
        String residueType1 = parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME1
                                                    ).trim();
        String residueType2 = parameter.getParameter(
                                              Parameter.AMINO_ACID_RESIDUE_NAME2
                                                    ).trim();
        double maxDist = Double.parseDouble(parameter.getParameter(
                                                      Parameter.MAXIMUM_DISTANCE
                                                                  ).trim());

        // Start the business

        if (doVerbose) {
            output.append(CrossLinkUtilities.getDistanceFileHeader(
                                                               doSolventPathDist
                                                                  )
                         );
        }
        if (crossLinkSet.size() != 0) {

            String crossLinkSetOutput = crossLinkSet.toString();

            for (String line : crossLinkSetOutput.split("\n")) {
                 String[] columns = line.split("\t");
                 output.append(columns[0] + "\t" + fileName);
                 for (int i = 1; i < columns.length; i++) {
                      output.append("\t" + columns[i]);
                }
                output.append(Constants.LINE_SEPERATOR);
            }
        } else {
            System.err.println("WARNING: " + fileName + "\tThere is no pair of "
                             + residueType1.replaceAll("#", "") + "\t"
                             + residueType2.replaceAll("#", "")
                             + " residues that has a distance smaller than "
                             + maxDist + ".");
        }
        return output.toString();
    }

    //--------------------------------------------------------------------------

    /**
     * Returns a PyMol script, which loads the complex, highlights all virtual
     * cross-linked atoms and draws dashed lines in-between them as indications
     * for their cross-link.
     * @param crossLinkSet
     *        - CrossLinkSet object holding all cross-links found for a protein
     *          complex.
     * @param parameter
     *        - CrossLinkParameter object holding all user set parameters for
     *          calculating cross-links
     * @return String object holding the PyMol script.
     */
    public static String outputPymolScript(final CrossLinkSet crossLinkSet,
                                           final CrossLinkParameter parameter
                                          ) {
        StringBuffer output = new StringBuffer();

        NumberFormat decFormat = new DecimalFormat("0.0");

        // get necessary values from CrossLinkParameter object.
        String infile = parameter.getParameter(Parameter.INFILE_PATH);
        String nl = Constants.LINE_SEPERATOR;

        String infileWithoutExtension = new File(infile).getName().replaceAll(
                                                                      ".pdb", ""
                                                                             );

        output.append("load " + infile + nl);
//        output.append("disable " + infileWithoutExtension + nl);
        output.append("hide everything, " + infileWithoutExtension + nl);
        output.append("set dash_radius, 1, " + infileWithoutExtension + nl);
        output.append("bg_color white" + nl);
        output.append("util.cbc" + nl);
        output.append("create het, hetatm and " + infileWithoutExtension + nl);
        output.append("show sticks, het" + nl);
        output.append("color grey, het" + nl);
        output.append("disable het" + nl);

        boolean emptyChainId = false;

        int i = 1;
        for (CrossLink crossLink : crossLinkSet) {
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

            String distName = i++ + "_";
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
            i = 1;

            WriteFile.deleteFile(pathFileName);

            for (CrossLink crossLink : crossLinkSet) {
                Atom atom1 = crossLink.getPreAtom();
                Atom atom2 = crossLink.getPostAtom();
                String distName = i + "_"
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
                           + crossLink.getPath().toString(i) + "END" + nl);
                i++;
            }
            output.append("load " + pathFileName + ", solventPaths" + nl);
            output.append("hide everything, *solvent*" + nl);
            output.append("show spheres, *solvent*" + nl);
            output.append("set sphere_scale, 0.2, *solvent*" + nl);
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
        output.append("for i in range(1," + crossLinkSet.size() + "): "
                    + "cmd.set(\"sphere_color\",\"auto\",\"*_*-*\",i)"
                    + nl
                    + "for i in range(1," + crossLinkSet.size() + "): "
                    + "cmd.set(\"dash_color\",\"auto\",\"*_*-*\",i)"
                    + nl
                    + "for i in range(1," + crossLinkSet.size() + "): "
                    + "cmd.set(\"label_color\",\"auto\",\"*_*-*\",i)"
                    + nl);
        output.append("reset" + nl);

    return output.toString();
    }
    //--------------------------------------------------------------------------


}
