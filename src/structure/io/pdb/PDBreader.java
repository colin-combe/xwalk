package structure.io.pdb;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.zip.DataFormatException;

import structure.constants.Constants;
import structure.io.ReadFile;
import structure.math.Point3d;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.hetgroups.SmallMolecule;
import structure.matter.pdb.AminoAcid;
import structure.matter.pdb.PolyPeptide;
import structure.matter.pdb.ProteinComplex;


/**
 * Generic PDB reader class for reading in protein complex, protein, residue and
 * atom information from the ATOM and HETATM entry lines in a PDB file.
 * The only checking done concerns the allowed format and range of the numbers.
 * Chemical validity must be checked elsewhere.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class PDBreader {
    /**
     * ReadFile object holding the content of a PDB file.
     */
    private ReadFile pdbFile;
    /**
     * path to the PDB file to be converted into Java Objects.
     */
    private String fileName = "";
    /**
     * AtomList object, which holds all atoms in the PDB file.
     * @see #readAllAtoms()
     */
    private ArrayList < AtomList > allAtoms = new ArrayList < AtomList > ();
    //--------------------------------------------------------------------------
    /**
     * Constructor; Reads in all ATOM and HETATM entries from a PDB file.
     * @param  fileName
     *         - Path to PDB file.
     * @throws IOException if error occurs while reading infile.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    public PDBreader(final String fileName) throws IOException,
                                                   DataFormatException {
        this.setFileName(fileName);
        this.pdbFile = new ReadFile(fileName);
        this.readAllAtoms(this.pdbFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor. Reads in all ATOM and HETATM entries from a PDB file hold
     * within a BufferedReader object.
     * @param  bufferedReader
     *         - BufferedReader object holding the PDB file.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    public PDBreader(final BufferedReader bufferedReader)
                                                    throws IOException,
                                                           DataFormatException {
        this.pdbFile = new ReadFile(bufferedReader);
        this.readAllAtoms(this.pdbFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor; Reads in all ATOM and HETATM entries from a PDB file.
     * @param  readFile
     *         - ReadFile object holding the content of a PDB file.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    public PDBreader(final ReadFile readFile) throws DataFormatException {
        this.pdbFile = readFile;
        this.readAllAtoms(this.pdbFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to the PDB file.
     * @return String object holding the path to the PDB file.
     */
    public final String getFilePath() {
        return fileName;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the path to the PDB file.
     * @param filePath
     *        - String object holding the path to the PDB file.
     */
    public final void setFileName(String filePath) {
        this.fileName = filePath;
    }
    //--------------------------------------------------------------------------
    /**
     * General method to read in all ATOM and HETATM lines in a PDB file.
     * @param fileContent
     *        - List of String object holding each a line of a PDB file to be
     *          read in.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     * @see #parseAtom(String)
     */
    private void readAllAtoms(final ArrayList < String > fileContent)
                                                    throws DataFormatException {
           int i = 0;
           AtomList atoms = new AtomList();
           for (String line : fileContent) {
               i++;
               if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
                   Atom atom = this.parseAtom(line);
                   atoms.add(atom);
                   // NMR files hold often more than one conformation of the
                   // same protein. For the moment it should suffice to read
                   // in only the first model.
                   if (line.startsWith("END")) {
                       if (atoms.size() > 0) {
                           this.allAtoms.add(atoms);
                           atoms = new AtomList();
                       }
                   }
               }
           }
           if (atoms.size() > 0) {
               this.allAtoms.add(atoms);               
           }
    }
    //--------------------------------------------------------------------------
    /**
     * Reads in all information of an atom in ATOM or HETATM entry lines within
     * a PDB file. Information beyond the temperature column are ignored.
     * @param  line
     *         - String object holding ATOM or HETATM line text in the PDB file
     * @return AtomRadius object that holds all information of the ATOM or
     *         HETATM entry line.
     * @throws DataFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
     */
    private Atom parseAtom(final String line) throws DataFormatException {
        Atom atom = new Atom();
        try {
            atom.setFlag(line.substring(0,6));
            atom.setSerialNumber(Integer.parseInt(line.substring(6,11).trim()));
            atom.setName(line.substring(12,16));
            atom.setAlternativeLocation(line.charAt(16));
            atom.setResidueName(line.substring(17,20));
            atom.setChainId(line.charAt(21));

            int residueNumber = Integer.parseInt(line.substring(22,26).trim());
            atom.setResidueNumber(residueNumber);

            atom.setICode(line.charAt(26));

            double x = Double.parseDouble(line.substring(30,38).trim());
            double y = Double.parseDouble(line.substring(38,46).trim());
            double z = Double.parseDouble(line.substring(46,54).trim());
            atom.setPoint3d(new Point3d(x, y, z));

            atom.setOccupancy(Double.parseDouble(line.substring(54,60).trim()));

            double tempFact = Double.parseDouble(line.substring(60,66).trim());
            atom.setTemperatureFactor(tempFact);
        } catch (Exception e) {
            throw new DataFormatException("ERROR: " + line + "; does not seem "
                                        + "to have PDB format: " + e
                                        + Constants.LINE_SEPERATOR);
        }
    return atom;
    }
    //--------------------------------------------------------------------------
    /**
     * Method to read in all ATOM entries in a PDB file and convert all ATOM
     * entries to AminoAcidType, Protein objects and return PDB file as
     * ProteinComplex object.
     * @param chainIds
     *        - String object of all chain IDs to be used to build up the
     *          ProteinComplex
     * @param alternativeLocations
     *        - String object of all alternative locations to be used to build
     *          up the ProteinComplex
     * @return A ProteinComplex object that consists of Protein objects, which
     *         consist themselves of AminoAcidType objects, which again consists
     *         themselves of AtomRadius objects.
     */
    public final ArrayList < ProteinComplex > getProteinComplex(
                                               final String chainIds,
                                               final String alternativeLocations
                                                               ) {
        ArrayList < ProteinComplex > complexes =
                                            new ArrayList < ProteinComplex > ();
        
        for (AtomList atoms : this.allAtoms) {
            // First put all atoms of user requested protein chains into a
            // selection hashtable.
            Hashtable < Character, AtomList > selection =
                                       new Hashtable < Character, AtomList > ();
            for (Atom atom : atoms) {
                if ((atom.getFlag().equals("ATOM  "))
                    &&
                    (chainIds.indexOf(atom.getChainId()) != -1)
                    &&
                    (alternativeLocations.indexOf(
                                                   atom.getAlternativeLocation()
                                                 ) != -1)) {
                    if (selection.get(atom.getChainId()) == null) {
                        AtomList list = new AtomList();
                        list.add(atom);
                        selection.put(atom.getChainId(), list);
                    } else {
                        selection.get(atom.getChainId()).add(atom);
                    }
                }
            }
            ProteinComplex complex = new ProteinComplex();
            // Get for each protein chain all amino acids.
            int rank = 1;
            for (Enumeration < Character > e = selection.keys();
                 e.hasMoreElements();) {
                AtomList atomList = selection.get(e.nextElement());
                ArrayList < AminoAcid > aminoAcids = this.getAllAminoAcids(
                                                                        atomList
                                                                          );
                // assign rank positions to amino acids.
                for (AminoAcid aa : aminoAcids) {
                     aa.setRank(rank++);
                }
                PolyPeptide polyPeptide = new PolyPeptide(aminoAcids);
                complex.add(polyPeptide);
            }
            complex.setName(new File(this.fileName).getName());
            complexes.add(complex);
        }
       return complexes;
    }
    //--------------------------------------------------------------------------
    /**
     * Method to read in all ATOM entries in a PDB file and convert all ATOM
     * entries to AminoAcidType, Protein objects and return PDB file as
     * ProteinComplex object.
     * @return ProteinComplex object that consists of Protein objects, which
     *         consist themselves of AminoAcidType objects, which again consists
     *         themselves of AtomRadius objects.
     */
    public final ArrayList < ProteinComplex > getEntireProteinComplex() {
        // First put all atoms of user requested protein chains into a selection
        // hashtable.
        return this.getProteinComplex(
                                      Constants.ALPHANUMERIC,
                                      Constants.ALPHANUMERIC
                                     );
    }
    //--------------------------------------------------------------------------
    /**
     * Method to read in all small molecule ligands in a PDB file defined by the
     * flag HETATM.
     * @return A MolecularGroup object of single SmallMolecule objects.
     */
    public final ArrayList < ArrayList < SmallMolecule >> 
                                                        getAllSmallMolecules() {
        ArrayList < ArrayList < SmallMolecule >> smallMolecules =
                                new ArrayList < ArrayList < SmallMolecule >> ();
        for (AtomList atoms : allAtoms) {
            AtomList selection = new AtomList();
           
            for (Atom atom : atoms) {
                if (atom.getFlag().equals("HETATM")) {
                    selection.add(atom);
                }
            }
            smallMolecules.add(this.getAllSmallMolecules(selection));
        }
        return smallMolecules;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts all amino acids as AminoAcid objects from an AtomList object.
     * An amino acid is defined as a list of atoms that have an unique
     * combination of residue number and residue name.
     * @param atomList
     *        - AtomList object holding all Atom objects of a peptide/protein
     *          structure.
     * @return An array of AminoAcid objects.
     */
    private ArrayList < AminoAcid > getAllAminoAcids(final AtomList atomList) {
        ArrayList < AminoAcid > aminoAcids = new ArrayList < AminoAcid > ();
        Atom preAtom = atomList.get(0);
        AtomList residue = new AtomList();
        for (Atom atom : atomList) {
             if (!MatterUtilities.equalsResidue(atom, preAtom)) {
                aminoAcids.add(new AminoAcid(residue));
                preAtom = atom;
                residue = new AtomList();
                residue.add(atom);
            } else {
                residue.add(atom);
            }
        }
        aminoAcids.add(new AminoAcid(residue));
    return aminoAcids;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts all single molecules as Molecule objects from an AtomList
     * object. A molecule is defined as a list of atoms that have an unique
     * combination of chain Id, residue number and residue name.
     * @param  atomList
     *         - AtomList object holding all Atom objects of all small
     *           molecules.
     * @return An array of SmallMolecule objects.
     */
    private ArrayList < SmallMolecule > getAllSmallMolecules(
                                                         final AtomList atomList
                                                            ) {
        ArrayList < SmallMolecule > molecules =
                                             new ArrayList < SmallMolecule > ();
        Atom preAtom = atomList.get(0);
        AtomList mol = new AtomList();
        for (Atom atom : atomList) {
             if (!MatterUtilities.equalsResidue(atom, preAtom)) {
                molecules.add(new SmallMolecule(mol));
                preAtom = atom;
                mol = new AtomList();
                mol.add(atom);
            } else {
                mol.add(atom);
            }
        }
        molecules.add(new SmallMolecule(mol));
    return molecules;
    }
    //--------------------------------------------------------------------------

}
