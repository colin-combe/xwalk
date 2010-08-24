package xwalk.io.pdb;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.zip.DataFormatException;

import xwalk.constants.Constants;
import xwalk.math.Point3d;
import xwalk.matter.Atom;
import xwalk.matter.AtomList;
import xwalk.matter.MatterUtilities;
import xwalk.matter.hetgroups.SmallMolecule;
import xwalk.matter.pdb.AminoAcid;
import xwalk.matter.pdb.PolyPeptide;
import xwalk.matter.pdb.ProteinComplex;

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
     * PDB file will be loaded into a BufferedReader object.
     */
    private BufferedReader reader;
    /**
     * path to the PDB file to be converted into Java Objects.
     */
    private String filePath;
    /**
     * AtomList object, which holds all atoms in the PDB file.
     * @see #readAllAtoms()
     */
    private AtomList allAtoms;
    //--------------------------------------------------------------------------
    /**
     * Constructor; Reads in all ATOM and HETATM entries from a PDB file.
     * @param  fileName
     *         - Path to PDB file.
     * @throws FileNotFoundException if file could not be found.
     */
    public PDBreader(final String fileName) throws FileNotFoundException {
        this.filePath = fileName;
        this.allAtoms = new AtomList();
        // check whether filename is within jar package or is user given
        // directory
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(
                                                                        fileName
                                                                              );
        if (ins != null) {
            this.reader = new BufferedReader(new InputStreamReader(ins));
        } else {
            this.reader = new BufferedReader(new FileReader(fileName));
           }
        this.readAllAtoms();
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor. Reads in all ATOM and HETATM entries from a PDB file hold
     * within a BufferedReader object.
     * @param  bufferedReader
     *         - BufferedReader object holding the PDB file.
     */
    public PDBreader(final BufferedReader bufferedReader) {
        this.reader = bufferedReader;
        this.readAllAtoms();
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
     * General method to read in all ATOM and HETATM lines in a PDB file.
     * @return {@code TRUE} if file reading is successful, {@code FALSE} if
     *         IOExeception or file does not conform to PDB standards.
     * @see #parseAtom(String)
     */
    private boolean readAllAtoms() {
           String line;
           int i = 0;
           try {
               while ((line = this.reader.readLine()) != null) {
                   i++;
                   if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
                       Atom atom = this.parseAtom(line);
                       this.allAtoms.add(atom);
                       // NMR files hold often more than one conformation of the
                       // same protein. For the moment it should suffice to read
                       // in only the first model.
                       if (line.startsWith("ENDMDL")
                           &&
                          (this.allAtoms.size() > 10)) {
                           break;
                       }
                   }
               }
           } catch (IOException e) {
               System.err.print("ERROR while reading line " + i + " in file \""
                              + this.filePath + Constants.LINE_SEPERATOR + e
                              + Constants.LINE_SEPERATOR);
               return false;
           } catch (DataFormatException e) {
               System.err.print("ERROR in line " + i + Constants.LINE_SEPERATOR
                               + e + Constants.LINE_SEPERATOR);
               return false;
           }
           return true;
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
    public final ProteinComplex getProteinComplex(
                                            final String chainIds,
                                            final String alternativeLocations
                                                 ) {
        // First put all atoms of user requested protein chains into a selection
        // hashtable.
        Hashtable < Character, AtomList > selection =
                                       new Hashtable < Character, AtomList > ();
        for (Atom atom : this.allAtoms) {
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
       return complex;
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
    public final ProteinComplex getEntireProteinComplex() {
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
    public final SmallMolecule[] getAllSmallMolecules() {
           AtomList selection = new AtomList();
           for (Atom atom : this.allAtoms) {
                if (atom.getFlag().equals("HETATM")) {
                   selection.add(atom);
               }
           }
           return this.getAllSmallMolecules(selection);
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
    private SmallMolecule[] getAllSmallMolecules(final AtomList atomList) {
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
    return (SmallMolecule[]) molecules.toArray();
    }
}
