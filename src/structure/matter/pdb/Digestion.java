package structure.matter.pdb;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import structure.constants.Constants;
import structure.matter.parameter.AminoAcidType;

/**
 * Class for digesting PolyPeptide object according to a particular protease.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 *
 */
public class Digestion {

    /**
     * Constructor.
     */
    protected Digestion() {
        // prevents calls from subclass
        throw new UnsupportedOperationException();
    }
    //--------------------------------------------------------------------------

    /**
     * PolyPeptide digestion according to a trypsin.
     * Following info was extracted from: <br> 
     * <a href="http://www.expasy.ch/tools/peptidecutter/
     *          peptidecutter_special_enzymes.html">here</a>. <br>
     * <pre>
     *                  Cleavage 
     *                     ||
     *  Pn-----P4-P3-P2-P1-||-P1'-P2'-P3'-P4'------Pm
     *   |     |  |  |  |  ||   |   |   |        |
     *  Sn-----S4-S3-S2-S1-||-S1'-S2'-S3'-S4'------Sm
     *                     ||
     *  </pre>
     *  Px = Peptide position <br>
     *  Sx = Protease position <br>
     *  <br>
     *  For Trypsin: 
     *  <ol>
     *     <li> Cleavage preferable at Arg, Lys at P1,
     *          but neighbouring AA have large impact on digest, in
     *          particular
     *     </li>
     *     <li> Pro at P1' has negative influence
     *     </li>
     *     <li> Arg, Lys at P1' induces inhibition
     *     </li>
     *     <br>     
     *     Additionally following exception might occur: 
     *     <li> Pro usually blocks the action when found in position P1', 
     *          but not when Lys is in position P1 and Trp is in position P2 at
     *          the same time. This blocking of cleavage exerted by Pro in
     *          position P1' is also negligible when Arg is in position P1 and
     *          Met is in position P2 at the same time
     *     </li>
     *     <li> Lys is found in position P1 the following situation considerably
     *          block the action of trypsin:
     *          <ul>
     *              <li>Asp in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and Asp in position P1'
     *              </li>
     *              <li>Cys in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Tyr in position P1'
     *              </li>
     *         </ul>
     *     <li>Arg is in P1 and the following situations are found:
     *         <ul>
     *              <li>Arg in position P2 and His in position P1'
     *              </li>
     *              <li>Cys in position P2 and Lys in position P1'
     *              </li>
     *              <li>Arg in position P2 and Arg in positionP1'
     *              </li>
     *         </ul>
     *     </li>
     *  </ol>
     *   
     * @param protein
     *        - PolyPeptide object to be digested.
     * @param useException
     *        - boolean value indicating to include exceptions to digestion. 
     * @return List of PolyPeptide object being the peptides remaining after
     *         digestion.
     */
    public static ArrayList < PolyPeptide > trypsinate(
                                                     final PolyPeptide protein,
                                                     final boolean useException
                                                  ) {

        ArrayList < PolyPeptide > trypticPeptides =
                                               new ArrayList < PolyPeptide > ();
        
        ArrayList < AminoAcid > sequence = protein.getAminoAcids();

        ArrayList < AminoAcid > peptide = new ArrayList < AminoAcid > ();
        
        for (int i = 0; i < sequence.size(); i++) {
            peptide.add(sequence.get(i));
            
            AminoAcid p1 = sequence.get(i);
            AminoAcidType p1type = p1.getType();
            AminoAcid p2 = null;
            AminoAcidType p2type = null;

            AminoAcid p1p = null;
            AminoAcidType p1pType = null;
                        
            if (i - 1 > 0) {
                p2 = sequence.get(i - 1);
                p2type = p2.getType();

            }
            if (i + 1 < sequence.size()) {
                p1p = sequence.get(i + 1);
                p1pType = p1p.getType();
            }
            
            
            if (// check for rule 1.)
                p1type == AminoAcidType.LYSINE
                ||
                p1type == AminoAcidType.ARGININE) {

                // digestion is inhibited if following cases occur:
                boolean inhibit = false;
                        
                // check for rule 2.)
                if (p1pType != null) {
                    if (p1pType == AminoAcidType.PROLINE) {
                        inhibit = true;

                        // check for rule 4.)
                        if (useException) {
                            if (p1type == AminoAcidType.LYSINE
                                &&
                                p2type == AminoAcidType.TRYPTOPHAN) {
                                inhibit = false;
                            }
                            if (p1type == AminoAcidType.ARGININE
                                &&
                                p2type == AminoAcidType.METHIONINE) {
                                inhibit = false;
                            }
                        }
                    }
                }

                // check for rule 3.)
                if (p1pType != null) {
                    if (p1pType == AminoAcidType.LYSINE
                        ||
                        p1pType == AminoAcidType.ARGININE) {
                        inhibit = true;
                    }
                }

                if (useException) {
                    // check for rule 5.)
                    if (p1type == AminoAcidType.LYSINE) {
                        
                        if (p2type == AminoAcidType.ASPARTIC_ACID 
                            &&
                            p1pType == AminoAcidType.ASPARTIC_ACID
                           ) {
                            inhibit = true;
                        }
                        if (p2type == AminoAcidType.CYSTEINE
                            &&
                              (p1pType == AminoAcidType.ASPARTIC_ACID
                               ||
                               p1pType == AminoAcidType.HISTIDINE
                               ||
                               p1pType == AminoAcidType.TYROSINE)
                           ) {
                            inhibit = true;
                        }
                    }
                    // check for rule 6.)
                    if (p1type == AminoAcidType.ARGININE) {

                        if (p2type == AminoAcidType.ARGININE 
                            &&
                            p1pType == AminoAcidType.HISTIDINE
                           ) {
                            inhibit = true;
                        }
                        if (p2type == AminoAcidType.CYSTEINE 
                            &&
                            p1pType == AminoAcidType.LYSINE
                           ) {
                            inhibit = true;
                        }
                        if (p2type == AminoAcidType.ARGININE 
                            &&
                            p1pType == AminoAcidType.ARGININE
                        ) {
                            inhibit = true;
                        }
                    }
                }
                
                if (!inhibit || i == sequence.size() - 1) {
                    if (peptide.size()
                        >=
                        Constants.MIN_PEPTIDE_LENGTH
                            &&
                        peptide.size()
                        <=
                        Constants.MAX_PEPTIDE_LENGTH
                       ) {
                        trypticPeptides.add(new PolyPeptide(peptide));
                    }
                    peptide = new ArrayList < AminoAcid > ();
                }
            }
        }

        Collections.sort(trypticPeptides, new Comparator<PolyPeptide>() {

            public int compare(final PolyPeptide p1, final PolyPeptide p2) {
                if (p1.getAminoAcids().size() < p2.getAminoAcids().size()) {
                    return 1;
                }
                if (p1.getAminoAcids().size() > p2.getAminoAcids().size()) {
                    return -1;
                }
                return 0;
            } });
    return trypticPeptides;
    }

}
