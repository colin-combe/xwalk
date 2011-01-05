package mm.hydrophobicity;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.math.Mathematics;
import structure.math.Point3d;
import structure.matter.Atom;
import structure.matter.parameter.AminoAcidType;
import structure.matter.parameter.AtomType;
import structure.matter.parameter.ParameterReader;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;

/**
 * Class providing functionalities to calculate hydrophobic properties on
 * molecules and molecular assemblies.
 * @author Abdullah Kahraman
 * @version 3.0
 * @since 3.0
 */
public class Hydrophobicity {

    //--------------------------------------------------------------------------
    /**
     * PolyPeptide object holding all amino acids to which atomic XlogP
     * values will be calculated.
     */
    private PolyPeptideList polyPeptideComplex;

    //--------------------------------------------------------------------------
    /**
     * Constructor. Reads in the XlogP parameter file and assigns the values
     * to the protein atoms.
     * @param complex
     *        Protein complex object holding all amino acids to which atomic
     *        XlogP values will be calculated.
     */
    public Hydrophobicity(final PolyPeptideList complex) {
        ParameterReader reader = null;
        try {
            reader = new ParameterReader(Constants.ParameterSets.XLOGP);
        } catch (IOException e) {
            System.err.print("ERROR: While reading the XlogP parameter file"
                             + Constants.LINE_SEPERATOR);
        }
        this.polyPeptideComplex = complex;
        this.setAtomicXlogP(reader);
    }
    //--------------------------------------------------------------------------
    /**
     * Sets to all amino acids atoms in a protein their associated XlogP values.
     * @param reader
     *        ParameterReader object holding all atomic XlogP parameter values.
     * @return double value representing the sum of XlogP values for the
     *         protein complex.
     */
    private double setAtomicXlogP(final ParameterReader reader) {

        Hashtable <AminoAcidType, Hashtable <AtomType, Double>> xlogPs =
                                             reader.getXlogPparameterSet();
        double sum = 0;
        for (PolyPeptide polyPeptide : this.polyPeptideComplex) {
            for (AminoAcid aa : polyPeptide) {
                for (Atom atom : aa.getAllAtoms()) {
                    if (!atom.getElement().getSymbol().equals("H")) {
                        Hashtable<AtomType, Double> atomicXlogPs =
                                                       xlogPs.get(aa.getType());
                        double xlogP = atomicXlogPs.get(atom.getType());
                        atom.setXlogP(xlogP);
                        sum += xlogP;
                    }
                }
            }
        }
        return sum;
    }

    //--------------------------------------------------------------------------
    /**
     * Maps the hydrophobicity property of the protein on a list of atom. Note,
     * that only atoms within
     * mm.constants.Constants.PHYSICOCHEMICAL_INFLUENCE_RADIUS will be included
     * in the potential calculation.
     * @param points
     *        List of points on which potential will be calculated.
     * @return List of double values where each double value represents the
     *         XlogP potential of its associated point in the points parameter.
     */
    public final ArrayList<Double> mapHydrophobicity(
                                                 final ArrayList<Point3d> points
                                               ) {

        ArrayList<Double> xlogPpotential = new ArrayList<Double>();
        double max = mm.constants.Constants.PHYSICOCHEMICAL_INFLUENCE_RADIUS;

        // Initialize xlogPpotential list
        for (int i = 0; i < points.size(); i++) { xlogPpotential.add(0.0); }

        // get all protein residues that are within 9 Angstroem.
        Hashtable<Atom, ArrayList<Point3d>> env =
                                      new Hashtable<Atom, ArrayList<Point3d>>();

        for (PolyPeptide polyPeptide : this.polyPeptideComplex) {
            for (AminoAcid aa : polyPeptide) {
                for (Atom atom : aa.getAllAtoms()) {
                    for (Point3d point : points) {
                        double dist = Mathematics.distance(
                                                           atom.getPoint3d(),
                                                           point
                                                          );
                        if (dist <= max) {
                            if (!atom.getElement().getSymbol().startsWith(
                                                                          "H"
                                                                         )) {
                                if (!env.containsKey(atom)) {
                                    ArrayList<Point3d> coords =
                                                       new ArrayList<Point3d>();
                                    coords.add(point);
                                    env.put(atom, coords);
                                } else {
                                    env.get(atom).add(point);
                                }
                            }
                        }
                    }
                }
            }
        }

        for (Atom atom : env.keySet()) {
            ArrayList<Point3d> coords = env.get(atom);

            for (Point3d point : coords) {
                double h = atom.getXlogP();
                double dist = Mathematics.distance(atom.getPoint3d(), point);
                double s = Mathematics.sigmoidFunction(dist, max);
                int index = points.indexOf(point);
                xlogPpotential.set(index, xlogPpotential.get(index) + (h * s));
            }
        }
    return xlogPpotential;
    }
}
