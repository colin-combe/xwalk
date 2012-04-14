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
 */

package mm.hydrophobicity;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import structure.constants.Constants;
import structure.math.Mathematics;
import structure.math.Point3f;
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
 * @version 0.1
 * @since 0.1
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
            System.err.print(e.getMessage() + Constants.LINE_SEPERATOR
                          + "ERROR: While reading the XlogP parameter file"
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
     * @return float value representing the sum of XlogP values for the
     *         protein complex.
     */
    private float setAtomicXlogP(final ParameterReader reader) {

        Hashtable <AminoAcidType, Hashtable <AtomType, Float>> xlogPs =
                                             reader.getXlogPparameterSet();
        float sum = 0;
        for (PolyPeptide polyPeptide : this.polyPeptideComplex) {
            for (AminoAcid aa : polyPeptide) {
                for (Atom atom : aa.getAllAtoms()) {
                    if (!atom.getElement().getSymbol().equals("H")) {
                        Hashtable<AtomType, Float> atomicXlogPs =
                                                       xlogPs.get(aa.getType());
                        float xlogP = atomicXlogPs.get(atom.getType());
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
     * @return List of float values where each float value represents the
     *         XlogP potential of its associated point in the points parameter.
     */
    public final ArrayList<Float> mapHydrophobicity(
                                                 final ArrayList<Point3f> points
                                               ) {

        ArrayList<Float> xlogPpotential = new ArrayList<Float>();
        float max = mm.constants.Constants.PHYSICOCHEMICAL_INFLUENCE_RADIUS;

        // Initialize xlogPpotential list
        for (int i = 0; i < points.size(); i++) { xlogPpotential.add(0.0f); }

        // get all protein residues that are within 9 Angstroem.
        Hashtable<Atom, ArrayList<Point3f>> env =
                                      new Hashtable<Atom, ArrayList<Point3f>>();

        for (PolyPeptide polyPeptide : this.polyPeptideComplex) {
            for (AminoAcid aa : polyPeptide) {
                for (Atom atom : aa.getAllAtoms()) {
                    for (Point3f point : points) {
                        float dist = Mathematics.distance(
                                                           atom.getXYZ(),
                                                           point
                                                          );
                        if (dist <= max) {
                            if (!atom.getElement().getSymbol().startsWith(
                                                                          "H"
                                                                         )) {
                                if (!env.containsKey(atom)) {
                                    ArrayList<Point3f> coords =
                                                       new ArrayList<Point3f>();
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
            ArrayList<Point3f> coords = env.get(atom);

            for (Point3f point : coords) {
                float h = atom.getXlogP();
                float dist = Mathematics.distance(atom.getXYZ(), point);
                double s = Mathematics.sigmoidFunction(dist, max);
                int index = points.indexOf(point);
                xlogPpotential.set(index, (float) (xlogPpotential.get(index)
                                           + (h * s)));
            }
        }
    return xlogPpotential;
    }
}
