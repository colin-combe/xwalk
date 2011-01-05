Xwalk (version 1.0)
-------------------

INTRODUCTION
------------
Chemical cross-linking of proteins or protein complexes and the mass spectrometry based localization of the cross-linked amino acids is a powerful method for generating distance restraints on the substrate’s topology. Here we introduce the algorithm Xwalk for predicting and validating these cross-links on existing protein structures. Xwalk calculates and displays non-linear distances between chemically cross-linked amino acids on protein surfaces, while mimicking the flexibility and non-linearity of cross-linker molecules. It returns a “Solvent Accessible Surface Distance” (SASD), which corresponds to the length of the shortest path between two amino acids, where the path leads through solvent occupied space without penetrating the protein surface.


TEST
-----
Test examples to execute Xwalk can be found in the test subdirectory.


OUTPUT
------
Distance information will be printed out to the STDOUT channel or via -out to a local file in the following tab delimeted format:
...
1	1brs.pdb	LYS-108-A-CB	LYS-98-B-CB	100	10.7	13.9	TTDHYQTFTKIR-ILYSSDWLIYKTTDHYQTFTK
...
1st column: index
2nd column: filename
3rd column: PDB information of the 1st amino acid in the format: PDBThreeLetterCode-ResidueId-ChainId-AtomName
4th column: PDB information of the 2nd amino acid in the format: PDBThreeLetterCode-ResidueId-ChainId-AtomName
5th column: Sequence distance within the PDB file as calculated from the residue identifiers of the 1st and 2nd amino acid.
5th column: Euclidean distacen between 1st and 2nd amino acid.
6th column: SASD between  1st and 2nd amino acid.
7th column: shortest tryptic peptide sequence.

Setting -p -out xxx.pml on the commandline will write a PyMOL script into the file xxx.pml, which can be load into the molecular viewer PyMOL to visualise the SASD paths (see NOTES).


ALGORITHM
---------
The SASD is calculated using a grid and the breadth-first search algorithm to search within the grid for the shortest path between two points on the protein surface using following algorithm:
1.)	Read in input data	a.	xyz.pdb, spatial coordinates of a protein or protein complex in PDB format.	b.	maxdist: maximum distance of the path (i.e. the length of the cross-linker)	c.	listXL, a list of experimentally determined cross-linked lysine residues.2.)	Remove all non-protein atoms in xyz.pdb and assign protein atoms a van der Waals radius sum of SURFNET atom radii + solvent radius3.)	Select a random lysine pair AAab from listXL, 4.)	Check Euclidean distance (Euc) of AAab. Continue, if Euc ≤ maxdist, disregard otherwise and go back to 3.) 5.)	Generate a grid of size maxdist and grid spacing 1Å centered at AAa6.)	Set Integer.MAX_VALUE as distance for all grid cells and label grid cells as residing in the	a.	protein	b.	solvent	c.	boundary between protein and solvent7.)	Label grid cells residing in AAab as solvent8.)	Set distance dist = 0.0 for central grid cell of AAa and store grid cell in the active list listactive 9.)	Start breadth-first search. Iterate through listactive 	a.	Check that grid cell i is labeled as solvent 	b.	Find all immediate neighbors listneighbour	c.	Iterate through listneighbour		i.	Check that grid cell j is labeled as solvent		ii.	Compute new distance for grid cell j as the sum of the distance in grid cell i and the Euclidean distance between grid cell i and j 		iii.	If distance sum in 8.c.ii is smaller than the current distance in grid cell j, store the distance sum as new distance for grid cell j and add grid cell j to the new active list listnew_active,		iv.	Break up iteration of 8.c.) if grid cell j == central grid of AAb10.)	Go back to step 9.) with listactive = listnew_active


NOTES
-----
- As the SASD is based on a grid calculation, the default heap size of the Java VM with 64MB is likely to be too small. You can increase the heap size with "java -Xmx512m"
- You can obtain PyMOL for free at the webpage: http://pymol.org/
- Beware that in order for PyMOL to recognize the script file, the file must have the name ending .pml
- You can load the script into PyMOL directly at the startup of PyMOL, i.e.
$> pymol 1TAB_c.pml


CONTACT
-------
abdullah@imsb.biol.ethz.ch


VERSION HISTORY
---------------
