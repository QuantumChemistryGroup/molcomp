Short manual to the MolComp program:

How to run: /path/to/program/a.out file1.xyz file2.xyz [options]

Options:

--alg=kabsch : 1-to-1 comparision of xyz structures, 
no attempts to renumber or generate enantio isomers.
--alg=std : (default), orientation along principal axes,
reflections (see also --keepenantio), renumbering (according to the Hungarian algorithm),
and then running Kabsch. The structure with the lowest Kabsch RMSD is used.
--alg=adv : as --alg=std, but no mapping between sorted principal
moment of inertia of the first and second molecule. (good when
all principal moments of both molecules are similar, i.e. for
sphere-like molecules.)

--keepenantio=true : - do not use reflections that
generate enantioisomers. I.e. enantiooisomers will be treated as different
structures.
--keepenantio=false : (default) also use reflections that can generegate enantioisomers.

--verbose=false : (default) only prints out (%f) RMSD or (%f %f %f) RMSD,
 average deviation, largest deviation (see --distance).
--verbose=true : detailed output, including the fitted first molecule 
with/without renumbering.

--connect=false : using no information on the molecular topology.
--connect=true : (default) Using connectivity information when comparing
two structures. If alligned atoms of the same element (e.g. C and C) different
in their bonding information, the distance between alligned atoms (in the cost matrix) 
is multiplied by 1000 comparing to the atoms with identical bonding patterns. The
identical bonding situations are determined based on tge BFS search + graph theory.

--bindh2heavy=true : (default), upon renumeration hydrogen and other univalent atoms
are bound to corresponding heavy atoms.
--bindh2heavy=false : upon renumation hydrogens and other univalent atoms are
renumbered according to the Hungarian algorithm (+graph theory).

--parallelidn=false : (default) populate the connectivity matrix using no parallel calculations
--parallelidn=true : parallel population of the connectivity matrix (good for large molecules)

--mol2=false : (default) use the connectivity only from the first molecule
--mol2=true : use the 1st molecule connectivity info for the first molecule,
and 2nd for the 2nd.

--writetopo=false : (default) do not write the connectivity to file
--writetopo=true : connectivity info will be written to file "topo.txt".
The file can be modified to preserve the correct connectivity. 

--readtopo=false : (default) do not read the connectivity info from the file "topo.txt"
--readtopo=true : read the connectivity info from the file "topo.txt"
Both writetopo & readtopo are especially relevant for the comparison of clusters

--dist=false : (default) print only RMSD, no average and largest absolute deviations
(distances) between corresponding atoms of the alligned structures.
--dist=true : print three numbers (%f %f %f), RMSD followed by the average and
the argest absolute deviations between the alligned atoms

--distratio=false : (default) do not use our own algorithm for comparison of structures
--distratio=true : use our own algorithm to compare the structures. The following
options are relevant:
--distlr=%f : (default: 1.0)
--distsm=%f : (default: 0.75)
--distrt=%f : (default: 4.0)
The result is printed (%li %f %f). The first digit (%li) is either 0 (for duplicates) or 1 (for different structures). The second number (%f) is the average deviation between the alligned atoms. The third number (%f) is the maximum deviation between the alligned atoms. 

The two conformers are considered to be different (%li=0) if either of the two criteria are true:
a) the largest devition between the alligned atoms is larger than a thershold (--distlr=)
b) the largest deviation is larger than another threshold (--distsm=), and the ratio between the largest deviation and the average deviation is larger than --distrt. 

-----
Duplicate removal
-----

Our python scripts can be used to remove duplicates from a given set of conformers.
If your conformers are merged into one file, first split it onto different file with 
the script split_xyz_E.py:
python3 /path/to/script/split_xyz_E.py your_file.xyz. It will generate multiple xyz
files. Then remove your original file with multiple xyz ("your_file.xyz") from this folder.
Each xyz file must have energy in the second line. 
After that run another python script molcomp.py in this folder:

IMPORTANT: keep your script exactly in the same folder as your executable
a.out

python3 /path/to/script/molcomp.py %i %f %i %i %s (%s)
first argument (%i) is the number of structures for comparision in the chunk. (12-20 is OK)
second argument (%f) is criterion below which the duplicate is recognized. If last argument 
(%s) == dist or distratio, then %f 1.0 is OK.
third argument (%i) - number of iterations. If not enough warning will be printed.
fourth argument (%i) - number of permutations. From 3 to 5 is OK.
fifth argument (%s):
	rmsd-kabsch - standard Kabsch RMSD, no renumbering
	rmsd-hungarian - renumbering followed by RMSD, as in charnley
	rmsd-graph - renumbering using BFS + graph + hungarian followed by RMSD
	dist - distance based criterion
	distratio
sixth argument - (optional) "readtopo"
If you want to read your topology from file "topo.txt", 
specify the last optional argument (%s) as readtopo.
The file "topo.txt" can be obtained by running standalone:
 a.out file1.xyz file2.xyz --connect=true writetopo=true 
 for any file1.xyz and file2.xyz from your folder.
 
 File "topo.txt" can be modified to get the correct topology, very useful for clusters.
 


 

