Expert evaluation was conducted by viewing automatically generated composite images (https://github.com/QuantumChemistryGroup/molcomp). Selecting the most informative structural orientation plays a key role in facilitating viewing. This is achieved by rotating the molecules to minimize overlapping of atoms and bonds in the sheet plane (in the xy projection).

Let's assume the current directory contains xyz conformer files and a text file DUPXYZ.txt, each line of which contains the filenames of a pair of duplicates. Then, running two bash scripts from the scripts.tar.gz archive generates a composite image all.png.
PAR.sh DUPXYZ.txt
all_png.sh *.png
These bash scripts use the Molcomp program, the xyz and align perl scripts from the scripts.tar.gz archive, and the publicly available jmol, mogrify and montage (ImageMagick), and pngquant programs.

