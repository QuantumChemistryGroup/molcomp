# Molcomp
Conformer deduplication, structural superposition

## Prerequisites
1) Linux operational system 
2) Python 3 (for scripts)

## How-to-use
### Ensemble deduplication
1. Download the Python3 script ```molcomp_11.py``` from folder *molcomp_script*
2. Download the molcomp binary from folder *molcomp_bin*. Alternatively, compile it yourself,
see folder *molcomp_source*.
3. Move both ```molcomp_11.py``` and ```molcomp``` to any folder on your computer
4. Go to the folder with conformers

   ```cd /folder/with/conformers```

5. Run python ```molcomp_11.py``` script to deduplicate your ensemble:

```python3 /path/to/your/script/molcomp_11.py angle``` (in case when all conformers
are given as separate xyz files)

```python3 /path/to/your/script/molcomp_11.py file.xyz angle``` (in case when all conformers
are given in one file)

6. If modified topology file (```topo.txt```) has to be read, add 'readtopo' option as the last argument

```python3 /path/to/your/script/molcomp_11.py file.xyz angle readtopo```

7. The resulting file ```DUPXYZ-angle.txt``` contains conformer pairs determined as duplicates
8. The resulting file ```dup-angle.txt``` contains structural characterization of a revealed conformer pair:

```gfn2_conformers_00000040.xyz gfn2_conformers_00000041.xyz 0.012 0.029 0.396 0.0```

Explanation: ```0.012``` - RMSD difference between the two conformers in Angstroem, ```0.029``` - maximum distance between
mapped atom pair in the two alligned structures (*MaxDist*), ```0.396``` - maximum difference in dihedral angle in the two conformers, 
(*MaxTors*), ```0.0``` - energy differnce (kcal/mol) between two structures (if energies are provided, otherwise - repulsion energy).

9. The resulting folder (```unique-angle```) or file (```filename-unique-angle.xyz```) contains only unique conformers
### Characterization of conformer pair (Molcomp standalone use)
Molcomp standard usage on a pair of conformers:

```/path/to/molcomp file1.xyz file2.xyz ```
Result: ```1.281062``` - RMSD between two structures
## Molcomp short manual
```/path/to/molcomp file1.xyz file2.xyz {options}```

{options}: ```--compare=angle```: prints out *RMSD*, *MaxDist*, and *MaxTors* 

{options}: ```--print=true```: print the coordinates of the first molecule providing the best fit with the second to file ```file1-mod.xyz``` (default: ```--print=false```)

{options}: ```--alg=kabsch```: standard Kabsch algorithm; ```--alg=std```: standard algorithm (orientation along inertia axes, Hungarian renumbering) (default), ```--alg=adv```: advanced algorithm, more computationally expensive, apply for sphere-like molecules

{options}: ```--connect=true```: build connectivity matrix based on covalent radii & use this info for renumering (default); ```connect=false```: do not try to build connectivity matrix, renumering is only pure Hungarian algorithm followed by Kabsch

{options}: ```--verbose=true```: prints out verbose info for debugging; ```--verbose=false``` no verbose printing (default);

{options}: ```--writetopo=true```: write topology to file *topo.txt*. Then this file can be modified to ensure right topology. ```--writetopo=false```: do not write topology file. (default)

{}





> [!NOTE]
> For any inquires, questions, suggestions, etc. do not hesitate to contact us: Yury.Minenkov"at"gmail.com (replace "at" with @)

> **When using this code please cite the following publications:**
> 1) "Molcomp: a software tool for the structural characterization and deduplication of conformers", A.M. Genaev, A.D. Moshchenkov, A.S. Ryzhako, A.A. Otlyotov, and Y. Minenkov Submitted, 2026 [**Generic reference**]
