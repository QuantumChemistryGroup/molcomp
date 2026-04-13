# Molcomp
Conformer deduplication, structural superposition

## Prerequisites
1) Linux operational system 
2) Python 3 (for scripts)

## How-to-use
### Ensemble deduplication
1. Download the Python3 script molcomp_11.py from folder molcomp_script
2. Download the molcomp binary from folder molcomp_bin. Alternatively, compile it yourself,
see folder molcomp_source.
3. Move both molcomp_11.py and molcomp to any folder on your computer
4. Go to the folder with conformers ```cd /folder/with/conformers```
5. Run python molcomp_11.py script to deduplicate your ensemble:
```python3 /path/to/your/script/molcomp_11.py angle``` (in case when all conformers
are given as separate xyz files)
```python3 /path/to/your/script/molcomp_11.py file.xyz angle``` (in case when all conformers
are given in one file)
6. 
### Characterization of conformer pair (Molcomp standalone use)

## Molcomp short manual

> [!NOTE]
> For any inquires, questions, suggestions, etc. do not hesitate to contact us: Yury.Minenkov"at"gmail.com (replace "at" with @)

> **When using this code please cite the following publications:**
> 1) "Molcomp: a software tool for the structural characterization and deduplication of conformers", A.M. Genaev, A.D. Moshchenkov, A.S. Ryzhako, A.A. Otlyotov, and Y. Minenkov Submitted, 2026 [**Generic reference**]
