import os
import sys
import numpy as np
import otmol as otm
from openbabel import pybel

def compute_otmol_rmsd(file_1, file_2):
    """
    Computes the optimal otmol molecular alignment using just the filenames.
    
    Parameters:
    file_1 (str): Filename of the first XYZ file (e.g., 'mol1.xyz').
    file_2 (str): Filename of the second XYZ file (e.g., 'mol2.xyz').
    
    Returns:
    float: The minimum calculated RMSD value.
    """
    # Automatically convert filenames to absolute system paths
    file_1_path = os.path.abspath(file_1)
    file_2_path = os.path.abspath(file_2)

    # Load the molecules using openbabel
    mol_1 = next(pybel.readfile("xyz", file_1_path))
    mol_2 = next(pybel.readfile("xyz", file_2_path))

    # Process molecular data using otmol
    X_A, T_A, B_A = otm.tl.process_molecule(mol_1)
    X_B, T_B, B_B = otm.tl.process_molecule(mol_2)

    # Set up the scaling parameter list
    alpha_list = np.arange(0.0, 1.0, 0.01)[1:]
    
    # Run the alignment routine
    assignment, rmsd_best, alpha_best, BCI = otm.tl.molecule_alignment(
        X_A, X_B, T_A, T_B, B_A, B_B, 
        cst_D=0.5,
        alpha_list=alpha_list, 
        minimize_mismatched_edges=True, 
        reflection=True
    )
    
    return rmsd_best


# --- Command-line execution block ---
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <file1.xyz> <file2.xyz>")
        sys.exit(1)
        
    rmsd = compute_otmol_rmsd(sys.argv[1], sys.argv[2])
    print(f'{rmsd:.3f}')
