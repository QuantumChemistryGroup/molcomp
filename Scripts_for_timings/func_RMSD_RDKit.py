#!/usr/bin/python
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

def robust_legacy_connectivity(mol, tolerance=0.4):
    """
    Determines connectivity using covalent radii of specific elements 
    to build an identical topology for XYZ files.
    """
    emol = Chem.EditableMol(mol)
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    pt = Chem.GetPeriodicTable()
    
    # Pre-calculate covalent radii for all atoms
    radii = [pt.GetRcovalent(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in range(num_atoms)]
    
    for i in range(num_atoms):
        pos_i = conf.GetAtomPosition(i)
        for j in range(i + 1, num_atoms):
            dist = pos_i.Distance(conf.GetAtomPosition(j))
            max_dist = radii[i] + radii[j] + tolerance
            
            if dist <= max_dist:
                emol.AddBond(i, j, Chem.BondType.SINGLE)
                
    return emol.GetMol()


def align_xyz_files(file1_path, file2_path, reflect=True, print_output=True):
    """
    Aligns two XYZ files, calculates the minimal RMSD, applies the 
    transformation, and optionally prints the aligned coordinates.
    
    Parameters:
    file1_path (str): Path to the mobile/probe XYZ file
    file2_path (str): Path to the reference/target XYZ file
    reflect (bool): Whether to allow reflection/inversion transformations
    print_output (bool): Whether to print the aligned coordinates to stdout
    
    Returns:
    tuple: (best_rmsd, transformed_mol1)
    """
    # Read first XYZ file
    with open(file1_path, 'r') as f:
        mol1 = Chem.MolFromXYZBlock(f.read())

    # Read second XYZ file
    with open(file2_path, 'r') as f:
        mol2 = Chem.MolFromXYZBlock(f.read())

    if mol1.GetNumAtoms() != mol2.GetNumAtoms():
        raise ValueError(f"Atom count mismatch: {file1_path} ({mol1.GetNumAtoms()}) vs {file2_path} ({mol2.GetNumAtoms()})")

    # Build topologies
    mol1_bonded = robust_legacy_connectivity(mol1)
    mol2_bonded = robust_legacy_connectivity(mol2)
    res_reflect = rdMolAlign.GetBestAlignmentTransform(mol1_bonded, mol2_bonded, reflect=True)
    res_no_reflect = rdMolAlign.GetBestAlignmentTransform(mol1_bonded, mol2_bonded, reflect=False)

    # 1. Get alignment transformation
    if res_no_reflect[0] < res_reflect[0]:
        best_rmsd, trans_matrix, best_match = res_no_reflect
    else:
        best_rmsd, trans_matrix, best_match = res_reflect


    # 2. Apply the structural transformation to mol1_bonded
    AllChem.TransformMol(mol1_bonded, trans_matrix)

    # 3. Handle optional console printout
    if print_output:
        print(f"Files: {file1_path} -> {file2_path}")
        print(f"Best Match RMSD: {best_rmsd:.4f} Å\n")
        print("M2_Idx | Element |   X    |   Y    |   Z   ")
        print("-" * 42)
        
        sorted_mapping = sorted(best_match, key=lambda x: x[1])
        conf = mol1_bonded.GetConformer()
        
        for prb_idx, ref_idx in sorted_mapping:
            atom = mol1_bonded.GetAtomWithIdx(prb_idx)
            pos = conf.GetAtomPosition(prb_idx)
            symbol = atom.GetSymbol()
            print(f" {ref_idx:<5} | {symbol:<7} | {pos.x:6.3f} | {pos.y:6.3f} | {pos.z:6.3f}")
            
    return best_rmsd, mol1_bonded


# --- Allows execution via command line arguments ---
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <file1.xyz> <file2.xyz>")
        sys.exit(1)
        
    align_xyz_files(sys.argv[1], sys.argv[2])
