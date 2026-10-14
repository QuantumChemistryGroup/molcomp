#!/usr/bin/python
import sys, string
from math import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

def robust_legacy_connectivity(mol, tolerance=0.4):
    """
    Simulates rdDetermineBonds by determining connectivity 
    using the actual covalent radii of the specific elements.
    """
    emol = Chem.EditableMol(mol)
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    pt = Chem.GetPeriodicTable()
    
    # Pre-calculate covalent radii for all atoms in the molecule
    radii = [pt.GetRcovalent(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in range(num_atoms)]
    
    for i in range(num_atoms):
        pos_i = conf.GetAtomPosition(i)
        for j in range(i + 1, num_atoms):
            dist = pos_i.Distance(conf.GetAtomPosition(j))
            
            # Max allowable bonding distance = sum of covalent radii + tolerance
            max_dist = radii[i] + radii[j] + tolerance
            
            if dist <= max_dist:
                emol.AddBond(i, j, Chem.BondType.SINGLE)
                
    return emol.GetMol()


# main program starts here

# first xyz
f = open (sys.argv[1], 'r')
mol = f.read()
mol1 = Chem.MolFromXYZBlock(mol)
f.close()
# second xyz
f = open (sys.argv[2], 'r')
mol = f.read()
mol2 = Chem.MolFromXYZBlock(mol)
f.close()

if mol1.GetNumAtoms() != mol2.GetNumAtoms():
  raise ValueError("The two XYZ files have a different number of atoms.")

'''
best_rmsd = AllChem.GetBestRMS(mol1, mol2)
print (best_rmsd)
'''
# 2. Automatically assign geometric connectivity (all bonds set to SINGLE)
# 2. Build identical covalent topologies
mol1_bonded = robust_legacy_connectivity(mol1)
mol2_bonded = robust_legacy_connectivity(mol2)


res_reflect = rdMolAlign.GetBestAlignmentTransform(mol1_bonded, mol2_bonded, reflect=True)
res_no_reflect = rdMolAlign.GetBestAlignmentTransform(mol1_bonded, mol2_bonded, reflect=False)

# 1. Get alignment transformation
if res_no_reflect[0] < res_reflect[0]:
	best_rmsd, trans_matrix, best_match = res_no_reflect
else:
	best_rmsd, trans_matrix, best_match = res_reflect

#best_rmsd, _, _ = rdMolAlign.GetBestAlignmentTransform(
#    mol1_bonded, 
#    mol2_bonded, 
#    reflect=True
#)

#best_rmsd = AllChem.GetBestRMS(mol1_bonded, mol2_bonded)
print (best_rmsd)
'''
# 1. Get the RMSD, the 3D Transformation Matrix, and the best symmetry Map
# best_match will be a tuple of pairs: ((prb_idx0, ref_idx0), (prb_idx1, ref_idx1), ...)
best_rmsd, trans_matrix, best_match = AllChem.GetBestAlignmentTransform(mol1_bonded, mol2_bonded)
print(f"Best Match RMSD: {best_rmsd:.4f} Å\n")

# 2. Physically apply the structural transformation to mol1_bonded
AllChem.TransformMol(mol1_bonded, trans_matrix)

# 3. Sort the mapping by the Reference Molecule (mol2) index 
# This ensures that row i in our printout maps directly to atom i in mol2.
sorted_mapping = sorted(best_match, key=lambda x: x[1])

conf = mol1_bonded.GetConformer()

print("M2_Idx | Element |   X    |   Y    |   Z   ")
print("-" * 42)
for prb_idx, ref_idx in sorted_mapping:
    atom = mol1_bonded.GetAtomWithIdx(prb_idx)
    pos = conf.GetAtomPosition(prb_idx)
    symbol = atom.GetSymbol()
    
    # ref_idx shows which row/atom of the second XYZ file this coordinate pairs with
    print(f" {symbol:<7}  {pos.x:6.3f}  {pos.y:6.3f}  {pos.z:6.3f}")
'''
