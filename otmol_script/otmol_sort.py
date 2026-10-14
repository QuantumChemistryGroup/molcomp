import os
import numpy as np
import otmol as otm
from openbabel import pybel
import sys

# from otmol.experiment.otmol_neon_clusters import assignment

pathway = os.getcwd()

file_1 = sys.argv[1]
file_2 = sys.argv[2]

mol_1 = next(pybel.readfile("xyz", os.path.join(pathway, file_1)))
mol_2 = next(pybel.readfile("xyz", os.path.join(pathway, file_2)))

X_A, T_A, B_A = otm.tl.process_molecule(mol_1)
X_B, T_B, B_B = otm.tl.process_molecule(mol_2)



alpha_list = np.arange(0.0, 1.0, 0.01)[1:]
assignment, rmsd_best, alpha_best, BCI = otm.tl.molecule_alignment(
    X_A, X_B, T_A, T_B, B_A, B_B, cst_D = 0.5,
    alpha_list = alpha_list, minimize_mismatched_edges = True, reflection=True)
# print('The BCI of the alignment is {}%'.format(BCI*100))
print(f'{rmsd_best:.3f}')

