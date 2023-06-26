from rdkit import Chem
from merge_topologies import constrained_two_mol
a_mol_p = 'mol1.mol'
b_mol_p = 'mol2.mol'

mola = Chem.MolFromMolFile(a_mol_p,removeHs=False)
molb = Chem.MolFromMolFile(b_mol_p,removeHs=False)
target = mola
constrained_two_mol(mola,molb,target)


