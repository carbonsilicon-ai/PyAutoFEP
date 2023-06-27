

from rdkit import Chem
import rdkit.Chem.PropertyMol
from merge_topologies_zcc import find_mcs, constrained_two_mol
import pickle
import os
import savestate_util
import networkx
from generate_perturbation_map_atom_map_wangwei import fill_thermograph

# read the progress.pkl and modify the perturbation map and corresponding atom maps
progress_data = savestate_util.SavableState('progress.pkl')

bias = progress_data['thermograph']['last_solution']['bias']
try:
    ref_pose = progress_data['superimpose_data']['reference_pose_superimpose']
    ref_mol = ref_pose
    print('tem_mol',Chem.MolToSmiles(ref_mol))
except:
    print('the reference_pose are not supplied')
    ref_mol = bias
perturbation_map = progress_data['perturbation_map']


# alignment between molecules of each edge
lig_data_path = 'lig_data'

for edge in perturbation_map:           
    a_mol_p = '{}.mol'.format(edge[0])
    b_mol_p = '{}.mol'.format(edge[1])
    #t_mol_p = '{}.mol'.format(ref_mol)
    
    a_mol_path = os.path.join(lig_data_path, a_mol_p)
    b_mol_path = os.path.join(lig_data_path, b_mol_p)
    #t_mol_path = os.path.join(lig_data_path, t_mol_p)

    if os.path.exists(a_mol_path) and os.path.exists(b_mol_path):
        mola = rdkit.Chem.MolFromMolFile(a_mol_path, removeHs=False)
        molb = rdkit.Chem.MolFromMolFile(b_mol_path, removeHs=False)
        
        #target = rdkit.Chem.MolFromMolFile(t_mol_path, removeHs=False)
        mol_al,mol_b,_ = constrained_two_mol(mola,molb,ref_mol)
        
        mol_file_names = ('{}.mol'.format(edge[0]), '{}.mol'.format(edge[1]))
        
        aligned_lig_data_path = 'aligned_lig_data'
        if not os.path.exists(aligned_lig_data_path):
            os.makedirs(aligned_lig_data_path)
            
        for mol, mol_file_name in zip((mol_al,mol_b), mol_file_names):
            mol_file_path = os.path.join(aligned_lig_data_path, mol_file_name)
            rdkit.Chem.MolToMolFile(mol, mol_file_path)

    os.system('cp lig_data/*.itp aligned_lig_data')
            




