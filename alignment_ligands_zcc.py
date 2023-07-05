

from rdkit import Chem
import rdkit.Chem.PropertyMol
from merge_topologies_zcc import find_mcs, constrained_two_mol
import pickle
import os
import savestate_util
import json
import argparse
def AlignmentMols(ref_mol_path, lig_data_path):
    progress_data = savestate_util.SavableState('progress.pkl')
    if ref_mol_path != None:
        ref_mol = rdkit.Chem.MolFromMolFile(ref_mol_path, removeHs=False)
    else:
        try:
            bias = progress_data['thermograph']['last_solution']['bias']
        except:
            bias = None
        try:
            ref_pose = progress_data['superimpose_data']['reference_pose_superimpose']
            ref_mol = ref_pose
            print('tem_mol', Chem.MolToSmiles(ref_mol))
        except:
            print('the reference_pose are not supplied')
            ref_mol = bias
    perturbation_map = progress_data['perturbation_map']

    # alignment between molecules of each edge


    for edge in perturbation_map:
        a_mol_p = '{}.mol'.format(edge[0])
        b_mol_p = '{}.mol'.format(edge[1])
        # t_mol_p = '{}.mol'.format(ref_mol)

        a_mol_path = os.path.join(lig_data_path, a_mol_p)
        b_mol_path = os.path.join(lig_data_path, b_mol_p)
        # t_mol_path = os.path.join(lig_data_path, t_mol_p)

        if os.path.exists(a_mol_path) and os.path.exists(b_mol_path):
            mola = rdkit.Chem.MolFromMolFile(a_mol_path, removeHs=False)
            molb = rdkit.Chem.MolFromMolFile(b_mol_path, removeHs=False)

            # target = rdkit.Chem.MolFromMolFile(t_mol_path, removeHs=False)
            mol_al, mol_b, _ = constrained_two_mol(mola, molb, ref_mol)

            mol_file_names = ('{}.mol'.format(edge[0]), '{}.mol'.format(edge[1]))

            aligned_lig_data_path = 'aligned_lig_data'
            if not os.path.exists(aligned_lig_data_path):
                os.makedirs(aligned_lig_data_path)

            for mol, mol_file_name in zip((mol_al, mol_b), mol_file_names):
                mol_file_path = os.path.join(aligned_lig_data_path, mol_file_name)
                rdkit.Chem.MolToMolFile(mol, mol_file_path)

        os.system(f'cp {lig_data_path}/*.itp aligned_lig_data')

def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data
def write_json(file_path):
    path_info = {
    'ref_mol_path':'lig_data/ejm_46.mol',
    'lig_data_path': 'lig_data_output'
    }
    data = json.dumps(path_info, indent=2, separators=(',', ': '))
    with open(file_path, 'w') as f:
        f.write(data)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='the path of align_ligands json file')
    parser.add_argument('--file_path', '-f', help='where is your json file?',
                        default='check_ligs_input.json')
    args = parser.parse_args()
    file_path = 'align_ligs_input.json'
    write_json(file_path)
    try:
        data = load_json(args.file_path)
    except:
        print('cannot find your json file')
        exit()
    AlignmentMols(data['ref_mol_path'], data['lig_data_path'])
            

