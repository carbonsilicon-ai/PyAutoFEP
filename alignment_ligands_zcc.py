from rdkit import Chem
import rdkit.Chem.PropertyMol
from merge_topologies_zcc import constrained_two_mol
import os
import savestate_util
import json
import argparse


def AlignmentMols(progress_file, ref_mol_path, lig_data_path, aligned_lig_data_path):
    progress_data = savestate_util.SavableState(progress_file)
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

    if os.path.exists(aligned_lig_data_path):
        print('Output aligned_lig_data_path Pathway already existed ')
        exit()
    os.makedirs(aligned_lig_data_path)

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

            mol_file_names = ('{}.mol'.format(
                edge[0]), '{}.mol'.format(edge[1]))

            for mol, mol_file_name in zip((mol_al, mol_b), mol_file_names):
                mol_file_path = os.path.join(
                    aligned_lig_data_path, mol_file_name)
                rdkit.Chem.MolToMolFile(mol, mol_file_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='the path of align_ligands json file')
    parser.add_argument('--file_path', '-f', help='where is your json file?',
                        default='check_ligs_input.json')
    args = parser.parse_args()

    try:
        with open(args.file_path, 'r') as f:
            data = json.load(f)
    except:
        print('cannot find your check_ligs_input json file')
        exit()

    AlignmentMols(data['progress_file'], data['ref_mol_path'],
                  data['lig_data_path'], data['aligned_lig_data_path'])
