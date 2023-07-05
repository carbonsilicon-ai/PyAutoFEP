import os
from rdkit import Chem
from rdkit.Chem import rdFMCS
# from merge_topologies_zcc import find_mcs
from itertools import combinations
import pandas as pd
import json
import argparse


def GetMCSInfo(mol_list):
    print('len_mol_list:', mol_list)
    mcs_result = rdFMCS.FindMCS(mol_list,
                                maximizeBonds=True,
                                completeRingsOnly=True,
                                ringMatchesRingOnly=True,
                                atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes,
                                )
    comm_num = mcs_result.numAtoms
    if comm_num < 3:
        return None
    else:
        tempmol = Chem.MolFromSmarts(mcs_result.smartsString)
        matches = mol_list[0].GetSubstructMatches(tempmol, uniquify=True,
                                                  useChirality=True,
                                                  useQueryQueryMatches=True,
                                                  maxMatches=1000)
        print('matches', matches)
        common_struct_smiles = Chem.MolFragmentToSmiles(mol_list[0], atomsToUse=matches[0], isomericSmiles=True,
                                                        canonical=False)
        print(common_struct_smiles)

        return common_struct_smiles


def compare_mcs(mols):
    groups = []
    other_groups = []
    stop_find = False
    for num in range(len(mols), 0, -1):
        print('the loop is:', num)
        if stop_find:
            break
        # sub_mols = mols[0:num]
        sub_mols_set = combinations(mols, num)
        for sub_mols in sub_mols_set:
            mcs_smi = GetMCSInfo(sub_mols)
            if mcs_smi == None:
                pass
            else:
                groups.append(sub_mols)
                if num < len(mols) - 1:
                    other_groups.append(list(set(mols) - set(sub_mols)))
                stop_find = True
                break
        if num == 1:
            other_groups.append(mols)
            print('Be careful, No mol has common structure!!!')
            break
    return groups, other_groups, mcs_smi


def CalMolFormalCharge(mol):
    mol_formal_charge = 0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        mol_formal_charge += charge
        # sym = atom.GetSymbol()
    return mol_formal_charge


def compare_charge(mols):
    valid_groups = []
    invalid_groups = []
    charges = [CalMolFormalCharge(mol) for mol in mols]
    maxTimes = max(charges, key=charges.count)
    for i, charge in enumerate(charges):
        if charge == maxTimes:
            valid_groups.append(mols[i])
        else:
            invalid_groups.append(mols[i])
    return valid_groups, invalid_groups, maxTimes
# contain .mols


def CheckMol(path, path_valid_mols, path_valid_info, path_invlid_info):
    """
    check mols validate; common structure; charge; check mol_name
    input: the path contains many "*.mol" filed
    """
    if not os.path.exists(path):
        print('No Valid Ligands Pathway was found ')
        exit()
    if os.path.exists(path_valid_mols):
        print('Output directory Pathway already existed ')
        exit()
    os.makedirs(path_valid_mols, exist_ok=True)

    if os.path.isfile(path_valid_info):
        os.remove(path_valid_info)
    if os.path.isfile(path_valid_info):
        os.remove(path_valid_info)

    valid_mols = {}
    invalid_mols = {}

    valid_Chemical_Rule = []
    mols = os.listdir(path)
    mols = [mol for mol in mols if ('.mol' in mol) or ('.sdf' in mol)]

    for mol_name in mols:   # check chemistry
        try:
            base_name = mol_name.split('.')[0]
            style = mol_name.split('.')[-1]
            mol_path = os.path.join(path, mol_name)
            if style == 'mol':
                # print('mol_path',mol_path)
                mol = Chem.MolFromMolFile(
                    mol_path, sanitize=True, removeHs=False)
            elif style == 'sdf':
                mol = Chem.SDMolSupplier(
                    mol_path, sanitize=True, removeHs=False)[0]
            else:
                print('We only support mol or sdf format file, Please check your file')
                # exit()
            mol.SetProp('_Name', base_name)
            valid_Chemical_Rule.append(mol)
            print('base-name', mol_name)
        except:
            print('Invalid Chemical Rules', mol_name)
            invalid_mols[base_name] = 'Invalid Chemical Rules'
    print('valid_Chemical_Rule:', valid_Chemical_Rule)
    charge_valid_groups, charge_invalid_groups, valid_charge = compare_charge(
        valid_Chemical_Rule)
    if len(charge_valid_groups) <= 0:
        print('The valid charge group is too small, Yhe number of valid charge group must large than two.')
        exit()
    if len(charge_invalid_groups):
        for in_m in charge_invalid_groups:
            name = in_m.GetProp('_Name')
            in_ch = CalMolFormalCharge(in_m)
            invalid_mols[name] = f'Charge Inconsistency, formal_charge={in_ch}'

    mcs_valid_groups, other_groups, mcs_smi = compare_mcs(charge_valid_groups)
    if len(other_groups) > 0:
        for in_m in other_groups[0]:
            name = in_m.GetProp('_Name')
            invalid_mols[name] = f'No Common Substructure'

    if len(mcs_valid_groups) > 0:
        valid_mols['msc'] = mcs_smi
        valid_mols['charge'] = valid_charge
        name_ob = {}
        for mol in mcs_valid_groups[0]:
            name = mol.GetProp('_Name')
            name_ob[name] = Chem.MolToSmiles(mol)
            mol_name = name+'.mol'
            path_valid_mol = os.path.join(path_valid_mols, mol_name)
            Chem.MolToMolFile(mol, path_valid_mol)
        valid_mols['mols'] = name_ob

        df_data = pd.DataFrame(valid_mols)
        df_data.to_csv(path_valid_info)
    else:
        print('Faild found purterbation map')
    if len(invalid_mols.keys()) > 0:
        df_data = pd.DataFrame(list(invalid_mols.items()))
        df_data.to_csv(path_invlid_info, index=False, header=False)


def load_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='the path of check_ligands json file')
    parser.add_argument('--file_path', '-f', help='where is your json file?',
                        default='check_ligs_input.json')
    args = parser.parse_args()

    try:
        data = load_json(args.file_path)
    except:
        print('cannot find your json file')
        exit()
    CheckMol(data['path'], data['path_valid_mols'],
             data['path_valid_info'], data['path_invalid_info'])
