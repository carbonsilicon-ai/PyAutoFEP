
import os
import argparse

#atom_type = 'gaff'
#lig_dir = 'lig_data'

def GenLigItp(lig_dir,atom_type):
    mols = os.listdir(lig_dir)
    mols = [mol for mol in mols if '.mol' in mol]
    curr_dir = os.getcwd()
    for mol in mols:
        mol_path = os.path.join(lig_dir, mol)
        cmd = f'acpype -i {mol_path} -a {atom_type}'
        mol_name = mol.split('.')[0]
        itp = f'{mol_name}.acpype/{mol_name}_GMX.itp'
        itp_path = os.path.join(curr_dir, itp)
        new_itp_path = os.path.join(lig_dir, f'{mol_name}.itp')
        try:
            os.system(cmd)
        except:
            print(f'{mol} generated itp failed')
        cmd1 = f'cp {itp_path} {new_itp_path}'
        os.system(cmd1)
        acpype_dir = os.path.join(curr_dir, f'{mol_name}.acpype')
        os.system(f'rm -rf {acpype_dir}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='the path of check_ligands json file')
    parser.add_argument('--lig_dir', '-input', help='where is your ligands file?',
                        default='/home/work/pyautofep/PyAutoFEP/aligned_lig_data')
    parser.add_argument('--atom_type', '-atom', help='which atom type do you choose?',
                        default='gaff2')
    args = parser.parse_args()
    lig_dir = args.lig_dir
    atom_type = args.atom_type
    GenLigItp(lig_dir, atom_type)




