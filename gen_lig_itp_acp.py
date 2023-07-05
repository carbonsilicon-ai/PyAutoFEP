import os
import argparse


def GenLigItp(atom_type, lig_dir, output_dir):
    mols = os.listdir(lig_dir)
    mols = [mol for mol in mols if '.mol' in mol]
    for mol in mols:
        mol_path = os.path.join(lig_dir, mol)
        cmd = f'acpype -i {mol_path} -a {atom_type}'
        mol_name = mol.split('.')[0]
        itp = f'{mol_name}.acpype/{mol_name}_GMX.itp'
        itp_path = os.path.join(output_dir, itp)
        new_itp_path = os.path.join(lig_dir, f'{mol_name}.itp')
        try:
            os.system(cmd)
        except:
            print(f'{mol} generated itp failed')
        cmd1 = f'cp {itp_path} {new_itp_path}'
        os.system(cmd1)
        acpype_dir = os.path.join(output_dir, f'{mol_name}.acpype')
        os.system(f'rm -rf {acpype_dir}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='the path of check_ligands json file')
    parser.add_argument('--lig_dir', '-input', help='where is your ligands file?',
                        default='aligned_lig_data')
    parser.add_argument('--atom_type', '-atom', help='which atom type do you choose?',
                        default='gaff2')
    parser.add_argument('--output_dir', '-output', help='where is your output dir?',
                        default='lig_itp_data')
    args = parser.parse_args()

    if not os.path.exists(args.lig_dir):
        print('No Valid Ligands Pathway was found ')
        exit()
    if os.path.exists(args.output_dir):
        print('Output directory Pathway already existed ')
        exit()
    os.makedirs(args.output_dir, exist_ok=True)

    GenLigItp(args.atom_type, args.lig_dir, args.output_dir)
