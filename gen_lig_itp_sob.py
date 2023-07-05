

import os
import subprocess
import json
import argparse
def make_lig_top(lig_path,sobtop_path):
    # ligand.mol2 file
    #lig_path = '/home/chenchengzhao/Auto_MMPBSA/ligands/CHEMBL3126182_0/ligand.mol2'
    lig_m2_path = lig_path.replace('.mol', '.mol2')
    os.system(f'obabel {lig_path} -O {lig_m2_path}')
    lig_gro_path = lig_path.replace('.mol', '.gro')
    lig_top_path = lig_path.replace('.mol', '.top')
    lig_itp_path = lig_path.replace('.mol', '.itp')
    #lig_atp_path = lig_path.replace('.mol2', '.atp')
    #sobtop_path = '/home/chenchengzhao/sobtop_1.0'
    bash_path = lig_path.replace('.mol', '_make_top.sh')
    #os.chdir(sobtop_path)
    l0 = '#!/bin/bash\n\n\n'
    #ll = f'echo -e "{lig_path}\\n2\\n{lig_gro_path}\\n1\\n2\\n4\\n{lig_top_path}\\n{lig_itp_path}"| ./sobtop'
    ll1 = f'{{(cd {sobtop_path} && echo -e "{lig_m2_path}\\n2\\n{lig_gro_path}\\n1\\n2\\n4\\n{lig_top_path}\\n{lig_itp_path}" | ./sobtop)}}\n\n'
    #os.system(ll)
    with open(bash_path,'w') as f:
        f.writelines([l0,ll1])
    os.system(f'bash {bash_path}')
    os.system(f'rm {lig_m2_path} {bash_path} {lig_top_path} {lig_gro_path}')

def GenLigsItp(lig_dir,sobtop_path):
    ligs = os.listdir(lig_dir)
    ligs = [lig for lig in ligs if '.mol' in lig ]
    for lig in ligs:
        lig_path = os.path.join(lig_dir,lig)
        make_lig_top(lig_path, sobtop_path)

if __name__ == '__main__':
    #lig_ra_path = 'lig_data/FXR_12.mol'
    #lig_path = os.path.join(os.getcwd(),f'{lig_ra_path}')
    parser = argparse.ArgumentParser(description='gen itp file')
    parser.add_argument('--lig_dir', '-input', help='where is your ligands file?',
                        default='aligned_lig_data')
    parser.add_argument('--sobtop_path', '-sob', help='which atom type do you choose?',
                        default='/home/work/sobtop_1.0')
    args = parser.parse_args()
    lig_dir = args.lig_dir
    sobtop_path = args.sobtop_path
    GenLigsItp(lig_dir, sobtop_path)

