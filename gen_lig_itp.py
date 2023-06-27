

import os
import subprocess

def make_lig_top(lig_path,sobtop_path):
    # ligand.mol2 file
    #lig_path = '/home/chenchengzhao/Auto_MMPBSA/ligands/CHEMBL3126182_0/ligand.mol2'
    lig_m2_path = lig_path.replace('.mol', '.mol2')
    os.system(f'obabel  {lig_path} -O {lig_m2_path}')
    lig_gro_path = lig_path.replace('.mol', '.gro')
    lig_top_path = lig_path.replace('.mol', '.top')
    lig_itp_path = lig_path.replace('.mol', '.itp')
    bash_path = lig_path.replace('.mol', '_make_top.sh')
    #os.chdir(sobtop_path)
    l0 = '#!/bin/bash\n\n\n'
    #ll = f'echo -e "{lig_m2_path}\\n2\\n{lig_gro_path}\\n1\\n2\\n4\\n{lig_top_path}\\n{lig_itp_path}"| ./sobtop'
    ll1 = f'{{(cd {sobtop_path} && echo -e "{lig_m2_path}\\n2\\n{lig_gro_path}\\n1\\n2\\n4\\n{lig_top_path}\\n{lig_itp_path}" | ./sobtop)}}\n\n'
    #os.system(ll)
    with open(bash_path,'w') as f:
        f.writelines([l0,ll1])
    os.system(f'bash {bash_path}')
    print(ll1)
    os.system(f'rm {bash_path}  {lig_top_path}  {lig_gro_path}')

if __name__ == '__main__':
    lig_ra_path = 'lig_data/FXR_12.mol'
    lig_path = os.path.join(os.getcwd(),f'{lig_ra_path}')
    print('lig_path',lig_path)
    sobtop_path = '/home/work/sobtop_1.0'

    make_lig_top(lig_path,sobtop_path)



