import pandas as pd
import numpy as np

core_pdbs = pd.read_csv('/mnt/research/woldring_lab/Members/Eaves/plip-plop/index/PDBbind_v2020_core.csv')['pdb_id'].to_list()

df_list = []
with open('PDBbind_v2020_index_PL_ORIGINAL.txt', 'r') as raw_f:
    for line in raw_f:
        if line.startswith('#'):
            continue
        pdb_id, resolution, year, target = line.split('  ')[:4]
        
        if pdb_id not in core_pdbs:
            continue
        try:
            target_type, target_value = target.split('=')
        except:
            target_value = target
        target_num = float(target_value[:-2])
        target_unit = target_value[-2:]
        
        if target_unit == 'uM':
            target_num = target_num*1e-6
        elif target_unit == 'mM':
            target_num = target_num*1e-3
        elif target_unit == 'nM':
            target_num = target_num*1e-9
        elif target_unit == 'pM':
            target_num = target_num*1e-12
        elif target_unit == 'fM':
            target_num = target_num*1e-15

        pk = -np.log10(target_num)
        
        df_list.append({'dataset': 'core', 'pdb_id': pdb_id, 'resolution': resolution, 'release_year': year, 'Kd_M': target_num, 'pK': pk})
        
df = pd.DataFrame(df_list)

df.to_csv('/mnt/research/woldring_lab/Members/Eaves/plip-plop/index/PDBbind_v2020_core-TARGET-FIX.csv', index=False)