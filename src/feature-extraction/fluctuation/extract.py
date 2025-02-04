import json
import sys
import os
import shutil
sys.path.append('../../utils/')
import springcraft_utils

DATASET = 'rigid-dataset'
DATASET_PATH = f'../../../datasets/{DATASET}'
OUTPUT_PATH = f'../../../data/features/fluctuation/{DATASET}'
OLD_OUTPUT_PATH = F'../../../data/features/fluctuation/old/rigid-dataset'

def load_apo_structures(path):
    with open(path) as f:
        dataset = json.load(f)
    
    ids = {}
    for apo_id, holo_structures in dataset.items():
        for holo_structure in holo_structures:
            # skip multichain structures
            if '-' in holo_structure['apo_chain']:
                continue
            id = apo_id + holo_structure['apo_chain']
            if id in ids: ids[id].update(holo_structure['apo_pocket_selection'])
            else: ids[id] = set([int(''.join(filter(str.isdigit, residue.split('_')[1]))) for residue in holo_structure['apo_pocket_selection']])
    return ids

ids = {}
for fold in [f'train-fold-{i}.json' for i in range(4)]:
    subset_ids = load_apo_structures(f'{DATASET_PATH}/folds/{fold}')
    ids = {**ids, **subset_ids}

for id in ids.keys():
    if f'{id}.npy' in os.listdir(f'{OUTPUT_PATH}/fluctuation'):
        continue
    if f'{id}.npy' in os.listdir(f'{OLD_OUTPUT_PATH}/fluctuation'):
        shutil.copy(f'{OLD_OUTPUT_PATH}/fluctuation/{id}.npy', f'{OUTPUT_PATH}/fluctuation')   
        shutil.copy(f'{OLD_OUTPUT_PATH}/indices/{id}.npy', f'{OUTPUT_PATH}/indices')   
    print(f'Processing {id} ...')
    # read file
    mmcif_filename = f'{id[:4]}.cif'
    springcraft_utils.generate_fluctuation(id, ids[id], OUTPUT_PATH)
