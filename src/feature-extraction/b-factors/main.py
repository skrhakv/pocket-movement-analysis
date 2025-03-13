from bfactor_utils import Chain3dStructure

import numpy as np
import os
import sys
sys.path.append('../../utils/')
import dataset_utils

DATASET = 'rigid-dataset'
DATASET_PATH = f'../../../datasets/{DATASET}'
OUTPUT_FOLDER = f'../../../data/features/B-factors/{DATASET}'
INPUT_FOLDER = '/home/vit/Projects/deeplife-project/data/cif_files'
PROTRUSION_RADIUS = 10


def get_annotations(structure, dataset, fold):
    def get_pdb_annotation(dataset, pdb_id):
        annotation = set()

        for holo_structure in dataset[pdb_id]:
            annotation.update(holo_structure['apo_pocket_selection'])
        return {i.split('_')[1] for i in annotation}

    auth_binding_residue_ids = get_pdb_annotation(
        dataset, structure.protein_id)

    binding_residue_indices = []
    for idx, potential_binding_residue in enumerate(structure.get_residue_list()):
        auth_binding_residue_id = potential_binding_residue.get_full_id()[3][1]
        if str(auth_binding_residue_id) in auth_binding_residue_ids:
            binding_residue_indices.append(idx)

    # with open(f'{OUTPUT_FOLDER}/annotation/{fold}.csv', 'a') as f:
    #     f.write(
    #         f'{structure.protein_id};{structure.chain_id};UNKNOWN;{" ".join([f"{structure.chain_id}_{i}" for i in sorted(binding_residue_indices)])};UNKNOWN\n')


def get_features(structure):
    # protrusion = structure.get_protrusion_vector(radius=PROTRUSION_RADIUS, protrusion_algorithm='atom_count_from_furthermost_atom_coords')
    # sequence = structure.compute_sequence()
    # sasa = structure.get_SASA_vector()
    bfactors = structure.get_bfactors(bfactor_algorithm='average_bfactor')

    # assert len(sasa) == len(sequence)
    # assert len(protrusion) == len(sequence)

    # save_features(f'{structure.protein_id}{structure.chain_id}',
    #               sequence, protrusion, sasa)
    save_features(f'{structure.protein_id}{structure.chain_id}',
                  None, None, None, bfactors)
    
def save_features(id, sequence, protrusion, sasa, bfactor):
    #with open(f'{OUTPUT_FOLDER}/sequence/{id}.txt', 'w') as f:
    #    f.write(sequence)
    # np.save(f'{OUTPUT_FOLDER}/protrusions/radius-10-furthermost-atom/{id}.npy', np.array(protrusion))
    #np.save(f'{OUTPUT_FOLDER}/sasa/{id}.npy', np.array(sasa))
    np.save(f'{OUTPUT_FOLDER}/{id}.npy', np.array(bfactor))


def get_fold(pdb_id, folds):
    for key, pdb_ids in folds.items():
        if pdb_id in pdb_ids:
            return key
    return -1
    # assert False, f"ERROR: {pdb_id} was not found in folds.json"


def main():
    apo_ids = [i.split('.')[0] for i in os.listdir('/home/vit/Projects/flexibility-analysis/data/features/fluctuation/cryptobench-dataset/fluctuation')]

    # with open(f'{CRYPTOBENCH_FOLDER}/folds.json') as f:
    #     folds = json.load(f)
# 
    # with open(f'{CRYPTOBENCH_FOLDER}/dataset.json') as f:
    #     cryptobench = json.load(f)

    for apo_id in apo_ids:
        print(f'Processing {apo_id} ...')

        if os.path.isfile(f'{OUTPUT_FOLDER}/{apo_id}.npy'):
            continue

        pdb_id = apo_id[:4]
        chain_id = apo_id[4:]
        structure = Chain3dStructure(
            pdb_id, chain_id, f'{INPUT_FOLDER}/{pdb_id}.cif', load=True)
        # fold = get_fold(pdb_id, folds)
# 
        # if fold == -1:
        #     continue

        get_features(structure)
        # get_annotations(structure, cryptobench, fold)
 

if __name__ == '__main__':
    main()
