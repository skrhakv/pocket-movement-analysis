import json


def get_annotations(dataset_path, whole_dataset=False):
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
                binding_residues = set([int(''.join(filter(str.isdigit, residue.split('_')[1]))) for residue in holo_structure['apo_pocket_selection']])
                if id in ids: ids[id].update(binding_residues)
                else: ids[id] = binding_residues
        return ids

    ids = {}
    if not whole_dataset:
        for fold in [f'train-fold-{i}.json' for i in range(4)]:
            subset_ids = load_apo_structures(f'{dataset_path}/folds/{fold}')
            ids = {**ids, **subset_ids}
    else:
        ids = load_apo_structures(f'{dataset_path}/dataset.json')
    return ids 


def load_subset(path):
    with open(path) as f:
        dataset = json.load(f)
    return dataset

def load_train_set(dataset_path):
    train_set = {}
    for fold in [f'train-fold-{i}.json' for i in range(4)]:
        train_subset = load_subset(f'{dataset_path}/folds/{fold}')
        train_set = {**train_set, **train_subset}
    return train_set

def load_main_apo_holo_pairs(dataset, multichain=True):
    """Retrieves 'main' HOLO structure for each APO structure

    Args:
        dataset (JSON): dataset
    """
    apo_to_holo = {}
    for apo_pdb_id, holo_structures in dataset.items():
        
        # rigid dataset versions don't have the 'is_main_holo_structure' - we need to sumpplement that
        if 'is_main_holo_structure' in holo_structures[0]:
            has_main_holo_structure = False
            for holo_structure in holo_structures:
                holo_pdb_id = holo_structure['holo_pdb_id']
                if holo_structure['is_main_holo_structure']:
                    apo_chain_id = holo_structure['apo_chain']
                    holo_chain_id = holo_structure['holo_chain']
                    if not multichain and ('-' in apo_chain_id or '-' in holo_chain_id):
                        has_main_holo_structure = True
                        break
                    apo_to_holo[f'{apo_pdb_id}{apo_chain_id}'] = f"{holo_pdb_id}{holo_chain_id}"
                    has_main_holo_structure = True
                    break
            assert has_main_holo_structure, f"No main holo structure found for {apo_pdb_id}"
        
        else:
            min_pRMSD = float('inf')
            min_pRMSD_holo_pdb_id = None
            min_pRMSD_apo_chain_id = None 
            min_pRMSD_holo_chain_id = None

            for holo_structure in holo_structures:
                apo_chain_id = holo_structure['apo_chain']
                holo_chain_id = holo_structure['holo_chain']

                if not multichain and ('-' in apo_chain_id or '-' in holo_chain_id):
                    continue
                if holo_structure['pRMSD'] < min_pRMSD:
                    min_pRMSD = holo_structure['pRMSD']
                    min_pRMSD_apo_chain_id = apo_chain_id
                    min_pRMSD_holo_chain_id = holo_chain_id
                    min_pRMSD_holo_pdb_id = holo_structure['holo_pdb_id']
            
            # can happen when all holos are multichain
            if min_pRMSD_holo_pdb_id is None:
                continue

            apo_to_holo[f'{apo_pdb_id}{min_pRMSD_apo_chain_id}'] = f"{min_pRMSD_holo_pdb_id}{min_pRMSD_holo_chain_id}"
    return apo_to_holo
