
import springcraft
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb
import numpy as np
import biotite_utils
import constants
import os

def generate_fluctuation(id, auth_binding_indices, output_path):
    protein = biotite_utils.load_structure(id)
    # some errors with non-standard residue types
    protein_backbone = biotite_utils.get_protein_backbone(protein, id[4:])
    ff = springcraft.TabulatedForceField.sd_enm(protein_backbone)
    eanm = springcraft.GNM(protein_backbone, ff)
    
    eig_values, eig_vectors = eanm.eigen()
    eig_vectors = np.square(eig_vectors)    
    fluctuation = np.zeros(eig_vectors.shape)

    for i in range(len(eig_vectors)):
        fluctuation[i] = eig_vectors[i] / eig_values[i]
    
    assert len(protein_backbone) == len(fluctuation), id

    print('\tSaving to file ...')

    binding_indices = [idx for idx in range(len(protein_backbone)) if protein_backbone[idx].res_id in auth_binding_indices]
    np.save(f'{output_path}/indices/{id}.npy', np.array(binding_indices))
    np.save(f'{output_path}/fluctuation/{id}.npy', np.array(fluctuation))
