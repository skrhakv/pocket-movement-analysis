
import springcraft
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb
import numpy as np
import biotite_utils
import os

def generate_fluctuation(id, auth_binding_indices, mmcif_file_location, output_path):
    mmcif_filename = f'{id[:4]}.cif'

    mmcif_filepath = f'{mmcif_file_location}/{mmcif_filename}'

    # download mmCIF file if doesn't exist
    if not os.path.isfile(mmcif_filepath):
        print(f'\tDownloading {mmcif_filepath} ...')
        rcsb.fetch(id[:4], 'cif', target_path=mmcif_file_location)

    mmcif_file = pdbx.CIFFile.read(mmcif_filepath)
    # load file to biotite object
    whole_structure = pdbx.get_structure(mmcif_file, model=1, include_bonds=True)
    protein = whole_structure[struc.filter_amino_acids(whole_structure)]

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
