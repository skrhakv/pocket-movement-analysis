import constants
import os
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb
from biotite.sequence import ProteinSequence
import biotite.sequence.align as align

def get_protein_backbone(protein, chain_id, indices=None):
    """Get the protein backbone atoms from a protein structure.
    
    Parameters
    ----------
    protein : The protein structure.
    chain_id : The chain ID of the protein.
    indices : The indices of the residues to return.
    Returns
    -------
    protein_backbone : AtomArray
        The protein backbone atoms.
    """
    protein_backbone = protein[(protein.atom_name == "CA") 
                       & (protein.element == "C") 
                       & (protein.chain_id == chain_id) 
                       & (
                             (protein.res_name == 'ALA')
                           | (protein.res_name == 'ARG')
                           | (protein.res_name == 'ASN')
                           | (protein.res_name == 'ASP')
                           | (protein.res_name == 'CYS')
                           | (protein.res_name == 'GLN')
                           | (protein.res_name == 'GLU')
                           | (protein.res_name == 'GLY')
                           | (protein.res_name == 'HIS')
                           | (protein.res_name == 'ILE')
                           | (protein.res_name == 'LEU')
                           | (protein.res_name == 'LYS')
                           | (protein.res_name == 'MET')
                           | (protein.res_name == 'PHE')
                           | (protein.res_name == 'PRO')
                           | (protein.res_name == 'SER')
                           | (protein.res_name == 'THR')
                           | (protein.res_name == 'TRP')
                           | (protein.res_name == 'TYR')
                           | (protein.res_name == 'VAL'))]
    
    if indices:
        return protein_backbone[indices]
    return protein_backbone

def get_sequence(protein, chain_id, from_backbone=False):
    """Retrieves the sequence of a protein chain.

    Args:
        protein : PDB ID of the protein
        chain_id : Chain ID of the protein
        from_backbone : If True, the step to retrieve backbone atoms is skipped
    """
    if not from_backbone:
        backbone = get_protein_backbone(protein, chain_id)
    else:
        backbone = protein
    sequence = ''.join([ProteinSequence.convert_letter_3to1(residue.res_name) for residue in backbone])
    return ProteinSequence(sequence)

def align_sequences(sequence1, sequence2):
    """Aligns two sequences.

    Args:
        sequence1 (ProteinSequence): sequence 1
        sequence2 (ProteinSequence): sequence 2
    """
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    alignments = align.align_optimal(
        sequence1, sequence2, matrix)
    return alignments

def load_structure(id, use_author_fields=True):
    """Downloads structure from PDB and returns the protein part of the structure.

    Args:
        id : ID of the protein (e.g. 1a2bA)
    """
    mmcif_filename = f'{id[:4]}.cif'

    mmcif_filepath = f'{constants.CIF_FILES_PATH}/{mmcif_filename}'

    # download mmCIF file if doesn't exist
    if not os.path.isfile(mmcif_filepath):
        print(f'\tDownloading {mmcif_filepath} ...')
        rcsb.fetch(id[:4], 'cif', target_path=constants.CIF_FILES_PATH)

    mmcif_file = pdbx.CIFFile.read(mmcif_filepath)
    # load file to biotite object
    whole_structure = pdbx.get_structure(mmcif_file, model=1, include_bonds=True, use_author_fields=use_author_fields)
    return whole_structure[struc.filter_amino_acids(whole_structure)]

    