import sys
sys.path.append('../../utils/')
import dataset_utils

import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc
import numpy as np
import os
import subprocess
from biotite.sequence import ProteinSequence
from pydca.msa_trimmer import msa_trimmer
from pydca.plmdca import plmdca

from Bio import AlignIO
import numpy as np

CIF_FILES_PATH = '/home/vit/Projects/deeplife-project/data/cif_files'
SKIPPED_IDS = ['5x6zD', '4ok2B', '5h61B', '1fl1B', '8hynA', '1sh0A', '2yx7A', '6lgyA', '4gu8H', '8j1kA', '1xhxB', '1pfzC', '4omgA', '1g24C', '6j35B', '6ialF', '4j2pA', '4x1oA', '5ytbC',
               '1havA', '4ekfA', '8iy0B', '3l28F', '2xc1A', '7wyoB', '3cbjA', '1of3B', '7gosB', '7poqA', '6spoA', '2wetA', '6tyoB', '7d48A', '2zf3F', '1mwkA', '5mf2D', '7oapEEE', '8dufA',
               '4pvrA', '4gf1B', '4u0mB', '1xt3B', '6ro0G', '3sebA', '3rmyD', '6twoAAA', '1nokA', '7oumC', '3ugjA', '2yglB', '4gzqA', '6d2cA', '6niwD']


def compute_coevolution(ids, dataset):
    fluctuation_path = f'/home/vit/Projects/flexibility-analysis/data/features/fluctuation/{dataset}/fluctuation'
    output_path = f'/home/vit/Projects/flexibility-analysis/data/features/coevolution/{dataset}'

    def get_coevolution_len(coevolution_pairs):
        min_node = float('inf')
        max_node = float('-inf')
        for i in coevolution_pairs:
            node1 = i[0][0]
            node2 = i[0][1]

            min_node = min(min_node, node1)
            min_node = min(min_node, node2)
            max_node = max(max_node, node1)
            max_node = max(max_node, node2)

        return len(range(min_node, max_node + 1))

    for id in ids.keys():
        if id in SKIPPED_IDS: continue
        if f'{id}.npy' in os.listdir(output_path):
            continue
        print(f'Processing {id} ...')
        # read file
        mmcif_filename = f'{id[:4]}.cif'
        mmcif_file = pdbx.CIFFile.read(f'{CIF_FILES_PATH}/{mmcif_filename}')

        # load file to biotite object
        whole_structure = pdbx.get_structure(mmcif_file, model=1, include_bonds=True)
        protein = whole_structure[struc.filter_amino_acids(whole_structure)]

        # some errors with MSE residue
        c_alphas = protein[(protein.atom_name == "CA") 
                           & (protein.element == "C") 
                           & (protein.chain_id == id[4:]) 
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

        print('len: ', len(c_alphas))

        # better check with the fluctuation data too:
        fluctuation_data_len = len(np.load(f'{fluctuation_path}/{id}.npy'))

        # sometimes the fluctuation length doesn't match
        # maybe because of different python version? who knows ..
        # anyway, skip such structures:
        if len(c_alphas) != fluctuation_data_len:
            print(f'Skipping {id}: Fluctuation data length doesn\'t match the sequence length! {len(c_alphas)} {fluctuation_data_len}')
            continue

        # get sequence
        sequence = ''.join([ProteinSequence.convert_letter_3to1(residue.res_name) for residue in c_alphas])

        # save it to tmp.fasta for the hmmer tool
        with open('data/tmp.fasta', 'w') as f:
            f.write(f'>{id}\n')
            f.write(sequence)

        # call hmmer to compute MSA
        subprocess.call(['sh', './create-msa.sh'])
        msa_file = 'data/msa.sto'

        if not os.path.isfile(msa_file) or os.stat(msa_file).st_size == 0:
            print(f'Skipping {id} - unsuccessful MSA calculation ...')
            continue
        
        # convert MSA from Stockholm to FASTA format
        AlignIO.convert('data/msa.sto', 'stockholm', 'data/msa.fasta', 'fasta')

        print('\tCalculating coevolution ...')
        refseq_file = 'data/tmp.fasta'
        msa_file = 'data/msa.fasta'

        trimmer = msa_trimmer.MSATrimmer(
            msa_file, biomolecule='protein', 
            refseq_file=refseq_file,
            max_gap=0.0
        )

        trimmed_data = trimmer.get_msa_trimmed_by_refseq()

        trimmed_data_outfile = 'data/tmp.output.trimmed.fasta'
        with open(trimmed_data_outfile, 'w') as fh:
            fh.write('>{}\n{}\n'.format(id, sequence))
            for seqid, seq in trimmed_data:
                fh.write('>{}\n{}\n'.format(seqid, seq))

        # compute coevolution
        coevolving_pairs_scores = plmdca.PlmDCA(
            trimmed_data_outfile,
            'protein',
            seqid = 0.8,
            lambda_h = 1.0,
            lambda_J = 20.0,
            num_threads = 10,
            max_iterations = 500,
        )
        # compute DCA scores summarized by Frobenius norm and average product corrected
        sorted_coevolving_pairs_scores = coevolving_pairs_scores.compute_sorted_FN_APC()

        coevolution_len = get_coevolution_len(sorted_coevolving_pairs_scores)

        assert len(c_alphas) == coevolution_len, f'{id}: coevolution data length doesn\'t match the sequence length! {len(c_alphas)} {coevolution_len}'

        a = np.array(sorted_coevolving_pairs_scores)

        np.save(f'{output_path}/{id}.npy', np.array(sorted_coevolving_pairs_scores))

if __name__ == "__main__":
    DATASET = 'rigid-dataset'
    DATASET_PATH = f'../../../datasets/{DATASET}'
    rigid_ids = dataset_utils.get_annotations(DATASET_PATH)
    compute_coevolution(rigid_ids, DATASET)