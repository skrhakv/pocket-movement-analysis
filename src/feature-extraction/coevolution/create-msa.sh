#!/usr/bin/env bash

rm data/tmp.sto data/tmp.hmm data/tmp.out.sto data/tmp.output.trimmed.fasta data/tmp.retrieved-sequences.sto data/tmp.retrieved-sequences.fasta data/msa.sto data/msa.fasta
# Obtain the profiles using the Uniprot proteome reference
phmmer -A data/tmp.sto data/tmp.fasta data/prot-ref.fasta

# Build the profiles using the `hmmbuild` command
hmmbuild --fragthresh 1 data/tmp.hmm data/tmp.sto

# Use the profile to search against the PDB sequences
hmmsearch -A data/tmp.retrieved-sequences.sto data/tmp.hmm data/pdb_seqres.txt

# Put the original sequence into the retrieved PDB sequences and convert to FASTA
esl-reformat -u fasta data/tmp.retrieved-sequences.sto > data/tmp.retrieved-sequences.fasta
cat data/tmp.fasta >> data/tmp.retrieved-sequences.fasta

# Align the original sequence with the retrieved PDB sequences using the HMM profile and convert to FASTA
hmmalign -o data/msa.sto data/tmp.hmm data/tmp.retrieved-sequences.fasta
