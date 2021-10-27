"""
PISACov is a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

## Full path to pdb_chain_uniprot.csv SIFTS database.
CSV_CHAIN_PATH = '/media/javier/databases/pdb_chain_uniprot.csv'

## Uncomment if Uniprot fasta database is to be read locally.
## Full path to Uniprot fasta database.
#UNICLUST_FASTA_PATH = ''

## Full path to PISA script
PISA_PATH = '/home/javier/Programs/ccp4-7.0/bin/pisa'

## NOTE: HHBLITS paths below should be the same provided in DMP_run.sh
## Full path to hhblits script
HHBLITS_PATH = '/home/javier/Programs/hh-suite/build/bin/hhblits'
## Name of sequence database to be used by HHBLITS (as requested by DeepMetaPSICOV)
HHBLITS_DATABASE_NAME = 'uniclust30_2018_08'  # e.g. uniclust30_2018_08
## Full path to directory containing sequence database to be used by HHBLITS
HHBLITS_DATABASE_DIR = '/media/javier/databases/uniclust30_2018_08_hhblits/'

## Input parameters for HHBLITS search. Uncomment if not using DeepMetaPSICOV default.
## (#iterations, E-value cutoff, Non-redundant seqs to keep, MinimumCoverageWithMasterSeq(%),MaxPairwiseSequenceIdentity)
## DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99]
## HHBLITS DEFAULT: [2, 0.001, 1000, 0, 90]
## Uncomment to use non-DeepMetaPSICOV-default values
#HHBLITS_PARAMETERS=[3, 0.001, 'inf', 50, 99]

## Path to DeepMetaPsicov script
DMP_PATH='/home/javier/Programs/DeepMetaPSICOV/run_DMP.sh'

## CONTACT PARAMETERS

## Remove intermolecular contacts that are intramolecular contacts too (bool)
## If True, NEIGHBOURS_MINDISTANCE will be set to 2; Intramolecular contacts determine this threshold)
REMOVE_INTRA_CONTACTS = True

## Minimum distance within the sequence for neighbours to be accounted for (int)
## |pos_res1-pos_res2|>=MinDistance; Minimum value: 2; Default: 2
## Uncomment if REMOVE_INTRA_CONTACTS==False;
#NEIGHBOURS_MINDISTANCE = 2