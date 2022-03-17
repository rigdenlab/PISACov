## Full path to pdb_chain_uniprot.csv SIFTS database.
SIFTS_PATH = '/media/javier/databases/pdb_chain_uniprot.csv'

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

## Path to DeepMetaPsicov script
DMP_PATH='/home/javier/Programs/DeepMetaPSICOV/run_DMP.sh'

## Input parameters for HHBLITS search. Uncomment if not using DeepMetaPSICOV default.

## (#iterations, E-value cutoff, Non-redundant seqs to keep, MinimumCoverageWithMasterSeq(%),MaxPairwiseSequenceIdentity)
## DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99]
## HHBLITS DEFAULT: [2, 0.001, 1000, 0, 90]
## Uncomment next line to use non-DeepMetaPSICOV-default values
#HHBLITS_PARAMETERS=[3, 0.001, 'inf', 50, 99]

## CONTACT PARAMETERS

## By default, pairs of residues that appear as intermolecular contacts
## and intramolecular contacts too are removed with NEIGHBOURS_MINDISTANCE set to 2;
## Intramolecular contacts determine this threshold

## Minimum distance within the sequence for neighbours to be accounted for (int)
## |pos_res1-pos_res2|>=MinDistance; Minimum value: 2; Default: 2
## Uncomment next line if you want to bypass the default behaviour and fix a custom minimum distance while ignoring intramolecular contacts.
#NEIGHBOURS_MINDISTANCE = 2
