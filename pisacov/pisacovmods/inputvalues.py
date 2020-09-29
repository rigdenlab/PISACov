#####################################################################
##### INPUT VALUES ##################################################

########## OUTPUT
# Full path to directory where output files will be created. Note that a /PDB_ID/ directory will be created for each protein
OUTPUT_DIR='/media/javier/moredata/javier/PD_RigdenLab/Task4/pythonCodes/python2/py3_output/'

########## Protein files
# Full path to PDB file
PDB_PATH='/media/javier/moredata/javier/PD_RigdenLab/Protein_Files/2POD/SourceFiles/2pod.pdb'

# Full path to protein sequence file (fasta)
SEQUENCE_PATH='/media/javier/moredata/javier/PD_RigdenLab/Protein_Files/2POD/SourceFiles/2pod.fasta'

# Use biological chain (True) or whole PDB chain, including artifacts (False) to produce MSA and subsequent covariance analysis?
USE_BIOCHAIN=True

# Full path to SIFTS libraries (CSV) containing information about biological N-terminal and C-terminal residues within the pdb protein sequence
CSV_CHAIN_PATH='/media/javier/databases/pdb_chain_uniprot.csv'
CSV_SEGMENTS_PATH='/media/javier/databases/uniprot_segments_observed.csv'

########## PISA files
# Full path to PISA script
PISA_PATH='/home/javier/Programs/ccp4-7.0/bin/pisa'

########## HH-suite files
# Full path to hh-suite script
HHSUITE_PATH='/home/javier/Programs/hh-suite/build/bin/hhblits'
# Name of sequence database to be used by HHBLITS
HHBLITS_DATABASE_NAME='uniclust30_2018_08'  # e.g. uniclust30_2018_08
# Full path to directory containing sequence database to be used by HHBLITS
HHBLITS_DATABASE_DIR='/media/javier/databases/uniclust30_2018_08_hhblits/'
# Input parameters for HHBLITS search: ( #iterations (def:2), E-value cutoff (def:0.001), Non-redundant seqs to keep (def=1000), MinimumCoverageWithMasterSeq(%)(def=0),MaxPairwiseSequenceIdentity (def=90))
HHBLITS_PARAMETERS=[ 2 , 0.001, 1000, 0, 90 ]
# Use DeepMetaPSICOV execution of HHBLITS instead ( #iterations=3, E-value cutoff=0.001, Non-redundant seqs to keep=inf,MinimumCoverageWithMasterSeq(%)=50,MaxPairwiseSequenceIdentity=99   )
HHBLITS_VIA_DMP=True

########## DMP files
# Path to DeepMetaPsicov script
DMP_PATH='/home/javier/Programs/DeepMetaPSICOV/run_DMP.sh'
########## CONTACT PARAMETERS
# Minimum distance within the sequence for neighbours to be accounted for (|pos_res1-pos_res2|>=MinDistance, recommended: 5, minimum value: 2, do not remove neighbours: 0 or 1)
NEIGHBOURS_MINDISTANCE=2
# Score threshold below which predicted contacts (DMP) will be ignored (def: 0.2, do not remove contacts: 0.0)
DMPSCORE_THRESHOLD=0.0
# Score threshold below which predicted contacts (PSICOV) will be ignored (def: 0.2 *doi: 10.1038/s41598-019-48913-8, do not remove contacts: 0.0)
PSICOVSCORE_THRESHOLD=0.2
# Remove intermolecular contacts that are also intramolecular contacts too (bool) (If True, NEIGHBOURS_MINDISTANCE will be set to 2; Intramolecular contacts determine this threshold)
REMOVE_INTRA_CONTACTS=True

########## LOG-STDOUT
# If stdout is to be stored in a log-file (True) or shown on screen (False)
LOG_STDOUT=True