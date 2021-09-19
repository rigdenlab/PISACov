#####################################################################
##### INPUT VALUES ##################################################

########## OUTPUT
# Full path to directory where output files will be created. Note that a /PDB_ID/ directory will be created for each protein
#OUTPUT_DIR='/media/javier/moredata/javier/PD_RigdenLab/Task4/pythonCodes/python2/py3_output/'
# ADD TO MAIN
########## Protein files
# Full path to PDB file
#PDB_PATH='/media/javier/moredata/javier/PD_RigdenLab/Protein_Files/2POD/SourceFiles/2pod.pdb'
# ADD TO MAIN
# Full path to protein sequence file (fasta)
#SEQUENCE_PATH='/media/javier/moredata/javier/PD_RigdenLab/Protein_Files/2POD/SourceFiles/2pod.fasta'
# ADD TO MAIN
# Use biological chain (True) or whole PDB chain, including artifacts (False) to produce MSA and subsequent covariance analysis?
#USE_BIOCHAIN=True
# ADD TO MAIN
# Full path to SIFTS libraries (CSV) containing information about biological N-terminal and C-terminal residues within the pdb protein sequence
CSV_CHAIN_PATH='/media/javier/databases/pdb_chain_uniprot.csv'
# KEEP
#CSV_SEGMENTS_PATH='/media/javier/databases/uniprot_segments_observed.csv'
# ADD TO MAIN? USEFUL? --> FASTA UNICLUST BETTER FOR UNIPROT LENGTHS
UNICLUST_FASTA_PATH=''
# KEEP


########## PISA files
# Full path to PISA script
PISA_PATH='/home/javier/Programs/ccp4-7.0/bin/pisa'
# KEEP

########## HH-suite files
# Full path to hh-suite script
HHSUITE_PATH='/home/javier/Programs/hh-suite/build/bin/hhblits'
# KEEP
# Name of sequence database to be used by HHBLITS
HHBLITS_DATABASE_NAME='uniclust30_2018_08'  # e.g. uniclust30_2018_08
# KEEP
# Full path to directory containing sequence database to be used by HHBLITS
HHBLITS_DATABASE_DIR='/media/javier/databases/uniclust30_2018_08_hhblits/'
# KEEP
# Input parameters for HHBLITS search:
# (#iterations, E-value cutoff, Non-redundant seqs to keep, MinimumCoverageWithMasterSeq(%),MaxPairwiseSequenceIdentity)
# DEEPMETAPSICOV & PISACOV DEFAULT: [3,0.001,'inf',50,99] - HHBLITS DEFAULT: [2,0.001,1000,0,90]
# HHBLITS_PARAMETERS=[ 2 , 0.001, 'inf', 50, 99 ] # UNCOMMENT FOR NON-DEFAULT VALUES
# KEEP
# HHBLITS_VIA_DMP=True
# REMOVE -> DEFAULT to DMP, NON-DEFAULT TO HHBLITS

########## DMP files
# Path to DeepMetaPsicov script
DMP_PATH='/home/javier/Programs/DeepMetaPSICOV/run_DMP.sh'
# KEEP
########## CONTACT PARAMETERS
# Minimum distance within the sequence for neighbours to be accounted for (|pos_res1-pos_res2|>=MinDistance, recommended: 5, minimum value: 2, do not remove neighbours: 0 or 1)
NEIGHBOURS_MINDISTANCE=2
# KEEP
# Score threshold below which predicted contacts (DMP) will be ignored (def: 0.2, do not remove contacts: 0.0)
DMPSCORE_THRESHOLD=0.0
# KEEP
# Score threshold below which predicted contacts (PSICOV) will be ignored (def: 0.2 *doi: 10.1038/s41598-019-48913-8, do not remove contacts: 0.0)
PSICOVSCORE_THRESHOLD=0.2
# KEEP
# Remove intermolecular contacts that are also intramolecular contacts too (bool) (If True, NEIGHBOURS_MINDISTANCE will be set to 2; Intramolecular contacts determine this threshold)
REMOVE_INTRA_CONTACTS=True
# ADD TO MAIN

########## LOG-STDOUT
# If stdout is to be stored in a log-file (True) or shown on screen (False)
# LOG_STDOUT=True
# ALWAYS TRUE (MAIN)