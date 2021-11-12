"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

import csv

def csvheader(outpath):
    with open(outpath, "w") as outcsv:
        outcsv.write('# PDBid, Mode, Interface, Lseq, SEQDepth, N_contacts, ' +
                     '1L_Avscore, 1L_Accscore, 1L_TPs, 1L_Jaccard, ' +
                     '2L_Avscore, 2L_Accscore, 2L_TPs, 2L_Jaccard\n')