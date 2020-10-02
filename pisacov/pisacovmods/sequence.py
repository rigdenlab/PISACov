# -*- coding: utf-8 -*-

import os
import gemmi
import csv
from conkit import io as ckio
from conkit import core as ckc


from .inputvalues import CSV_CHAIN_PATH
from .init2 import pdbid
from .init2 import ressymbol
from .output import output_tmpdir
from .output import printout
#########################################################
##### SEQUENCE-RELATED FUNCTIONS ########################

# INTERFACE PDB CROPPING AND RENUMBERING
def renumberpdbs(PDBINT_PATH, fastaseq,itisbio=False,outdir=output_tmpdir("pisacov")): # CROPS
    """
    INTERFACE PDB CROPPING AND RENUMBERING

    Parameters
    ----------
    PDBINT_PATH : str
        Input path.
    fastaseq : ConKit sequence # CROPS
        Full fasta sequence from PDB database
    itisbio : bool
        Use limits given by SIFTS (True) or whole sequence (False) (def: False)

    Returns
    -------
    None

    """
    n_chains = 0
    n_resmax = 0
    pdb_structure = gemmi.read_structure(PDBINT_PATH)

    for model in pdb_structure:
        n_chains += len(model)
        for chain in model:
            if len(chain) > n_resmax:
                n_resmax = len(chain)

    pos = [[0 for j in range(n_resmax)] for i in range(n_chains)]

    n_chains = 0

    ninterface=(os.path.splitext(os.path.splitext(PDBINT_PATH)[0])[1])[1:]
    printout('INTERFACE: ' + str(ninterface))
    for model in pdb_structure:
        for chain in model:
            solved = False
            for shift in range(int(len(chain)/2)):
                cnt=0
                gap=0
                score=0
                newseq=''
                newseq += '-'*shift
                for residue in chain:
                    if residue == chain[0]:
                        if ressymbol(residue.name) == fastaseq.seq[shift]:# CROPS
                            score += 1
                            pos[n_chains][cnt]=1+shift
                            newseq += ressymbol(residue.name)
                    else:
                        if (chain[cnt].seqid.num-chain[cnt-1].seqid.num > 1):
                            gap += (chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                            newseq += '-'*(chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                        pos[n_chains][cnt]=cnt+1+gap+shift
                        if ressymbol(residue.name) == fastaseq.seq[cnt+gap+shift]:# CROPS
                            score += 1
                            newseq += ressymbol(residue.name)
                        if residue==chain[-1]:
                            if cnt+gap+shift+1 < len(fastaseq.seq):# CROPS
                                newseq += '-'*(len(fastaseq.seq)-(cnt+gap+shift+1))# CROPS
                    cnt += 1
                if score == len(chain):
                    solved = True
                    printout('ALIGNMENT OF CHAIN ' + str(n_chains+1))
                    printout(fastaseq.seq)# CROPS
                    printout(newseq, extraline=True)
                    break
            if solved:
                cnt=0
                for residue in chain:
                    residue.seqid.num = pos[n_chains][cnt]
                    cnt += 1
            n_chains += 1
            solved = False

    PDBINT_PATH_OUT=os.path.join(outdir, os.path.splitext(os.path.basename(PDBINT_PATH))[0]+".fasta.pdb")
    pdb_structure.write_pdb(PDBINT_PATH_OUT)

    if itisbio:
        PDBINT_PATH_OUT=os.path.join(outdir, os.path.splitext(os.path.basename(PDBINT_PATH))[0]+".bio.fasta.pdb")

        n_chains = 0
        n_resmax = 0
        for model in pdb_structure:
            n_chains += len(model)
            for chain in model:
                if len(chain) > n_resmax:
                    n_resmax = len(chain)
        delres = [[False for j in range(n_resmax)] for i in range(n_chains)]
        n_chains = 0
        for model in pdb_structure:
            for chain in model:
                r=0
                newresnum=1
                fastaends=biochain_ends('fasta')
                for residue in chain:
                    if residue.seqid.num < fastaends[0]:
                        delres[n_chains][r] = True
                    elif residue.seqid.num > fastaends[1]:
                        delres[n_chains][r] = True
                    elif residue.seqid.num == fastaends[0]:
                        prev_res_num = residue.seqid.num
                        residue.seqid.num = 1
                    else:
                        if residue == chain[0]:
                            prev_res_num = residue.seqid.num
                            residue.seqid.num = residue.seqid.num - fastaends[0] + 1
                        else:
                            if (residue.seqid.num - prev_res_num == 0 ): # SEQUENCE NUMBERS INSERTED (1A, 1B, 1C, ...)
                                newresnum += 1
                                prev_res_num = residue.seqid.num
                                residue.seqid.num = newresnum
                            elif (residue.seqid.num - prev_res_num == 1 ): # SEQUENCE NUMBERS ARE CONSECUTIVE
                                newresnum += 1
                                prev_res_num = residue.seqid.num
                                residue.seqid.num = newresnum
                            elif (residue.seqid.num - prev_res_num > 1 ): # SEQUENCE NUMBERS SHOW GAPS
                                newresnum += residue.seqid.num- prev_res_num
                                prev_res_num = residue.seqid.num
                                residue.seqid.num = newresnum
                            else:
                                print('residue : ' + str(residue.seqid.num))
                                print('previous: ' + str(prev_res_num))
                                printout(' ERROR: Sequence numbers not sorted',extraline=True, errorlog=True)
                    r += 1
                n_chains += 1

        n_chains = 0
        for model in pdb_structure:
            for chain in model:
                for res in reversed(range(len(chain))):
                    if delres[n_chains][res]:
                        del chain[res]
                n_chains += 1

        pdb_structure.write_pdb(PDBINT_PATH_OUT)

def biochain_ends(which): # CROPS
    """
    RETRIEVE SEQUENCE POSITION OF BIOLOGICAL PROTEIN CHAIN ENDS FROM CSV FILE (SIFTS)

    Parameters
    ----------
    which : str
        Either 'fasta' or 'pdb'

    Returns
    -------
    resends : list(int), dim=2
        The position (sequence numbers) of the biological ends at the pdb file

    """

    #pdbid=os.path.splitext(os.path.basename(PDB_PATH))[0]
    csv_chain_file = open(CSV_CHAIN_PATH)
    csv_chain = csv.reader(csv_chain_file)

    resends = [0 for i in range(2)]

    for entry in csv_chain:
        if entry[0] == pdbid() and entry[1]=="A":
            if which == 'fasta':
                resends[0]=int(entry[3])
                resends[1]=int(entry[4])
            elif which == 'pdb':
                resends[0]=int(entry[7])
                resends[1]=int(entry[8])
            break
    return resends


def crop_fasta(fastaseq, outdir=output_tmpdir()): # CROPS
    """
    FASTA SEQUENCE CROPPING

    Parameters
    ----------
    seqpath : ConKit Sequence# CROPS
        Source fasta sequence
    outdir : str, optional
        Directory where results biological sequence will be printed out. The default is output_tmpdir().

    Returns
    -------
    bioseq : ConKit Sequence
        Cropped (biological) sequence
    newseqpath : str
        Sequence path (bio.fasta file)

    """
    # Obtain new chain ends (residue number)
    fastaends = biochain_ends('fasta')

    # Check that the sequence is consistent with the limits retrieved from the database
    if fastaends[1]-fastaends[0] + 1 > fastaseq.seq_len:# CROPS
        isitbio=False
        printout('WARNING: The biological sequence limits include a section greater than the input sequence.', errorlog=True)
        printout('         Skipping cropping. Returning input values.', errorlog=True,extraline=True)
        bioseq = fastaseq# CROPS
        newseqfile = pdbid() + '.fasta'
        newseqpath = os.path.join(outdir, newseqfile)
    elif fastaends[1]-fastaends[0] + 1 == fastaseq.seq_len:# CROPS
        isitbio=True
        printout('         Biological and input sequences have the same length. Skipping cropping. Returning input values.',extraline=True)
        bioseq = fastaseq# CROPS
        newseqfile = pdbid() + '.bio.fasta'
        newseqpath = os.path.join(outdir, newseqfile)
        ckio.write(newseqpath,"fasta",hierarchy=bioseq)
    else:
        if fastaends[0] > fastaseq.seq_len:# CROPS
            isitbio=False
            printout('WARNING: The sequence upper limit imported from the database is higher than the upper limit from the fasta file.',errorlog=True) #LOGGING
            printout('         Skipping cropping. Returning input values.', errorlog=True,extraline=True)
            bioseq = fastaseq# CROPS
            newseqfile = pdbid() + '.fasta'
            newseqpath = os.path.join(outdir, newseqfile)
        else:
            isitbio=True
            # Append new info to sequence
            newid=fastaseq.id# CROPS
            newid = newid +"|NO_CLONING_ARTIFACTS"

            # Create new sequence
            newseq=fastaseq.seq[fastaends[0]-1:fastaends[1]-1]# CROPS
            bioseq=ckc.Sequence(newid, newseq) # CROPS

            # Write new sequence to file
            #locpdbid=os.path.splitext(os.path.basename(PDB_PATH))[0]
            #outdir = os.path.join(OUTPUT_DIR, locpdbid,"")

            newseqfile = pdbid() + '.bio.fasta'
            newseqpath = os.path.join(outdir, newseqfile)

            ckio.write(newseqpath,"fasta",hierarchy=bioseq)

    return bioseq, newseqpath, isitbio

def biofile(itisbio):
    """
    Returns string containing ".bio" or empty string depending on fasta sequence employed

    Parameters
    ----------
    itisbio : bool
        Contains information about the nature of fasta sequence.

    Returns
    -------
    bio_path : str
        Either ".bio" or empty string.

    """
    if itisbio:
        bio_path='.bio'
    else:
        bio_path=''

    return bio_path