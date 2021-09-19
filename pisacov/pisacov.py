#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#############################################
## Python 2.7 and 3.x versions compatible
#############################################

import pisacovmods.inputvalues as pmin
import pisacovmods.init1 as pmi1
import pisacovmods.init2 as pmi2
import pisacovmods.output as pmo
import pisacovmods.pisacovmod as pmp
import pisacovmods.sequence as pms
import pisacovmods.contacts as pmc

from pisacov.about import __prog__, __description__, __version__
from pisacov.about import  __author__, __date__, __copyright__

from pisacov import command_line as pcl

import argparse
import conkit
import os
import shutil
import matplotlib.pyplot as plt
# from string import format

import time

logger = None

def create_argument_parser():
    """Create a parser for the command line arguments used in crops-renumber"""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)
    parser.add_argument("input_seqpath",nargs=1, metavar="Original_RCSB_Sequence_filepath",
                        help="Input sequence filepath as downloaded from RCSB.")
    parser.add_argument("input_strpath",nargs=1, metavar="Original_RCSB_Structure_filepath",
                        help="Input structure filepath or dir as downloaded from RCSB. If a directory is inserted, it will act on all structure files in such directory.")

    parser.add_argument("-o","--outdir",nargs=1,metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")

    parser.add_argument("-i","--intervals",nargs=1, metavar="Intervals_database",
                        help="Override input intervals database filepath from installation file.")


    sections=parser.add_mutually_exclusive_group(required=False)
    sections.add_argument("-t","--terminals",action='store_true',default=False,
                          help="Ignore interval discontinuities and only crop the ends off.")
    sections.add_argument("-u","--uniprot_threshold", nargs=2, metavar=("Uniprot_ratio_threshold","Sequence_database"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. [MIN,MAX)=[0,100) uniprot_sprot.fasta-path')

    parser.add_argument('--version', action='version', version='%(prog)s '+ __version__)

def main():

    starttime=time.time()
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.crops_logger(level="info")
    logger.info(pcl.welcome())

    inseq=check_path(args.input_seqpath[0],'file')
    indb=check_path(args.input_database[0],'file')
    insprot=check_path(args.uniprot_threshold[1]) if args.uniprot_threshold is not None else None

    minlen=float(args.uniprot_threshold[0]) if args.uniprot_threshold is not None else 0.0
    targetlbl=ctg.target_format(indb,terms=args.terminals, th=minlen)
    infixlbl=ctg.infix_gen(indb,terms=args.terminals)

    if args.outdir is None:
        outdir=check_path(os.path.dirname(inseq),'dir')
    else:
        outdir=check_path(os.path.join(args.outdir[0],''),'dir')

    if args.sort is not None:
        if (args.sort[0].lower()!='ncrops' and args.sort[0].lower()!='percent' and
            args.sort[0].lower()!='ncropsin' and args.sort[0].lower()!='percentin'):
            raise ValueError("Arguments for sorting option can only be either 'ncrops' or 'percent'.")
        else:
            sorter=args.sort[0].lower()

#############################################################

    sources=["deepmetapsicov", "psicov"]
    confiledir=["deepmetapsicov", "deepmetapsicov"]
    confilesuffix=["deepmetapsicov.con","psicov"]
    n_sources=len(sources)

# BEGIN: DEBUGGING ONLY

# Skip the execution of PISA, HHBLITS and DeepMetaPSICOV (?) [bool]
#    SKIP_EXEC=[True, True, True]
    SKIP_EXEC=[False, False, False]

# Skip the creation of new dir (?) [bool]
#    SKIP_MKDIR=True
    SKIP_MKDIR=False

# END: DEBUGGING ONLY


# Create output directory
    if not SKIP_MKDIR:
        pmo.mkdirout()
        pmo.printout('Output directory created',extraline=True)

    starttime=pmi2.gettime()
    pmo.printout('Starting Time: '+ starttime[1].strftime("%-d %B %Y, %X") + ' UTC\n\n')

# Check that input values are correct

    hhparameters=pmi1.hhparam('logit')
    dumvar=pmi1.minneigh('logit')
    dumvar=pmi1.scorethreshold("deepmetapsicov",'logit')
    dumvar=pmi1.scorethreshold("psicov",'logit')

#####################################################################
##### MSA generator #################################################

    pmo.printout('*********************************************************')
    pmo.printout('*** MULTIPLE SEQUENCE ALIGNMENT *************************',extraline=True)

    pmo.printout('Importing RCSB PDB sequence from fasta file...')
    fasta_seq=conkit.io.read(pmin.SEQUENCE_PATH,"fasta")[0]
    pmo.printout(fasta_seq)

    if pmin.USE_BIOCHAIN:
        pmo.printout('Obtaining Biological sequence from fasta file...')
        if not SKIP_EXEC[1]:
            biological_seq, seqpath, bio = pms.crop_fasta(fasta_seq,pmo.output_tmpdir("deepmetapsicov"))
            if not bio:
                shutil.copyfile(pmin.SEQUENCE_PATH, seqpath)
        else:
            tmpfile = pmi2.pdbid() + '.bio.fasta'
            seqpath = os.path.join(pmo.output_tmpdir("deepmetapsicov"), tmpfile)
            try:
                biological_seq=conkit.io.read(seqpath,"fasta")[0]
                bio=True
            except:
                bio=False
                tmpfile = pmi2.pdbid() + '.fasta'
                seqpath = os.path.join(pmo.output_tmpdir("deepmetapsicov"), tmpfile)
                biological_seq=conkit.io.read(seqpath,"fasta")[0]
        seq=biological_seq
        tmpfile = pmi2.pdbid() + pms.biofile(bio)+ '.fasta'
        newseqpath = os.path.join(pmo.output_dir(), tmpfile)
        shutil.copyfile(seqpath, newseqpath)
        if bio:
            pmo.printout('Biological sequence:')
            pmo.printout(biological_seq)
        else:
            pmo.printout('WARNING: Biological Sequence not available. Using RCSB PDB fasta sequence for analysis (including cloning artifacts)',errorlog=True, extraline=True) # LOGGING
    else:
        bio=False
        pmo.printout('Using RCSB PDB fasta sequence for analysis (including cloning artifacts)', extraline=True)
        seqpath = os.path.join(pmo.output_dir(), os.path.splitext(os.path.basename(pmin.SEQUENCE_PATH))[0])+ os.path.splitext(os.path.basename(pmin.SEQUENCE_PATH))[1]
        shutil.copyfile(pmin.SEQUENCE_PATH, seqpath)
        seqpath = os.path.join(pmo.output_tmpdir("deepmetapsicov"), os.path.splitext(os.path.basename(pmin.SEQUENCE_PATH))[0])+ os.path.splitext(os.path.basename(pmin.SEQUENCE_PATH))[1]
        shutil.copyfile(pmin.SEQUENCE_PATH, seqpath)
        seq=fasta_seq

    pmo.printout('')
    if pmin.HHBLITS_VIA_DMP:
        pmo.printout("Using DeepMetaPSICOV's execution of HHblits...")
        msapath=''
    else:
        if not SKIP_EXEC[1]:
            pmo.printout('Creating Multiple Sequence Alignment with HHblits...')
            msa, msapath = pmp.runhhblits(bio,param=hhparameters,spath=seqpath)
        else:
            pmo.printout('Reading Multiple Sequence Alignment...')
            msafile = pmi2.pdbid()+pms.biofile(bio)+".msa.aln"
            msapath = os.path.join(pmo.output_dir(), msafile)
            msa = conkit.io.read(msapath,'jones')

    msacovpath=os.path.splitext(msapath)[0]+".coverage.png"
    msaformat='jones'
    sit=0.7

#####################################################################
##### DeepMetaPSICOV: Contact Prediction ############################
    pmo.printout('*********************************************************')
    pmo.printout('*** CONTACT PREDICTION **********************************',extraline=True)

    if not SKIP_EXEC[2]:
        pmo.printout('Running DeepMetaPSICOV for Contact prediction list...')
    else:
        pmo.printout('Skipping DeepMetaPSICOV execution.')

    if pmin.HHBLITS_VIA_DMP:
        msa, msapath = pmp.rundmp(seqpath, msapath,skiphhblits=SKIP_EXEC[1],skipdmp=SKIP_EXEC[2]) # Input Sequence fasta file path, MSA file path
    else:
        if not SKIP_EXEC[2]:
            pmp.rundmp(seqpath, msapath) # Input Sequence fasta file path, MSA file path


#####################################################################
##### PISA : OBTAIN INTERFACE PDB FILES #############################
    pmo.printout('*********************************************************')
    pmo.printout('*** INTERFACE IDENTIFICATION ****************************',extraline=True)

# Obtain Number of interfaces in PDB and produce Interface PDB files (including renumbered)
    if not SKIP_EXEC[0]:
        pmo.printout('Running PISA...')
        n_interfaces = pmp.runpisa(bio)
        for nif in range (n_interfaces):
            pdbintpath = os.path.join(pmo.output_tmpdir("pisa"), pmi2.pdbid()+".interface."+str(nif+1)+".pdb")
            pms.renumberpdbs(pdbintpath, fasta_seq, bio)
    else:
        pmo.printout('Reading PISA xml file...')
        n_interfaces = pmp.n_int_xml()


#####################################################################
##### Interface Confidence Scores ###################################
    pmo.printout('*********************************************************')
    pmo.printout('*** SCORING INTERFACES **********************************',extraline=True)

    pmo.printout('Importing contact prediction file ...',extraline=True)
    conpred=[]
    conpred_id=[]
    conpredpath=[]

    for n in range(n_sources):
        conpredfile=pmi2.pdbid()+pms.biofile(bio)+"."+confilesuffix[n]

        conpredpath.append(os.path.join(pmo.output_tmpdir(confiledir[n]), conpredfile))

        conpred.append(conkit.io.read(conpredpath[n], 'psicov')[0])
        conpred[n].sequence=seq
        conpred[n].set_sequence_register()
        conpred_id.append(conpred[n].id)

    pmo.printout('Remove excluded contact scores ...',extraline=True)

    for n in range(n_sources):
        conpred[n]=pmc.filter_contacts(conpred[n],sources[n])

    cntint=0
    cntlig=0
    cntunk=0
    interfacetype=[0 for i in range(n_interfaces+1)]
    scores = [[[0 for k in range(4)] for j in range(n_sources)] for i in range(n_interfaces)]
    n_contacts_all = [[[0 for k in range(4)] for j in range(n_sources)] for i in range(n_interfaces)]
    ndec=6

    for nif in range(n_interfaces):
        pmo.printout('Importing Interface '+str(nif+1)+' PDB file...')
        maps=pmc.import_interfacepdb(nif,bio)
        conpredInt=[]
        pmo.printout('  Number of maps : ' + str(len(maps)) )
        if (len(maps)==4):
            for n in range(n_sources):
                interfacetype[nif]="MM"
                conpredInt.append(conpred[n].deepcopy())
                if pmin.REMOVE_INTRA_CONTACTS:
                    pmo.printout("  Removing those contacts from interface that also appear as intramolecular contacts ...")
                    for element in [0,3]:
                        conpredInt[n]=pmc.remove_intra_contacts(conpredInt[n],maps[element])

                n_contacts_all[nif][n][0]=conpredInt[n].ncontacts # Total number of contacts predicted for interface
                n_contacts_all[nif][n][1]=maps[1].ncontacts # Total number of contacts from interface PDB file

            #pmo.printout(conpredInt)
            pmo.printout(maps[1])
            cntint += 1

            # Write contact list
            pmo.printout('  Writing contact lists ...' )
            if pmin.REMOVE_INTRA_CONTACTS:
                for element in [0,3]:
                    chainid = str(1) if element == 0 else str(2)
                    confile=pmi2.pdbid()+".pdb.interface."+str(nif+1)+".intrachain"+chainid+".conkit.con"
                    tmppath = os.path.join(pmo.output_tmpdir("pisacov"), confile)
                    conkit.io.write(tmppath, 'psicov', hierarchy=maps[element])
            confile=pmi2.pdbid()+".pdb.interface."+str(nif+1)+".conkit.con"
            tmppath = os.path.join(pmo.output_tmpdir("pisacov"), confile)
            conkit.io.write(tmppath, 'psicov', hierarchy=maps[1])

            pmo.printout('  Matching contact prediction and interface contact lists ...',extraline=True )
            contactpath=os.path.join(pmo.output_tmpdir("pisacov"), confile)
            for n in range(n_sources):
                conpredInt_pdb=pmc.pred_pdb_matchlist(conpredInt[n],contactpath)
                conpredInt_pdb.id="     Matching contacts of " + conpred_id[n]
                conpredInt_pdb.sequence = seq.deepcopy()
                conpredInt_pdb.set_sequence_register()

                if conpredInt_pdb.ncontacts == 0:
                    pmo.printout("  WARNING: No contacts found in the filtered interface "+str(nif)+" pdb file. SKIP this result type.",errorlog=True)
                    for sc in range(4):
                        scores[nif][n][sc]="***"
                    n_contacts_all[nif][n][3] = "***"
                else:
                    map_matched, n_contacts_all[nif][n][2], scores[nif][n][0], scores[nif][n][1], scores[nif][n][2], scores[nif][n][3] = pmc.match_maps(conpredInt_pdb, maps[1])
                    for sc in range(4):
                        scores[nif][n][sc]=round(scores[nif][n][sc],ndec)

                    n_contacts_all[nif][n][3] = n_contacts_all[nif][n][2] / n_contacts_all[nif][n][1]
                    ## Plot matched map
                    pngfile=pmi2.pdbid()+".interface."+str(nif+1)+"."+sources[n]+".con.png"
                    confile=pmi2.pdbid()+"."+sources[n]+".interface."+str(nif+1)+".conkit.con"

                    pngpath = os.path.join(pmo.output_dir(), pngfile)
                    fig = conkit.plot.ContactMapFigure(map_matched, reference=maps[1])
                    fig.savefig(pngpath, overwrite=True)

                    tmppath = os.path.join(pmo.output_dir(), confile)
                    conkit.io.write(tmppath, 'psicov',hierarchy=conpredInt_pdb)
                    plt.close('all')

        elif (len(maps)==2):
            interfacetype[nif]="Unk"
            pmo.printout("  WARNING: Unexpected number of maps (2)", errorlog=True)
            pmo.printout('----> REVISE CONTACT MAPS FOR THIS INTERFACE', errorlog=True,extraline=True)

            cntunk += 1
            cnt=0
            for score in range(4):
                for n in range(n_sources):
                    scores[nif][n][score]="***"
            for mapn in maps:
                confile=pmi2.pdbid()+".pdb."+str(cnt+1)+".interface."+str(nif+1)+".conkit.con"  #### CHECK FILES CREATED WHEN FIXING THIS SECTION.
                tmppath = os.path.join(pmo.output_dir(), confile)
                conkit.io.write(tmppath, 'psicov', hierarchy=maps[1])
                cnt += 1
        elif (len(maps)==1):
            interfacetype[nif]="LM"
            pmo.printout("  ... Ligand-monomer interface detected")
            pmo.printout('      SKIP', extraline=True)

            cntlig += 1
            for score in range(4):
                for n in range(n_sources):
                    scores[nif][n][score]="***"
        else:
            interfacetype[nif]="Unk"
            pmo.printout("  WARNING: Unexpected number of maps", errorlog=True)
            pmo.printout('      SKIP', extraline=True, errorlog=True)

            cntunk += 1
            for score in range(4):
                for n in range(n_sources):
                    scores[nif][n][score]="***"



#####################################################################
##### OUTPUT ########################################################

    pmo.printout('*********************************************************')
    pmo.printout('*** FINAL OUTPUT ****************************************',extraline=True)

    pmo.printout('Printing final data to file...', extraline=True)
    datfile=pmi2.pdbid()+".pisacov.out.dat"
    datpath = os.path.join(pmo.output_dir(), datfile)
    now = pmi2.gettime()
    space=str(ndec+4)
    with open(datpath, 'w') as out:
        if bio:
            out.write( "# PDB id: "+ pmi2.pdbid() + ' - Contact def: 0<d<8 Angstroms - Cloning artifacts removed from fasta file sequence - ' )
        else:
            out.write( "# PDB id: "+ pmi2.pdbid() + ' - Contact def: 0<d<8 Angstroms - Original fasta file sequence used - ' )
        out.write( now[1].strftime("%-d %B %Y, %X") + ' UTC')
        out.write('\n')
        datahead='# IF_id IF_type n_pdb  '
        for n in range(n_sources):
            datahead += 'n_'+sources[n]+'  '
            datahead += 'n_'+sources[n]+'/n_pdb  '
            datahead += 'Av_score_'+sources[n]+'  '
            datahead += 'Acc_score_'+sources[n]+'  '
            datahead += 'P_'+sources[n]+'  '
        out.write(datahead+'\n')
        for nif in range(n_interfaces):
            outformat='{:>4s}{:>4s}{:>4s}'
            outstring = str(outformat.format(str(nif+1),interfacetype[nif],str(n_contacts_all[nif][0][1])))
            for n in range(n_sources):
                #outformat='{:>4s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}'
                outformat='{:>4s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}{:>'+space+'s}'
                outstring += str(outformat.format(str(n_contacts_all[nif][n][2]),str(scores[nif][n][0]),str(scores[nif][n][1]),str(scores[nif][n][2]),str(scores[nif][n][3])))
            out.write(outstring+'\n')
    sumlogfile=pmi2.pdbid()+".pisacov.out.summary.log"
    sumlogpath = os.path.join(pmo.output_dir(), sumlogfile)
    now = pmi2.gettime()

    with open(sumlogpath, 'w') as out:
        out.write('*********************************************************\n')
        out.write('*** S U M M A R Y  **  P I S A C O V ********************\n')
        out.write('*********************************************************\n\n')
        out.write(now[1].strftime("%-d %B %Y, %X") + ' UTC\n\n')
        out.write('Protein PDB ID:  '+pmi2.pdbid()+'\n\n')
        out.write('*********************************************************\n')
        out.write('--- Multiple Sequence Alignment ---\n\n')
        outformat='{:<55s}{:<80s}'
        out.write(outformat.format('   MSA file: ', msapath))
        out.write('\n')
        out.write(outformat.format('   MSA format: ', msaformat))
        out.write('\n')
        out.write(outformat.format('   Sequence Identity Threshold: ', str(sit) ))
        out.write('\n')
        out.write(outformat.format('   Length of the Target Sequence: ',str(msa.top_sequence.seq_len)))
        out.write('\n')
        out.write(outformat.format('   Total number of sequences: ', str(msa.nseq)))
        out.write('\n')
        out.write(outformat.format('   Number of Effective Sequences: ', str(msa.meff)))
        out.write('\n')
        out.write(outformat.format('   Proportion of Effective Sequences: ', str(round(100*msa.meff/msa.nseq,2))+' %'))
        out.write('\n')
        out.write(outformat.format('   Sequence Coverage Plot: ', msacovpath))
        out.write('\n\n')
        out.write(outformat.format('   MSA created with: ',pmin.HHSUITE_PATH))
        out.write('\n')
        out.write(outformat.format('   Reference database name: ', pmin.HHBLITS_DATABASE_NAME))
        out.write('\n')
        out.write(outformat.format('   Reference database path: ', pmin.HHBLITS_DATABASE_DIR))
        out.write('\n')
        out.write(outformat.format('   Number of iterations: ', str(pmi1.hhparam()[0])))
        out.write('\n')
        out.write(outformat.format('   E-value cutoff for inclusion in result alignment: ', str(pmi1.hhparam()[1])))
        out.write('\n')
        out.write(outformat.format('   Non-redundant sequences to keep: ', str(pmi1.hhparam()[2])))
        out.write('\n')
        out.write(outformat.format('   Minimum coverage with master sequence (%): ', str(pmi1.hhparam()[3])))
        out.write('\n')
        out.write(outformat.format('   Maximum pairwise sequence identity: ', str(pmi1.hhparam()[4])))
        out.write('\n\n')
        out.write('*********************************************************\n')
        out.write('--- Contact Prediction ---\n\n')
        outformat='{:<55s}{:<80s}'
        out.write(outformat.format('   Contact prediction list(s) created with: ', pmin.DMP_PATH))
        out.write('\n')
        for n in range(n_sources):
            out.write(outformat.format('   Contact prediction file ('+sources[n]+'): ', conpredpath[n]))
            out.write('\n')
        out.write(outformat.format('   Minimum distance within sequence (neigh. cutoff): ', str(pmi1.minneigh()) ))
        out.write('\n')
        for n in range(n_sources):
            out.write(outformat.format('   Contact prediction score threshold ('+sources[n]+'): ', str(pmi1.scorethreshold(sources[n]))))
            out.write('\n')
        out.write('\n')
        out.write('*********************************************************\n')
        out.write('--- Interfaces ---\n\n')
        outformat='{:<55s}{:<80s}'
        out.write(outformat.format('   Total Number of Interfaces: ', str(n_interfaces) ))
        out.write('\n')
        out.write(outformat.format('   Number of Intramolecular interfaces: ', str(cntint) ))
        out.write('\n')
        out.write(outformat.format('   Number of Ligand-monomer interfaces: ', str(cntlig) ))
        out.write('\n')
        out.write(outformat.format('   Number of Unidentified interfaces: ', str(cntunk) ))
        out.write('\n\n')
        for nif in range(n_interfaces):
            nif1=nif+1
            out.write('   --- Interface '+str(nif1)+' ---\n')
            if (interfacetype[nif] == "MM"):
                out.write(outformat.format('      Interface type:','monomer - monomer (MM)'))
                out.write('\n')
                out.write(outformat.format('      Total number of intermolecular contacts (pdb): ', str(n_contacts_all[nif][0][1])))
                out.write('\n')
                for n in range(n_sources):
                    extral='\n\n' if n==n_sources-1 else '\n'
                    out.write('      + '+sources[n]+' Scores +\n')
                    out.write(outformat.format('        Number of True Positives : ', str(n_contacts_all[nif][n][2]) ))
                    out.write('\n')
                    out.write(outformat.format('        Proportion of True positives: ', str(n_contacts_all[nif][n][3])))
                    out.write('\n')
                    out.write(outformat.format('        Jaccard Index (|PDB ⋂ pred| / |PDB U pred|): ', str(scores[nif][n][0])))
                    out.write('\n')
                    out.write(outformat.format('        Average value of the scores of True Positives: ', str(scores[nif][n][1])))
                    out.write('\n')
                    out.write(outformat.format('        Sum of all the scores of True Positives: ', str(scores[nif][n][2])))
                    out.write('\n')
                    out.write(outformat.format('        Probabilistic score for whole interface: ', str(scores[nif][n][3])))
                    out.write(extral)

            elif (interfacetype[nif] == "LM"):
                out.write(outformat.format('      Interface type:','ligand - monomer (LM)'))
                out.write('\n\n')
            elif (interfacetype[nif] == "Unk"):
                out.write(outformat.format('      Interface type:','unidentified (Unk)'))
                out.write('\n\n')


    pmo.printout('*********************************************************')
    pmo.printout('*** S U M M A R Y ***************************************')
    pmo.printout('*********************************************************', extraline=True)

    pmo.printout('Protein PDB ID: %s' % pmi2.pdbid(), extraline=True)

    pmo.printout('--- Multiple Sequence Alignment ---')
    pmo.printout('   MSA file: ' + msapath)
    pmo.printout('   MSA format: ' + msaformat)
    pmo.printout('   Sequence Identity Threshold: '+ str (sit))
    pmo.printout('   Length of the Target Sequence: ' + str(msa.top_sequence.seq_len))
    pmo.printout('   Total number of sequences: ' + str(msa.nseq))
    pmo.printout('   Number of Effective Sequences: ' + str(msa.meff))
    pmo.printout('   Proportion of Effective Sequences: ' + str( round(100*msa.meff/msa.nseq,2))+' %')
    pmo.printout('   Sequence Coverage Plot: ' + msacovpath, extraline=True)

    pmo.printout('   MSA created with: ' + pmin.HHSUITE_PATH)
    pmo.printout('   Reference database name: ' + pmin.HHBLITS_DATABASE_NAME)
    pmo.printout('   Reference database path: ' + pmin.HHBLITS_DATABASE_DIR)
    pmo.printout('   Number of iterations: ' + str(pmi1.hhparam()[0]))
    pmo.printout('   E-value cutoff for inclusion in result alignment: ' + str(pmi1.hhparam()[1]))
    pmo.printout('   Non-redundant sequences to keep: ' + str(pmi1.hhparam()[2]))
    pmo.printout('   Minimum coverage with master sequence (%): ' + str(pmi1.hhparam()[3]))
    pmo.printout('   Maximum pairwise sequence identity: ' + str(pmi1.hhparam()[4]), extraline=True)

    pmo.printout('--- Contact Prediction ---')
    pmo.printout('   Contact prediction list created with: ' + pmin.DMP_PATH)
    for n in range(n_sources):
        pmo.printout('   Contact prediction file ('+sources[n]+'): ' + conpredpath[n])
    pmo.printout('   Minimum distance within sequence (neigh. cutoff): '+ str(pmi1.minneigh()))
    for n in range(n_sources):
        extral=True if n==n_sources-1 else False
        pmo.printout('   Contact prediction score threshold ('+sources[n]+'): ' + str(pmi1.scorethreshold(sources[n])),extraline=extral)

    pmo.printout('--- Interfaces ---')
    pmo.printout('   Number of Interfaces: ' + str(n_interfaces) + ", of which:")
    pmo.printout('     Number of Intramolecular interfaces: ' + str(cntint))
    pmo.printout('     Number of Ligand-monomer interfaces: ' + str(cntlig))
    pmo.printout('     Number of Unidentified interfaces: ' + str(cntunk),extraline=True)

    for nif in range(n_interfaces):
        nif1=nif+1
        pmo.printout('  --- Interface '+str(nif1)+' ---',extraline=True)

        if (interfacetype[nif] == "MM"):
            pmo.printout('   Interface type: monomer - monomer (MM)')
            pmo.printout('       Total number of intermolecular contacts (pdb): ' + str(n_contacts_all[nif][0][1]))
            for n in range(n_sources):
                extral=True if n==n_sources-1 else False
                pmo.printout('       + Scores ('+sources[n]+'):')
                pmo.printout('         Number of True Positives : ' + str(n_contacts_all[nif][n][2]))
                pmo.printout('         Proportion of True positives: ' + str(n_contacts_all[nif][n][3]) )
                pmo.printout('         Jaccard Index (|PDB ⋂ pred| / |PDB U pred|): ' + str(scores[nif][n][0]) )
                pmo.printout('         Average value of the scores of True Positives: ' + str(scores[nif][n][1] ))
                pmo.printout('         Sum of all the scores of True Positives: ' + str(scores[nif][n][2] ))
                pmo.printout('         Probabilistic score for whole interface: ' + str(scores[nif][n][3] ),extraline=extral)

        elif (interfacetype[nif] == "LM"):
            pmo.printout('   Interface type: ligand - monomer (LM)',extraline=True)
        elif (interfacetype[nif] == "Unk"):
            pmo.printout('   Interface type: unidentified (Unk)',extraline=True)

#####################################################################
##### TIMEIT ########################################################

    pmo.printout('*********************************************************')
    pmo.printout('*** T I M I N G *****************************************')
    pmo.printout('*********************************************************', extraline=True)

    with open(sumlogpath, 'a') as out:
        out.write('*********************************************************\n')
        out.write('*** T I M I N G *****************************************\n')
        out.write('*********************************************************\n\n')

        endtime=pmi2.gettime()
        totaltime=pmi2.readabletime(endtime[0]-starttime[0])

        out.write('Starting Time: '+ starttime[1].strftime("%-d %B %Y, %X") + ' UTC\n\n')
        out.write('Ending Time: '+ endtime[1].strftime("%-d %B %Y, %X") + ' UTC\n\n')
        out.write('Process Wallclock time: '+str(endtime[0]-starttime[0])+' s\n')
        out.write('      or, equivalently: ' + totaltime+'\n\n')

    pmo.printout('Starting Time: '+ starttime[1].strftime("%-d %B %Y, %X")+' UTC', extraline=True)
    pmo.printout('Ending Time: '+ endtime[1].strftime("%-d %B %Y, %X")+' UTC', extraline=True)
    pmo.printout('Process Wallclock time: '+str(endtime[0]-starttime[0])+' s')
    pmo.printout('      or, equivalently: ' + totaltime, extraline=True)


if __name__ == "__main__":
    import sys
    #import traceback

    try:
        main()
        sys.exit(0)
    except:
        sys.exit(1)
    #except Exception as e:
    #    if not isinstance(e, SystemExit):
    #        msg = "".join(traceback.format_exception(*sys.exc_info()))
    #        logger.critical(msg)
    #    sys.exit(1)
