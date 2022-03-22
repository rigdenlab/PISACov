"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

# import pisacovmods.inputvalues as pmin
# import pisacovmods.init1 as pmi1
# import pisacovmods.init2 as pmi2
# import pisacovmods.output as pmo
# import pisacovmods.pisacovmod as pmp
# import pisacovmods.sequence as pms
# import pisacovmods.contacts as pmc

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov import command_line as pcl
from pisacov import io as pio
from pisacov.io import paths as ppaths
from pisacov.io import conf as pco
from pisacov import sys as psys

from crops.elements import sequences as csq
from crops import io as cio
from crops.io import parsers as cps
from crops.core import ops as cops

from conkit import io as ckio
from conkit import core as ckc
from conkit import plot as ckplot

import argparse
import os
from shutil import copyfile
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
# from string import format

import time

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in crops-renumber."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)
    # MUTUALLY EXCLUSIVE: PDBID or FASTA path + STR path
    main_args = parser.add_mutually_exclusive_group(required=True)
    main_args.add_argument("-i", "--initialise", nargs=2,
                           metavar=("Seqfile", "Structure"),
                           help="Input sequence and structure filepaths.")
    main_args.add_argument("-s", "--skip_conpred", nargs=2,
                           metavar=("Seqfile", "Structure"),
                           help="If HHBLITS and DMP files have already been generated in pdbid/deepmetapsicov, they will be read and those processeses bypassed. Input sequence and structure filepaths.")
    main_args.add_argument("-d", "--skip_default_conpred", nargs=2,
                           metavar=("Seqfile", "Structure"),
                           help="If combined with --add_noncropped, non-default sequences will run through HHBLITS and DMP while default ones will not. If --add_noncropped not used, this function is equivalent to --skip_conpred. Input sequence and structure filepaths.")

    parser.add_argument("-a", "--add_nondefault", action='store_true', default=False,
                        help="Provide results for non-default sequences too. Default sequences are cropped if possible. Non-default sequences are the original sequences if and only if a cropped version exists.")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence. If -s option is used, be aware that this directory must already contain the pdbid/deepmetapsicov directory and its files.")
    parser.add_argument("-c", "--collection_file", nargs=1, metavar="Collection_File",
                        help="Path to CSV file where pisacov signals will be appended. Default: outdir/pisacov_data.csv.")

    parser.add_argument("-u", "--uniprot_threshold", nargs=1, metavar=("Uniprot_ratio_threshold"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. [MIN,MAX)=[0,100).')

    main_args.add_argument("-a", "--hhblits_arguments", nargs=5,
                           metavar=("n_iterations", "evalue_cutoff",
                                    "n_nonreduntant", "mincoverage",
                                    "maxpairwaise_seqidentity"),
                           help=("Introduce HHBLITS arguments." + os.linesep +
                           "DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99]. " +
                           "HHBlits DEEFAULT: [2, 0.001, 1000, 0, 90]"))

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser


def main():
    starttime = time.time()
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    logger.info(pcl.welcome())

    # PARSE CONFIGURATION FILE:
    invals = pco._initialise_inputs()

    invals['INSEQ'] = None
    invals['INSTR'] = None
    invals['ALTDB'] = None
    invals['OUTROOT'] = None
    invals['OUTCSVPATH'] = None
    invals['UPTHRESHOLD'] = None

    # READ INPUT ARGUMENTS
    if args.initialise is not None:
        invals['INSEQ'] = ppaths.check_path(args.initialise[0], 'file')
        invals['INSTR'] = ppaths.check_path(args.initialise[1], 'file')
        skipexec = [False, True] if args.add_nondefault is False else [False, False]
        scoring = [True, False] if args.add_nondefaul is False else [True, True]
    elif args.skip_conpred is not None:
        invals['INSEQ'] = ppaths.check_path(args.skip_conpred[0], 'file')
        invals['INSTR'] = ppaths.check_path(args.skip_conpred[1], 'file')
        skipexec = [True, True]
        scoring = [True, False] if args.add_nondefault is False else [True, True]
    elif args.skip_default_conpred is not None:
        invals['INSEQ'] = ppaths.check_path(args.skip_default_conpred[0], 'file')
        invals['INSTR'] = ppaths.check_path(args.skip_default_conpred[1], 'file')
        skipexec = [True, True] if args.add_nondefault is False else [True, False]
        scoring = [True, False] if args.add_nondefault is False else [True, True]

    if args.outdir is None:
        invals['OUTROOT'] = ppaths.check_path(os.path.dirname(invals['INSEQ']))
    else:
        invals['OUTROOT'] = ppaths.check_path(os.path.join(args.outdir[0], ''))
    ppaths.mdir(invals['OUTROOT'])

    if args.collection_file is None:
        invals['OUTCSVPATH'] = ppaths.check_path(os.path.join(
                                        invals['OUTROOT'], "pisacov_data.csv"))
    else:
        invals['OUTCSVPATH'] = ppaths.check_path(args.collection_file[0])

    if os.path.isfile(invals['OUTCSVPATH']) is False:
        pio.outcsv.csvheader(invals['OUTCSVPATH'])

    if args.uniprot_threshold is not None:
        try:
            invals['UPTHRESHOLD'] = float(args.uniprot_threshold)
        except ValueError:
            logger.critical('Uniprot threshold given not valid.')
        if invals['UNICLUST_FASTA_PATH'] is None:
            invals['UNICLUST_FASTA_PATH'] = pco._uniurl
    else:
        pass

    if args.hhparams is not None:
        invals['HHBLITS_PARAMETERS'] = pco._check_hhparams(args.hhparams)
    else:
        pass

    # Define formats used
    sources = pco._sources()
    n_sources = len(sources)

    # Parse sequence and structure files
    logger.info('Parsing sequence file...')
    seqs = cps.parseseqfile(invals['INSEQ'])

    logger.info('Parsing structure file...')
    strs, filestrs = cps.parsestrfile(invals['INSTR'])

    if len(seqs) == 1 or len(strs) == 1:
        if len(seqs) == 1:
            for key in seqs:
                pdbid = key.lower()
        elif len(seqs) > 1 and len(strs) == 1:
            for key in strs:
                pdbid = key.lower()
    else:
        raise Exception('More than one pdbid in sequence and/or structure set.')

    seq = seqs[pdbid]
    structure = strs[pdbid]

    # logger.info('Parsing SIFTS database file...')
    # sifts = cps.import_db(indb, pdb_in=pdbid)

    # CROPPING AND RENUMBERING
    if skipexec[0] is False:
        logger.info('Cropping and renumbering sequences, ' +
                    'structures according to SIFTS database.')
        psys.crops.runcrops(invals['INSEQ'],
                            invals['INSTR'],
                            invals['SIFTS_PATH'],
                            invals['UPTHRESHOLD'],
                            invals['UNICLUST_FASTA_PATH'],
                            invals['OUTROOT'])

        outpdbdir = os.path.join(invals['OUTROOT'], pdbid, "")
        pio.mdir(outpdbdir)

        inseqc = os.path.join(invals['OUTROOT'],
                              os.path.basename(invals['INSEQ']))
        instrc = os.path.join(invals['OUTROOT'],
                              os.path.basename(invals['INSTR']))

        copyfile(invals['INSEQ'], inseqc)
        copyfile(invals['INSTR'], instrc)

        cmappath = os.path.join(os.path.splitext(inseqc), '.cropmap')

    # MSA GENERATOR # ALL THIS MUST CHANGE
    setch = set()
    cseqpath = os.path.join(outpdbdir, pdbid + '.crops.to_uniprot.fasta')
    hhdir = os.path.join(outpdbdir, 'hhblits','')
    pio.mdir(hhdir)
    neff = {}
    if not skipexec[0] or not skipexec[1]:
        if hhparameters == ['3', '0.001', 'inf', '50', '99']:
            logger.info('Generating Multiple Sequence Alignment using DeepMetaPSICOV default parameters... [AS RECOMMENDED]')
        elif hhparameters == ['2', '0.001', '1000', '0', '90']:
            logger.info('Generating Multiple Sequence Alignment using HHBlits default parameters...')
        else:
            logger.info('Generating Multiple Sequence Alignment using user-custom parameters...')

        if os.path.isfile(cmappath) and not skipexec[0]:
            psys.msagen.runhhblits(cseqpath, hhparameters, hhdir)
            cmsaa3mfile = os.path.splitext(os.path.basename(cseqpath))[0] +".msa.a3m"
            cmsaa3mpath = os.path.join(hhdir, cmsaa3mfile)
            neff['cropped'] = psys.msagen.msafilesgen(cmsaa3mpath)
            neff['original'] = None
            if not skipexec[1]:
                logger.info('    Repeating process for non-default sequence...')
                neff['cropped'] = None
        else:
            logger.info('    No cropped sequence found, using original sequence instead...')

        if not os.path.isfile(cmappath) or not skipexec[1]:
            psys.msagen.runhhblits(inseq, hhparameters, hhdir)
            msaa3mfile = os.path.splitext(os.path.basename(inseq))[0] +".msa.a3m"
            msaa3mpath = os.path.join(hhdir, msaa3mfile)
            neff['original'] = psys.msagen.msafilesgen(msaa3mpath)

    # DEEP META PSICOV RUN
    if not skipexec[0] or not skipexec[1]:
        logger.info('Generating contact prediction lists via DeepMetaPSICOV...')
        dmpdir = os.path.join(outpdbdir, 'dmp','')
        pio.mdir(dmpdir)
        if os.path.isfile(cmappath) and not skipexec[0]:
            psys.dmp.rundmp(cseqpath, cmsaa3mpath, dmpdir)
            if not skipexec[1]:
                logger.info('    Repeating process for non-default sequence...')
        else:
            logger.info('    No cropped sequence found, using original sequence instead...')

        if not os.path.isfile(cmappath) or not skipexec[1]:
            psys.dmp.rundmp(inseq, msaa3mpath, dmpdir)

    # INTERFACE GENERATION, PISA
    cstrpath = os.path.join(outpdbdir, pdbid+'.oldids.crops.to_uniprot.pdb')
    pisadir = os.path.join(outpdbdir, 'pisa','')
    n_ifaces={}
    if not skipexec[0] or not skipexec[1]:
        logger.info('Generating interface files via PISA...')

        if os.path.isfile(cmappath) and not skipexec[0]:
            n_ifaces['cropped']=psys.pisa.runpisa(cstrpath, pisadir)
            if not skipexec[1]:
                logger.info('    Repeating process for non-default sequence...')
        else:
            n_ifaces['cropped']=None
            logger.info('    No cropped sequence found, using original sequence instead...')

        if not os.path.isfile(cmappath) or not skipexec[1]:
            n_ifaces['original']=psys.pisa.runpisa(instr, pisadir)
        else:
            n_ifaces['original']=None

    # CONTACT ANALYSIS AND MATCH
    resultdir = os.path.join(outpdbdir, 'results','')
    logger.info('Opening output csv files...')
    pdbcsvfile=os.path.join(resultdir,pdbid+".evcovsignal.csv")
    pio.outcsv.csvheader(pdbcsvfile)
    if not csvexists:
        pio.outcsv.csvheader(outcsvfile)

    logger.info('Parsing contact predictions lists...')
    ckseq = ckio.read(inseq, 'fasta')
    conpred={}
    for source, attribs in sources.items():
        for mode in ['cropped', 'original']:
            seqfile = cseqpath if mode == 'cropped' else inseq
            confile = (os.path.splitext(os.path.basename(seqfile))[0] +
                       attribs[1])
            conpath = os.path.join(outpdbdir, attribs[0], confile)
            cropmapping = cps.parsemapfile(cmappath)
            if os.path.isfile(conpath):
                conpred[mode][source] = ckio.read(conpath, attribs[2])[0]
                for contact in conpred[mode][source]:
                    contact.res1_seq = cropmapping[pdbid][chid]['cropbackmap'][contact.res1_seq]
                    contact.res2_seq = cropmapping[pdbid][chid]['cropbackmap'][contact.res2_seq]
            else:
                conpred[mode][source] = None

                # NOT SURE IT IS WORKING WITH SEVERAL CHAINS OF SAME SEQUENCE. CHECK EVERYTHING.




    logger.info('    Parsing PISA interfaces...')
    interfaces={}

    for mode in ['cropped', 'original']:
        if n_ifaces[mode] is not None:
            interfaces[mode]=[]
            for i in range(int(n_ifaces[mode])):
                strfile = cstrpath if mode=='cropped' else instr
                pdbfilei = os.path.splitext(strfile)[0]+".interface."+str(i+1)+".pdb"
                interfaces[mode].append(ckio.read(pdbfilei, 'pdb'))
        else:
            interfaces[mode]=None





    # OUTPUT


# CODE TO PRINT NON-REPEATED LINES
#    import csv
#    rows = csv.reader(open("file.csv", "rb"))
#    newrows = []
#    for row in rows:
#        if row not in newrows:
#            newrows.append(row)
#    writer = csv.writer(open("file.csv", "wb"))
#    writer.writerows(newrows)



    return


if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        logger.info(pcl.ok())
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
