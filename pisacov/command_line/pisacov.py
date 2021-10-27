"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

#import pisacovmods.inputvalues as pmin
#import pisacovmods.init1 as pmi1
#import pisacovmods.init2 as pmi2
#import pisacovmods.output as pmo
#import pisacovmods.pisacovmod as pmp
#import pisacovmods.sequence as pms
#import pisacovmods.contacts as pmc

from pisacov.about import __prog__, __description__, __version__
from pisacov.about import  __author__, __date__, __copyright__

from pisacov import command_line as pcl
from pisacov import io as pio
from pisacov import sys as psys

from crops.elements import sequence as csq
from crops import io as cio
from crops.io import parsers as cps
from crops.core import ops as cops

from conkit import io as ckio
from conkit import core as ckc
from conkit import plot as ckplot

import argparse
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
    # MUTUALLY EXCLUSIVE: PDBID or FASTA path + STR path
    main_args = parser.add_mutually_exclusive_group(required=True)
    main_args.add_argument("-i", "--initialise", nargs=2, metavar="Initial_files",
                           help="Input sequence and structure filepaths.")
    main_args.add_argument("-s", "--skip_conpred", nargs=2, metavar="Initial_files",
                           help="If HHBLITS and DMP files have already been generated in pdbid/deepmetapsicov, they will be read and those processeses bypassed. Input sequence and structure filepaths.")

    parser.add_argument("-o","--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence. If -s option is used, be aware that this directory must already contain the pdbid/deepmetapsicov directory and its files.")
    parser.add_argument("-c","--collection_file", nargs=1, metavar="Collection_File",
                        help="Path to CSV file where pisacov signals will be appended. Default: oudir/pisacov_data.csv.")

    parser.add_argument("-u","--uniprot_threshold", nargs=1, metavar=("Uniprot_ratio_threshold"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. [MIN,MAX)=[0,100).')

    parser.add_argument("-h","--hhparams", nargs=5, metavar=("hhblits_new_parameters"),
                          help='Override default HHBLITS parameters in config file by introducing new ones: #iterations, E-value cutoff, Non-redundant seqs to keep, MinimumCoverageWithMasterSeq(%), MaxPairwiseSequenceIdentity.')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser

def main():
    starttime=time.time()
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    logger.info(pcl.welcome())

    # READ INPUT ARGUMENTS
    if args.initialise is not None:
        inseq = pio.check_path(args.initialise[0], 'file')
        instr = pio.check_path(args.initialise[1], 'file')
        indb = pio.check_path(pio.conf.CSV_CHAIN_PATH, 'file')
        skipexec = False
    elif args.skip_conpred is not None:
        inseq = pio.check_path(args.skip_conpred[0], 'file')
        instr = pio.check_path(args.skip_conpred[1], 'file')
        skipexec = True

    if args.outdir is None:
        outrootdir = pio.check_path(os.path.dirname(inseq))
    else:
        outrootdir = pio.check_path(os.path.join(args.outdir[0], ''))

    if args.collection_file is None:
        outcsv = pio.check_path(os.path.join(outrootdir,"pisacov_data.csv"))
    else:
        outcsv = pio.check_path(args.collection_file[0])
    try:
        pio.check_path(outcsv, 'file')
        csvexists = True
    except:
        csvexists = False

    if args.uniprot_threshold is not None:
        thuprot, dbuprot = pio.check_uniprot(args.uniprot_threshold[0])
    else:
        thuprot = 0.0
        dbuprot = None

    if args.hhparams is not None:
        hhparameters = pio.check_hhparams(args.hhparams)
    else:
        try:
            hhparameters = pio.check_hhparams(pio.conf.HHBLITS_PARAMETERS)
        except:
            hhparameters = pio.check_hhparams('dmp')

    # Define formats used
    sources = pio.paths.sources()
    n_sources = len(sources)

    # Parse sequence and structure files
    logger.info('Parsing sequence file...')
    seq = cps.parseseqfile(inseq)
    if len(seq)==1:
        for key in seq.keys():
            pdbid=key.lower()
    else:
        raise Exception('More than one pdbid in sequence set.')

    logger.info('Parsing structure file...')
    structure = cps.parsestrfile(instr)[0][pdbid]

    logger.info('Parsing SIFTS database file...')
    sifts = cps.import_db(indb, pdb_in=pdbid)

    # CROPPING AND RENUMBERING
    if not skipexec:
        logger.info('Cropping and renumbering sequences, structures according to SIFTS database.')
        psys.crops(inseq, instr, indb, thuprot, dbuprot, outrootdir)



    # MSA GENERATOR
    if not skipexec:
        if hhparameters == 'dmp':
            logger.info('Generating Multiple Sequence Alignment using DeepMetaPSICOV default parameters... [AS RECOMMENDED]')
        elif hhparameters == [2, 0.001, 1000, 0, 90]:
            logger.info('Generating Multiple Sequence Alignment using HHBlits default parameters...')
        else:
            logger.info('Generating Multiple Sequence Alignment using user-custom parameters...')

    path=os.path.join(os.getcwd(),pdbid,pdbid+'.crops.oldids.to_uniprot.a3m')
    seqalign=ckio.read(os.path.join(os.getcwd(),pdbid,pdbid+'.crops.oldids.to_uniprot.a3m'),'a3m')
    neff=seqalign.meff

    if not skipexec:
        fig = ckplot.SequenceCoverageFigure(seqalign)
        fig.savefig(path+".png", overwrite=True)

    # DEEP META PSICOV RUN

    try:
        inxml=cio.check_path(args.input_interfaces[0],'file')
        xml=ET.parse(inxml)
    except:
        raise argparse.ArgumentError()

    # CONTACT ANALYSIS AND MATCH



    # OUTPUT


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