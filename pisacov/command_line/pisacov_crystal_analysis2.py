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
__script__ = 'PISACov Crystal Analysis script'

from pisacov import command_line as pcl

from pisacov.io import paths as ppaths
from pisacov.io import _conf_ops as pco
from pisacov.io import outcsv as pic
from pisacov.sys import crop as psc
from pisacov.sys import dmp as psd
from pisacov.sys import msagen as psm
from pisacov.sys import pisa as psp
from pisacov.core import contacts as pcc
from pisacov.core import scores as pcs
from pisacov.core import interfaces as pci

from crops.io import parsers as cps

import numpy as np
import argparse
import os
import datetime
from shutil import copyfile

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in pisacov_crystal_analysis."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__ + os.linesep + __doc__,
                                     epilog="Check pisacov.rtfd.io for more information.")
    parser.add_argument("seqpath", nargs=1, metavar=("Seqfile"),
                        help="Input original sequence filepath.")
    parser.add_argument("crystalpath", nargs=1,
                         metavar=("Crystal_structure"),
                         help="Input original crystal structure filepath.")

    # Use CROPS
    parser.add_argument("-r", "--remove_insertions", action='store_true', # DO IT FOR SIFTS AND ALSO OPTIONALLY FOR CUSTOM CSV
                             default=False,
                             help="Use CROPS and SIFTS database to remove insertions from crystal structure.")
    parser.add_argument("-u", "--uniprot_threshold", nargs=1,
                            metavar=("Uniprot_ratio_threshold"),
                            help='If Uniprot sequence contribution to crystal sequence is below threshold ( %% ), all its residues are removed. [MIN,MAX)=[0,100).')


    # HHBLITS modification
    parser.add_argument("-a", "--hhblits_arguments", nargs=5,
                        metavar=("n_iterations", "evalue_cutoff",
                                 "n_nonreduntant", "mincoverage",
                                 "maxpairwaise_seqidentity"),
                        help=("Introduce HHBLITS arguments." + os.linesep +
                              "DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99]. " +
                              os.linesep + "HHBlits DEEFAULT: [2, 0.001, 1000, 0, 90]"))

    # SKIP CONPRED
    parser.add_argument("-s", "--skip_conpred", action='store_true',
                        default=False,
                        help="Skip execution of external programs by importing already generated files from default filepath. This option ignores any values given of HHBLITS parameters or UniProt threshold.")

    # OUTPUT OPTIONS
    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence. If -s option is used, be aware that this directory must already contain the pdbid/deepmetapsicov directory and its files.")
    parser.add_argument("-c", "--collection_file", nargs=1,
                        metavar=("Collection_CSV_Path"),
                        help="Path to CSV file where pisacov signals will be appended. Default: outdir/evcovsignal.cropped.pisacov.csv and outdir/evcovsignal.full.pisacov.csv.")

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser


def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    welcomemsg, starttime = pcl.welcome(command=__script__)
    logger.info(welcomemsg)

    # PARSE CONFIGURATION FILE:
    invals = pco._initialise_inputs()

    invals['INSEQ'] = None
    invals['INSTR'] = None
    invals['ALTDB'] = None
    invals['OUTROOT'] = None
    invals['OUTCSVPATH'] = None
    invals['UPTHRESHOLD'] = None

    # READ INPUT ARGUMENTS
    invals['INSEQ'] = ppaths.check_path(args.seqpath[0], 'file')
    invals['INSTR'] = ppaths.check_path(args.crystalpath[0], 'file')

    if args.hhblits_arguments is not None:
        invals['HHBLITS_PARAMETERS'] = pco._check_hhparams(args.hhblits_arguments)
    else:
        pass

    if args.uniprot_threshold is not None:
        try:
            invals['UPTHRESHOLD'] = float(args.uniprot_threshold[0])
        except ValueError:
            logger.critical('Uniprot threshold given not valid.')
        if invals['UNICLUST_FASTA_PATH'] is None:
            invals['UNICLUST_FASTA_PATH'] = pco._uniurl
    else:
        pass

    if args.skip_conpred is True:
        skipexec = True
        if (args.hhblits_arguments is not None or
                args.uniprot_threshold is not None):
            logger.info('HHblits, UniProt threshold parameters given bypassed by --skip_conpred')
    else:
        skipexec = False
    cropping = args.remove_insertions
    scoring = [cropping, not cropping]

    if args.outdir is None:
        invals['OUTROOT'] = ppaths.check_path(os.path.dirname(invals['INSEQ']))
    else:
        invals['OUTROOT'] = ppaths.check_path(os.path.join(args.outdir[0], ''))
    ppaths.mdir(invals['OUTROOT'])

    invals['OUTCSVPATH'] = []
    if args.collection_file is None:
        invals['OUTCSVPATH'].append(ppaths.check_path(os.path.join(
                                    invals['OUTROOT'], ("evcovsignal" +
                                                        os.extsep + "cropped" +
                                                        os.extsep + "pisacov" +
                                                        os.extsep + "csv"))))
        invals['OUTCSVPATH'].append(ppaths.check_path(os.path.join(
                                    invals['OUTROOT'], ("evcovsignal" +
                                                        os.extsep + "full" +
                                                        os.extsep + "pisacov" +
                                                        os.extsep + "csv"))))
    else:
        if cropping is True:
            invals['OUTCSVPATH'].append(ppaths.check_path(args.collection_file[0]))
            invals['OUTCSVPATH'].append(ppaths.check_path(
                os.path.splitext(args.collection_file[0])[0] + os.extsep +
                'full' + os.extsep + os.path.splitext(args.collection_file[0])[1]))
        else:
            invals['OUTCSVPATH'].append(None)
            invals['OUTCSVPATH'].append(ppaths.check_path(args.collection_file[0]))

    # Define formats used
    sources = pco._sources()

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
                for key2 in seqs:
                    if key == key2:
                        pdbid = key.lower()
                    else:
                        if key2.lower() in key.lower():
                            pdbid = key2.lower()
    else:
        raise Exception('More than one pdbid in sequence and/or structure set.')

    seq = seqs[pdbid]
    #structure = strs[pdbid]

    # CROPPING AND RENUMBERING
    outpdbdir = os.path.join(invals['OUTROOT'], pdbid, "")
    instrc = os.path.join(invals['OUTROOT'], pdbid,
                          os.path.basename(invals['INSTR']))

    fseq = {}
    fmsa = {}
    if skipexec is False:
        if cropping is True:
            logger.info('Cropping and renumbering sequences, ' +
                        'structures according to SIFTS database.')
            logger.info(pcl.running('CROPS-cropstr'))
            itime = datetime.datetime.now()
            psc.runcrops(invals['INSEQ'],
                         invals['INSTR'],
                         invals['SIFTS_PATH'],
                         invals['UPTHRESHOLD'],
                         invals['UNICLUST_FASTA_PATH'],
                         invals['OUTROOT'])
            logger.info(pcl.running('CROPS-cropstr', done=itime))
        else:
            logger.info('Renumbering structure ' +
                        'according to position in sequence.')
            logger.info(pcl.running('CROPS-renumber'))
            itime = datetime.datetime.now()
            psc.renumcrops(invals['INSEQ'],
                           invals['INSTR'],
                           invals['OUTROOT'])
            logger.info(pcl.running('CROPS-renumber', done=itime))


        copyfile(invals['INSTR'], instrc)
        ppaths.mdir(outpdbdir)

    for i, iseq in seq.imer.items():
        fiseq = pdbid + '_' + i + '.fasta'
        fseq[i] = os.path.join(invals['OUTROOT'], pdbid, fiseq)
        fiseq = pdbid + '_' + i + '.msa.aln'
        fmsa[i] = os.path.join(invals['OUTROOT'], pdbid, 'hhblits', fiseq)
        if skipexec is False:
            iseq.dump(fseq[i])

    # Parse cropped sequences and maps
    if cropping is True:
        amap = {}
        fcropseq = {}
        fcropmsa = {}
        for i, iseq in seq.imer.items():
            fprefix = pdbid + '_' + i + '.crops.to_uniprot'
            fmap = os.path.join(invals['OUTROOT'], pdbid,
                                fprefix + os.extsep + 'cropmap')
            amap.update(cps.parsemapfile(fmap)[pdbid])
            fcropseq[i] = os.path.join(invals['OUTROOT'], pdbid,
                                       fprefix + os.extsep + 'fasta')
            fcropmsa[i] = os.path.join(invals['OUTROOT'], pdbid, 'hhblits',
                                       (fprefix + os.extsep +
                                        'msa' + os.extsep + 'aln'))
            seq.set_cropmaps(amap, cropmain=True)

    # EXECUTION OF EXTERNAL PROGRAMS
    hhdir = os.path.join(invals['OUTROOT'], pdbid, 'hhblits', '')
    dmpdir = os.path.join(invals['OUTROOT'], pdbid, 'dmp', '')
    pisadir = os.path.join(invals['OUTROOT'], pdbid, 'pisa', '')
    fstr = os.path.join(invals['OUTROOT'], (pdbid + os.extsep + 'crops' +
                                            os.extsep + 'seq' +
                                            os.extsep + 'pdb'))
    if cropping:
        fcropstr = os.path.join(invals['OUTROOT'], pdbid,
                                (pdbid + os.extsep + 'crops' +
                                 os.extsep + 'oldids' +
                                 os.extsep + 'to_uniprot' +
                                 os.path.splitext(invals['INSTR'])[1]))
    if skipexec is False:
        # MSA GENERATOR
        ppaths.mdir(hhdir)
        if invals['HHBLITS_PARAMETERS'] == ['3', '0.001', 'inf', '50', '99']:
            logger.info('Generating Multiple Sequence Alignment using DeepMetaPSICOV default parameters... [AS RECOMMENDED]')
        elif invals['HHBLITS_PARAMETERS'] == ['2', '0.001', '1000', '0', '90']:
            logger.info('Generating Multiple Sequence Alignment using HHBlits default parameters...')
        else:
            logger.info('Generating Multiple Sequence Alignment using user-custom parameters...')

        for i, iseq in seq.imer.items():
            sfile = fcropseq[i] if cropping is True else fseq[i]
            afile = fcropmsa[i] if cropping is True else fmsa[i]
            logger.info(pcl.running('HHBlits'))
            itime = datetime.datetime.now()
            themsa = psm.runhhblits(sfile,
                                    invals['HHBLITS_PARAMETERS'],
                                    hhdir)
            logger.info(pcl.running('HHBlits', done=itime))
            if cropping is True:
                iseq.cropmsa = themsa
                if iseq.ncrops() == 0:
                    iseq.msa = iseq.cropmsa
                    logger.info('    Cropped sequence ' +
                                iseq.oligomer_id + '_' + iseq.name +
                                ' is identical to original sequence.')
                    continue
                else:
                    pass
            else:
                iseq.msa = themsa

    # DEEP META PSICOV RUN
    ppaths.mdir(dmpdir)
    if skipexec is False:
        logger.info('Generating contact prediction lists via DeepMetaPSICOV...')
        for i, iseq in seq.imer.items():
            sfile = fcropseq[i] if cropping is True else fseq[i]
            afile = fcropmsa[i] if cropping is True else fmsa[i]
            nsfile = os.path.join(dmpdir, os.path.basename(sfile))
            if sfile != nsfile:
                copyfile(sfile, nsfile)
            logger.info(pcl.running('DeepMetaPSICOV'))
            itime = datetime.datetime.now()
            psd.rundmp(nsfile, afile, dmpdir)
            logger.info(pcl.running('DeepMetaPSICOV', done=itime))

    # INTERFACE GENERATION, PISA
    ppaths.mdir(pisadir)
    if skipexec is False:
        logger.info('Generating interface files via PISA...')
        sfile = fcropstr if cropping is True else fstr
        logger.info(pcl.running('PISA'))
        itime = datetime.datetime.now()
        iflist = psp.runpisa(sfile, pisadir, sessionid = pdbid)
        logger.info(pcl.running('PISA', done=itime))

    endmsg = pcl.ok(starttime, command=__script__)
    logger.info(endmsg)

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
