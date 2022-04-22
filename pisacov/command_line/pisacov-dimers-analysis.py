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
from pisacov.io import _conf_ops as pco
from pisacov.sys import crops as psc
from pisacov.sys import dmp as psd
from pisacov.sys import msagen as psm
from pisacov.sys import pisa as psp
from pisacov.core import contacts as pcc
from pisacov.core import scores as pcs
from pisacov.core import interfaces as pci

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
                                     description=__description__+' ('+__prog__+')  v.'+__version__ + os.linesep + __doc__,
                                     epilog="Check pisacov.rtfd.io for more information.")
    parser.add_argument("seqpath", nargs=1, metavar=("Seqfile"),
                        help="Input original sequence filepath.")
    parser.add_argument("dimers", nargs='+',  # IF "*" is used and not parsed ok, do it later.
                         metavar=("Dimer_structures"),
                         help="Input the paths to each interface pdb file (this option will skip execution of PISA, and CROPS-related options will be ignored).")
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
                        help="Skip execution of external programs by importing already generated cropped sequence, structure, contacts, etc with default filepaths. This option ignores any values given of HHBLITS parameters or UniProt threshold.")

    # OUTPUT OPTIONS
    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence. If -s option is used, be aware that this directory must already contain the pdbid/deepmetapsicov directory and its files.")
    parser.add_argument("-c", "--collection_file", nargs=1,
                        metavar=("Collection_CSV_Path"),
                        help="Path to CSV file where pisacov signals will be appended. Default: outdir/evcovsignal.full.pisacov.csv.")

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
    invals['INIFS'] = None
    invals['OUTROOT'] = None
    invals['OUTCSVPATH'] = None

    # READ INPUT ARGUMENTS
    invals['INSEQ'] == ppaths.check_path(args.seqpath[0], 'file')


    invals['INIFS'] = []
    args.remove_insertions = False
    for fp in args.dimers:
        if '*' in fp:
            invals['INIFS'] += ppaths.check_wildcard(fp)
        else:
            invals['INIFS'].append(ppaths.check_path(fp, 'file'))
    invals['INIFS'] = list(dict.fromkeys(invals['INIFS']))

    if args.hhblits_arguments is not None:
        invals['HHBLITS_PARAMETERS'] = pco._check_hhparams(args.hhblits_arguments)
    else:
        pass

    if args.skip_conpred is True:
        skipexec = True
        if args.hhblits_arguments is not None:
            logger.info('HHblits parameters given bypassed by --skip_conpred')
    else:
        skipexec = False

    if args.outdir is None:
        invals['OUTROOT'] = ppaths.check_path(os.path.dirname(invals['INSEQ']))
    else:
        invals['OUTROOT'] = ppaths.check_path(os.path.join(args.outdir[0], ''))
    ppaths.mdir(invals['OUTROOT'])

    if args.collection_file is None:
        invals['OUTCSVPATH'] = ppaths.check_path(os.path.join(
                                    invals['OUTROOT'], "evcovsignal.full.pisacov.csv"))
    else:
        invals['OUTCSVPATH'] = ppaths.check_path(args.collection_file[0])

    if os.path.isfile(invals['OUTCSVPATH']) is False:
        pio.outcsv.csvheader(invals['OUTCSVPATH'], cropped=False)

    # Define formats used
    sources = pco._sources()

    # Parse sequence and structure files
    logger.info('Parsing sequence file...')
    seqs = cps.parseseqfile(invals['INSEQ'])

    if len(seqs) == 1:
        if len(seqs) == 1:
            for key in seqs:
                pdbid = key.lower()
    else:
        raise Exception('More than one pdbid in sequence set.')

    seq = seqs[pdbid]

    outpdbdir = os.path.join(invals['OUTROOT'], pdbid, "")

    # RENUMBERING
    fseq = {}
    fmsa = {}

    if skipexec is False:
        if invals['INIFS'] is not None:
            logger.info('Renumbering interfaces provided ' +
                        'according to position in sequence.')
            for path in invals['INIFS']:
                instrc = os.path.join(invals['OUTROOT'],
                                      os.path.basename(path))
                psc.renumcrops(invals['INSEQ'],
                               path,
                               invals['OUTROOT'])
                copyfile(path, instrc)

        pio.mdir(outpdbdir)

    for i, iseq in seq.imer.items():
        fiseq = pdbid + '_' + i + '.fasta'
        fseq[i] = os.path.join(invals['OUTROOT'], pdbid, fiseq)
        fiseq = pdbid + '_' + i + '.msa.aln'
        fmsa[i] = os.path.join(invals['OUTROOT'], pdbid, 'hhblits', fiseq)
        if skipexec is False:
            iseq.dump(fseq[i])

    # EXECUTION OF EXTERNAL PROGRAMS
    hhdir = os.path.join(invals['OUTROOT'], pdbid, 'hhblits', '')
    dmpdir = os.path.join(invals['OUTROOT'], pdbid, 'dmp', '')
    fstr = []
    for file in invals['INIFS']:
        fstr.append(os.path.join(invals['OUTROOT'],
                            (os.path.splitext(os.path.basename(file))[0] +
                             '.crops.seq.pdb')))

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
            sfile = fseq[i]
            afile = fmsa[i]
            themsa = psm.runhhblits(sfile,
                                    invals['HHBLITS_PARAMETERS'],
                                    hhdir)
            iseq.msa = themsa

    # DEEP META PSICOV RUN
        logger.info('Generating contact prediction lists via DeepMetaPSICOV...')

        ppaths.mdir(dmpdir)
        for i, iseq in seq.imer.items():
            sfile = fseq[i]
            afile = fmsa[i]
            psd.rundmp(sfile, afile, dmpdir)

    # GENERATE INTERFACE LIST
    iflist = []
    for filepath in fstr:
        ifname = os.path.splitext(os.path.basename(filepath))[0]
        iflist.append(pci.interface(name = ifname))

    # CONTACT ANALYSIS AND MATCH
    logger.info('Opening output csv files...')
    resultdir = os.path.join(invals['OUTROOT'], pdbid, 'pisacov', '')
    ppaths.mdir(resultdir)
    csvfile = os.path.join(resultdir, pdbid + ".evcovsignal.full.pisacov.csv")

    pio.outcsv.csvheader(csvfile, cropped=False, pisascore=False)

    logger.info('Parsing sequence files...')
    for i, fpath in fseq.items():
        seq.imer[i].seqs['conkit'] = ckio.read(fpath, 'fasta')[0]
        seq.imer[i].biotype = csq.guess_type(seq.imer[i].seqs['mainseq'])

    logger.info('Parsing contact predictions lists...')
    conpred = {}
    matches = []
    for s in seq.imer:
        if s not in conpred:
            conpred[s] = {}
        for source, attribs in sources.items():
            fc = os.path.splitext(os.path.basename(fseq[s]))[0]
            fc += attribs[1]
            confile = os.path.join(invals['OUTROOT'], pdbid, attribs[0], fc)
            conpred[s][source] = ckio.read(confile, attribs[2])[0]

    logger.info('Parsing crystal structure contacts...')
    for i in range(len(iflist)):
        inputmap = ckio.read(fstr[i], 'pdb')
        if len(inputmap.structure) == 4:
            ch = tuple(inputmap.id)
            if seq.whatseq(ch[0]) != seq.whatseq(ch[1]):
                logger.info('Interface ' + str(i) +
                            ' is not a homodimer. Ignoring.')
                iflist[i].structure = None
                matches.append(None)
                continue
            else:
                if seq.imer[seq.whatseq(ch[0])].biotype != "Protein":
                    logger.info('Interface ' + str(i) +
                                ' is not a Protein-Protein interface. Ignoring.')
                    iflist[i].structure = None
                    matches.append(None)
                    continue
                s = seq.whatseq(ch[0])
            try:
                iflist[i].structure = []
                for m in range(len(inputmap)):
                    iflist[i].structure.append(inputmap[m].as_contactmap())
                    iflist[i].structure[m].id = inputmap[m].id
            except Exception:
                for m in range(len(inputmap)):
                    iflist[i].structure.append(inputmap[m])  # ConKit LEGACY.

            matches.append({})
            for source, attribs in sources.items():
                matches[i][source] = pcc.contact_atlas(
                                            name=pdbid+'_'+str(s),
                                            conpredmap=conpred[s][source],
                                            strmap=iflist[i].structure,
                                            sequence=seq.imer[s],
                                            removeintra=True)
        else:
            iflist[i].structure = None
            matches.append(None)
            continue

    logger.info('Computing results and writing them to file...')
    for i in range(len(iflist)):
        if matches[i] is None:
            continue
        results=[pdbid, str(i+1)]
        results.append(matches[i]['psicov'].chain1)
        results.append(matches[i]['psicov'].chain2)
        sid = seq.whatseq(matches[i]['psicov'].chain1)
        results.append(str(sid))
        results.append(str(seq.imer[sid].length()))
        results.append(str(seq.imer[sid].cropmsa.meff))
        results.append(str(seq.imer[sid].ncrops()))
        results.append(str(seq.imer[sid].full_length()))
        for source, attribs in sources.items():
            appresults = pcs.list_scores(matches[i][source], tag=source)
            results += appresults

        pio.outcsv.lineout(results, csvfile)
        pio.outcsv.lineout(results, invals['OUTCSVPATH'])

    logger.info('Analysis finished. Exiting.')

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
