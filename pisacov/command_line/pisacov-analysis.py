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
from pisacov.core import contacts as pcc
from pisacov.core import scores as pcs

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

    # MUTUALLY EXCLUSIVE: Crystal PDB structure or PDB dimers from any source
    structure_args = parser.add_mutually_exclusive_group(required=True)
    structure_args.add_argument("-c", "--crystal_structure", nargs=1,
                                metavar=("Crystal_structure"),
                                help="Input original crystal structure filepath.")
    structure_args.add_argument("-d", "--dimers", nargs='+',  # IF "*" is used and not parsed ok, do it later.
                                metavar=("Dimer_structures"),
                                help="Input the paths to each interface pdb file (this option will skip execution of PISA, and CROPS-related options will be ignored).")

    # Use CROPS
    parser.add_argument("-r", "--remove_insertions", action='store_true', # DO IT FOR SIFTS AND ALSO OPTIONALLY FOR CUSTOM CSV
                             default=False,
                             help="Use CROPS and SIFTS database to remove insertions from crystal structure. ONLY if crystal structure provided, otherwise this option is ignored.")
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
                        help="Skip execution of external programs by importing already generated cropped sequence, structure, contacts, etc with default filepaths. This option ignores any values given of HHBLITS parameters or UniProt threshold.")

    # OUTPUT OPTIONS
    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence. If -s option is used, be aware that this directory must already contain the pdbid/deepmetapsicov directory and its files.")
    parser.add_argument("-C", "--collection_file", nargs=1,
                        metavar=("Collection_CSV_Path"),
                        help="Path to CSV file where pisacov signals will be appended. Default: outdir/evcovsignal.cropped.pisacov.csv and outdir/evcovsignal.full.pisacov.csv.")

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
    invals['INIFS'] = None
    invals['ALTDB'] = None
    invals['OUTROOT'] = None
    invals['OUTCSVPATH'] = None
    invals['UPTHRESHOLD'] = None

    # READ INPUT ARGUMENTS
    invals['INSEQ'] == ppaths.check_path(args.seqpath, 'file')

    if args.crystal_structure is not None:
        invals['INSTR'] = ppaths.check_path(args.crystal_structure, 'file')
    elif args.dimers is not None:
        invals['INIFS'] = []
        args.remove_insertions = False
        for fp in args.dimers:
            if '*' in fp:
                invals['INIFS'] += ppaths.check_wildcard(fp)
            else:
                invals['INIFS'].append(ppaths.check_path(fp, 'file'))
        invals['INIFS'] = list(dict.fromkeys(invals['INIFS']))
        if args.uniprot_threshold is not None:
            logger.info('Uniprot threshold given bypassed by --dimers')
        if args.remove_insertions is not None:
            logger.info('Insertion removal bypassed by --dimers')

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

    if args.skip_conpred is True:
        skipexec = [True, True]
        if args.uniprot_threshold is not None:
            logger.info('Uniprot threshold given bypassed by --skip_conpred')
        if args.hhblits_arguments is not None:
            logger.info('HHblits parameters given bypassed by --skip_conpred')
    else:
        skipexec = [not args.remove_insertions, args.remove_insertions]
    scoring = [args.remove_insertions, not args.remove_insertions]

    if args.outdir is None:
        invals['OUTROOT'] = ppaths.check_path(os.path.dirname(invals['INSEQ']))
    else:
        invals['OUTROOT'] = ppaths.check_path(os.path.join(args.outdir[0], ''))
    ppaths.mdir(invals['OUTROOT'])

    if args.collection_files is None:
        invals['OUTCSVPATH'] = []
        invals['OUTCSVPATH'].append(ppaths.check_path(os.path.join(
                                    invals['OUTROOT'], "evcovsignal.cropped.pisacov.csv")))
        invals['OUTCSVPATH'].append(ppaths.check_path(os.path.join(
                                    invals['OUTROOT'], "evcovsignal.full.pisacov.csv")))
    else:
        invals['OUTCSVPATH'] = ppaths.check_path(args.collection_files)

    if os.path.isfile(invals['OUTCSVPATH'][0]) is False:
        pio.outcsv.csvheader(invals['OUTCSVPATH'], cropped=True)
    if os.path.isfile(invals['OUTCSVPATH'][1]) is False:
        pio.outcsv.csvheader(invals['OUTCSVPATH'], cropped=False)

    # Define formats used
    sources = pco._sources()
    n_sources = len(sources)

    # Parse sequence and structure files
    logger.info('Parsing sequence file...')
    seqs = cps.parseseqfile(invals['INSEQ'])

    if invals['INSTR'] is not None:
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
    else:
        if len(seqs) == 1:
            if len(seqs) == 1:
                for key in seqs:
                    pdbid = key.lower()
        else:
            raise Exception('More than one pdbid in sequence set.')
        strs = None
    seq = seqs[pdbid]

    outpdbdir = os.path.join(invals['OUTROOT'], pdbid, "")

    # CROPPING AND RENUMBERING
    if invals['INSTR'] is not None:
        instrc = os.path.join(invals['OUTROOT'],
                              os.path.basename(invals['INSTR']))
    fseq = {}
    fmsa = {}

    if skipexec[1] is False:
        if invals['INIFS'] is not None:
            logger.info('Renumbering interfaces provided ' +
                        'according to position in sequence.')
            for path in invals['INIFS']:
                instrc = os.path.join(invals['OUTROOT'],
                                      os.path.basename(path))
                psys.crops.runcrops(invals['INSEQ'],
                                    path,
                                    invals['OUTROOT'])
                copyfile(path, instrc)
        else:
            logger.info('Renumbering structures ' +
                        'according to position in sequence.')
            psys.crops.renumcrops(invals['INSEQ'],
                                invals['INSTR'],
                                invals['OUTROOT'])
            copyfile(invals['INSTR'], instrc)

        pio.mdir(outpdbdir)

    if skipexec[0] is False:
        logger.info('Cropping and renumbering sequences, ' +
                    'structures according to SIFTS database.')
        psys.crops.runcrops(invals['INSEQ'],
                            invals['INSTR'],
                            invals['SIFTS_PATH'],
                            invals['UPTHRESHOLD'],
                            invals['UNICLUST_FASTA_PATH'],
                            invals['OUTROOT'])

        pio.mdir(outpdbdir)
        copyfile(invals['INSTR'], instrc)

    for i, iseq in seq.imer.items():
        fiseq = pdbid + '_' + i + '.fasta'
        fseq[i] = os.path.join(invals['OUTROOT'], pdbid, fiseq)
        fiseq = pdbid + '_' + i + '.msa.aln'
        fmsa[i] = os.path.join(invals['OUTROOT'], pdbid, 'hhblits', fiseq)
        if skipexec[0] is False or skipexec[1] is False:
            iseq.dump(fseq[i])

    # Parse cropped sequences and maps
    amap = {}
    fcropseq = {}
    fcropmsa = {}
    for i, iseq in seq.imer.items():
        fprefix = pdbid + '_' + i + '.crops.to_uniprot'
        fmap = os.path.join(invals['OUTROOT'], pdbid, fprefix + '.cropmap')
        amap.update(cps.parsemapfile(fmap)[pdbid])
        fcropseq[i] = os.path.join(invals['OUTROOT'], pdbid, fprefix + '.fasta')
        fcropmsa[i] = os.path.join(invals['OUTROOT'], pdbid,
                                   'hhblits', fmap + '.msa.aln')
    seq.set_cropmaps(amap, cropmain=True)

    # EXECUTION OF EXTERNAL PROGRAMS
    hhdir = os.path.join(invals['OUTROOT'], pdbid, 'hhblits', '')
    dmpdir = os.path.join(invals['OUTROOT'], pdbid, 'dmp', '')
    pisadir = os.path.join(invals['OUTROOT'], pdbid, 'pisa', '')
    fstr = os.path.join(invals['OUTROOT'], pdbid + '.crops.seq.pdb')
    fcropstr = os.path.join(invals['OUTROOT'],
                            pdbid + '.crops.oldids.crops.to_uniprot.pdb')
    n_ifaces=[]
    if skipexec[0] is False or skipexec[1] is False:
        # MSA GENERATOR
        ppaths.mdir(hhdir)

        if invals['HHBLITS_PARAMETERS'] == ['3', '0.001', 'inf', '50', '99']:
            logger.info('Generating Multiple Sequence Alignment using DeepMetaPSICOV default parameters... [AS RECOMMENDED]')
        elif invals['HHBLITS_PARAMETERS'] == ['2', '0.001', '1000', '0', '90']:
            logger.info('Generating Multiple Sequence Alignment using HHBlits default parameters...')
        else:
            logger.info('Generating Multiple Sequence Alignment using user-custom parameters...')

        for i, iseq in seq.imer.items():
            for n in range(2):
                if skipexec[n] is False:
                    sfile = fcropseq[i] if n == 0 else fseq[i]
                    afile = fcropmsa[i] if n == 0 else fmsa[i]
                    themsa = psys.msagen.runhhblits(sfile,
                                           invals['HHBLITS_PARAMETERS'],
                                           hhdir)
                    if n == 0:
                        iseq.cropmsa = themsa
                    else:
                        iseq.msa = themsa
                    if skipexec[1] is False and n == 0:
                        if iseq.ncrops() == 0:
                            iseq.msa = iseq.cropmsa
                            logger.info('    No cropping was performed for ' +
                                        iseq.oligomer_id + '_' + iseq.name +
                                        ', only original sequence considered.')
                            break
                        else:
                            logger.info('    Repeating process for non-default sequence...')

    # DEEP META PSICOV RUN
        logger.info('Generating contact prediction lists via DeepMetaPSICOV...')

        ppaths.mdir(dmpdir)
        for i, iseq in seq.imer.items():
            for n in range(2):
                if skipexec[n] is False:
                    sfile = fcropseq[i] if n == 0 else fseq[i]
                    afile = fcropmsa[i] if n == 0 else fmsa[i]
                    psys.dmp.rundmp(sfile, afile, dmpdir)
                    if skipexec[1] is False and n == 0:
                        if iseq.ncrops() == 0:
                            logger.info('    No contact prediction was performed for ' +
                                        iseq.oligomer_id + '_' + iseq.name +
                                        ', only original sequence considered.')
                            break
                        else:
                            logger.info('    Repeating process for non-default sequence...')

    # INTERFACE GENERATION, PISA
        logger.info('Generating interface files via PISA...')
        for n in range(2):
            if skipexec[n] is False:
                if n == 0:
                    n_ifaces.append(psys.pisa.runpisa(fcropstr, pisadir))
                else:
                    n_ifaces.append(psys.pisa.runpisa(fstr, pisadir))
                if skipexec[1] is False and n == 0:
                    if iseq.ncrops() == 0:
                        if iseq.ncrops() == 0:
                            logger.info('    No PISA analysis was performed for ' +
                                        iseq.oligomer_id +
                                        ', only original structure considered.')
                            break
                        else:
                            logger.info('    Repeating process for non-default structure...')

    # READ DATA IF SKIPEXEC USED:
    if skipexec[0] is True:
        logger.info('Parsing already generated files...')
        for i, iseq in seq.imer.items():
            iseq.cropmsa = ckio.read(fcropmsa[i], 'jones')
        xfile = os.path.join(pisadir,
                             os.path.basename(fcropstr)[0] + '.interface.xml')
        n_ifaces.append(psys.pisa.n_int_xml(xfile))

    if skipexec[1] is True and scoring[1] is True:
        for i, iseq in seq.imer.items():
            if iseq.ncrops() == 0:
                iseq.msa = iseq.cropmsa # pass?
            else:
                iseq.msa = ckio.read(fmsa[i], 'jones')
        xfile = os.path.join(pisadir,
                             os.path.basename(fstr)[0] + '.interface.xml')
        n_ifaces.append(psys.pisa.n_int_xml(xfile))


    # CONTACT ANALYSIS AND MATCH
    logger.info('Opening output csv files...')
    resultdir = os.path.join(invals['OUTROOT'], pdbid, 'pisacov', '')
    ppaths.mdir(resultdir)
    csvfile=[]
    csvfile.append(os.path.join(resultdir,
                                pdbid + ".evcovsignal.full.pisacov.csv"))
    csvfile.append(os.path.join(resultdir,
                                pdbid + ".evcovsignal.cropped.pisacov.csv"))

    for n in range(2):
        if scoring[n] is True:
            pio.outcsv.csvheader(csvfile[n])

    logger.info('Parsing sequence files...')
    for i, fpath in fseq.items():
        seq.imer[i].seqs['conkit'] = ckio.read(fpath, 'fasta')[0]

    logger.info('Parsing contact predictions lists...')
    conpred = [{}, {}]
    pdbmaps = [[], []]
    matches = [[], []]
    for n in range(2):
        if n == 0 or scoring[1] is True:

            for s in seq.imer:
                if s not in conpred[n]:
                    conpred[n][s] = {}
                fs = fcropseq[s] if n == 0 else fseq[s]
                for source, attribs in sources.items():
                    fc = os.path.splitext(os.path.basename(fs))[0]
                    fc += attribs[1]
                    confile = os.path.join(invals['OUTROOT'], pdbid,
                                            attribs[0], fc)
                    conpred[n][s][source] = ckio.read(confile, attribs[2])[0]

            for i in range(n_ifaces[n]):
                matches[n].append([])
                fs = fcropstr if n == 0 else fstr
                fs = (os.path.splitext(os.path.basename(fs))[0] +
                           ".interface."+str(i+1)+".pdb")
                inputmap = ckio.read(fs, 'pdb')
                if len(inputmap) == 4:
                    ch = tuple(inputmap.id)
                    if seq.whatseq(ch[0]) != seq.whatseq(ch[1]):
                        logger.info('Interface ' + str(i) +
                                    ' is not a homodimer. Ignoring.')
                        continue
                    else:
                        s = seq.whatseq(ch[0])
                    try:
                        pdbmaps[n].append([])
                        for m in range(len(inputmap)):
                            pdbmaps[n][i].append(inputmap[m].as_contactmap())
                            pdbmaps[n][i][m].id = inputmap[m].id
                    except Exception:
                        for m in range(len(inputmap)):
                            pdbmaps[n][i].append(inputmap[m])  # ConKit LEGACY.

                    matches[n][i].append({})
                    for source, attribs in sources.items():
                        matches[n][i][source] = pcc.contact_atlas(
                                                    name=pdbid+'_'+str(s),
                                                    conpredmap=conpred[n][s][source],
                                                    strmap=pdbmaps[n][i],
                                                    sequence=seq.imer[s],
                                                    removeintra=True)

    logger.info('Computing results and writing them to file...')
    for n in range(2):
        for i in range(n_ifaces[n]):
            results=[pdbid, str(i+1)]
            results.append(matches[n][i]['psicov'].chain1)
            results.append(matches[n][i]['psicov'].chain2)
            sid = seq.whatseq(matches[n][i]['psicov'].chain1)
            results.append(str(sid))
            results.append(str(seq.imer[sid].length()))
            results.append(str(seq.imer[sid].cropmsa.meff))
            results.append(str(seq.imer[sid].ncrops()))
            results.append(str(seq.imer[sid].full_length()))
            for source, attribs in sources.items():
                appresults = pcs.list_scores(matches[n][i][source], tag=source)
                results += appresults

            pio.outcsv.lineout(results, csvfile[n])
            pio.outcsv.lineout(results, invals['OUTCSVPATH'][n])


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
