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
from pisacov import sys as psys
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
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)
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
                args.uniprot_threshold[0] is not None):
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
                                    invals['OUTROOT'], "evcovsignal.cropped.pisacov.csv")))
        invals['OUTCSVPATH'].append(ppaths.check_path(os.path.join(
                                    invals['OUTROOT'], "evcovsignal.full.pisacov.csv")))
    else:
        if cropping is True:
            invals['OUTCSVPATH'].append(ppaths.check_path(args.collection_file[0]))
            invals['OUTCSVPATH'].append(ppaths.check_path(
                os.path.join(os.path.splitext(args.collection_file[0])[0],
                             'full', os.path.splitext(args.collection_file[0])[1])))
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
    instrc = os.path.join(invals['OUTROOT'],
                          os.path.basename(invals['INSTR']))

    fseq = {}
    fmsa = {}
    if skipexec is False:
        if cropping is True:
            logger.info('Cropping and renumbering sequences, ' +
                        'structures according to SIFTS database.')
            psys.crops.runcrops(invals['INSEQ'],
                                invals['INSTR'],
                                invals['SIFTS_PATH'],
                                invals['UPTHRESHOLD'],
                                invals['UNICLUST_FASTA_PATH'],
                                invals['OUTROOT'])
        else:
            logger.info('Renumbering structure ' +
                        'according to position in sequence.')
            psys.crops.renumcrops(invals['INSEQ'],
                                invals['INSTR'],
                                invals['OUTROOT'])

        copyfile(invals['INSTR'], instrc)
        pio.mdir(outpdbdir)

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
    if cropping:
        fcropstr = os.path.join(invals['OUTROOT'],
                                pdbid + '.crops.oldids.crops.to_uniprot.pdb')
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
            themsa = psys.msagen.runhhblits(sfile,
                                   invals['HHBLITS_PARAMETERS'],
                                   hhdir)
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
            psys.dmp.rundmp(sfile, afile, dmpdir)

    # INTERFACE GENERATION, PISA
    if skipexec is False:
        logger.info('Generating interface files via PISA...')
        sfile = fcropstr if cropping is True else fstr
        iflist = psys.pisa.runpisa(sfile, pisadir)

    # READ DATA IF SKIPEXEC USED:
    if skipexec is True:
        logger.info('Parsing already generated files...')
        for i, iseq in seq.imer.items():
            sfile = fcropstr if cropping is True else fstr
            afile = fcropmsa[i] if cropping is True else fmsa[i]
            if cropping is True:
                iseq.cropmsa = ckio.read(afile, 'jones')
                if iseq.ncrops() == 0:
                    scoring[1] = True
                    iseq.msa = ckio.read(fmsa[i], 'jones')
            else:
                iseq.msa = ckio.read(afile, 'jones')
        ixml = os.path.join(pisadir,
                             os.path.basename(sfile)[0] + '.interface.xml')
        axml = os.path.join(pisadir,
                            os.path.basename(sfile)[0] + '.assembly.xml')

        iflist = pci.parse_interface_xml(ixml, axml)

    # CONTACT ANALYSIS AND MATCH
    logger.info('Opening output csv files...')
    resultdir = os.path.join(invals['OUTROOT'], pdbid, 'pisacov', '')
    ppaths.mdir(resultdir)
    csvfile = []
    csvfile.append(os.path.join(resultdir,
                                pdbid + ".evcovsignal.cropped.pisacov.csv"))
    csvfile.append(os.path.join(resultdir,
                                pdbid + ".evcovsignal.full.pisacov.csv"))

    for n in range(2):
        if scoring[n] is True:
            cpd = True if n == 0 else False
            pio.outcsv.csvheader(csvfile[n], cropped=cpd, pisascore=True)

    logger.info('Parsing sequence files...')
    for i, fpath in fseq.items():
        seq.imer[i].seqs['conkit'] = ckio.read(fpath, 'fasta')[0]

    logger.info('Parsing contact predictions lists...')
    conpred = {}
    matches = []
    for s in seq.imer:
        if s not in conpred:
            conpred[s] = {}
        fs = fcropseq[s] if n == 0 else fseq[s]
        for source, attribs in sources.items():
            fc = os.path.splitext(os.path.basename(fs))[0]
            fc += attribs[1]
            confile = os.path.join(invals['OUTROOT'], pdbid, attribs[0], fc)
            conpred[s][source] = ckio.read(confile, attribs[2])[0]

    logger.info('Parsing crystal structure contacts...')
    for i in range(len(iflist)):
        fs = fcropstr if n == 0 else fstr
        fs = (os.path.splitext(os.path.basename(fs))[0] +
              ".interface." + str(i+1) + ".pdb")
        inputmap = ckio.read(fs, 'pdb')

        if len(inputmap.structure) == 4:
            ch = tuple(inputmap.id)
            if seq.whatseq(ch[0]) != seq.whatseq(ch[1]):
                logger.info('Interface ' + str(i) +
                            ' is not a homodimer. Ignoring.')
                iflist[i].structure = None
                matches.append(None)
                continue
            else:
                chtypes = list(iflist[i].chains.values())
                if chtypes[0] != "Protein" or chtypes[1] != "Protein":
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
        results.append(str(iflist[i].stable))
        for n in range(2):
            if scoring[n] is True:
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
