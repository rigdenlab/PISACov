"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__
__script__ = 'PISACov Crystal Analysis script'

from pisacov import command_line as pcl

from pisacov import iomod as pio
from pisacov.iomod import paths as ppaths
from pisacov.iomod import _conf_ops as pco
from pisacov.iomod import outcsv as pic
from pisacov.sys import crop as psc
from pisacov.sys import dmp as psd
from pisacov.sys import msagen as psm
from pisacov.sys import pisa as psp
from pisacov.core import contacts as pcc
from pisacov.core import scores as pcs
from pisacov.core import interfaces as pci

from crops.iomod import parsers as cps

from conkit import io as ckio

# import numpy as np
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

    #THRESHOLD ABOVE WHICH CONTACTS COUNT
    parser.add_argument("-l", "--lower_threshold", nargs=3,
                        metavar=("psicov", "ccmpred", "dmp"),
                        help=("Remove predicted contacts scored below these values. Use '-inf' for no lower cutoff."))

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
    parser.add_argument("-p", "--plot_formats", nargs='+',
                        metavar=("Plot file format(s)"),
                        help="One or more formats of 'png', 'eps' and 'agr' of figures/data to be produced.")

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

    if args.lower_threshold is not None:
        invals['LOWTHRESHOLD'] = pco._check_lowth(args.lower_threshold)
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

    if args.plot_formats is None:
        plotformats={'png'}
    else:
        plotformats=set()
        for element in args.plot_formats:
            if element.lower() in {'png', 'eps', 'agr'}:
                plotformats.add(element.lower())

    # Define formats used
    sources = pco._sources(lowth=invals['LOWTHRESHOLD'])

    # Parse sequence and structure files
    logger.info('Parsing sequence file...')
    # seqs = cps.parseseqfile(invals['INSEQ'])
    seqs = pio.read(invals['INSEQ'], 'fasta')

    logger.info('Parsing structure file...')
    # strs, filestrs = cps.parsestrfile(invals['INSTR'])
    strs, filestrs = pio.read(invals['INSTR'], 'pdb')

    if len(seqs) == 1 or len(strs) == 1:
        if len(seqs) == 1:
            for key in seqs:
                pdbid = key
        elif len(seqs) > 1 and len(strs) == 1:
            for key in strs:
                for key2 in seqs:
                    if key.upper() == key2.upper():
                        pdbid = key.upper()
                    else:
                        if key2.upper() in key.upper():
                            pdbid = key2.upper()
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

        ppaths.mdir(outpdbdir)
        copyfile(invals['INSTR'], instrc)

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
            if iseq.ncrops() == 0:
                logger.info('    Cropped sequence ' +
                            iseq.oligomer_id + '_' + iseq.name +
                            ' is identical to the original sequence.')
            else:
                logger.info('    Cropped sequence ' +
                            iseq.oligomer_id + '_' + iseq.name +
                            ' is ' + str(iseq.ncrops()) + ' residues ' +
                            'shorter than the original sequence.')

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
            logger.info('Generating Multiple Sequence Alignment using DeepMetaPSICOV default parameters... [RECOMMENDED OPTION]')
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

    # READ DATA IF SKIPEXEC USED
    if skipexec is True:
        logger.info('Parsing already generated files...')
        for i, iseq in seq.imer.items():
            sfile = fcropstr if cropping is True else fstr
            afile = fcropmsa[i] if cropping is True else fmsa[i]
            if cropping is True:
                # iseq.cropmsa = ckio.read(afile, 'jones')
                iseq.cropmsa = pio.read(afile, 'jones')
                if iseq.ncrops() == 0:
                    scoring[1] = True
                    # iseq.msa = ckio.read(afile, 'jones')
                    iseq.msa = ckio.read(afile, 'jones')
            else:
                # iseq.msa = ckio.read(afile, 'jones')
                iseq.msa = pio.read(afile, 'jones')
        ixml = os.path.join(pisadir,
                            (os.path.splitext(os.path.basename(sfile))[0] +
                             os.extsep + 'interface' +
                             os.extsep + 'xml'))
        axml = os.path.join(pisadir,
                            (os.path.splitext(os.path.basename(sfile))[0] +
                             os.extsep + 'assembly' +
                             os.extsep + 'xml'))

        iflist = pci.parse_interface_xml(ixml, axml)

    # CONTACT ANALYSIS AND MATCH
    logger.info('Opening output csv files...')
    resultdir = os.path.join(invals['OUTROOT'], pdbid, 'pisacov', '')
    ppaths.mdir(resultdir)
    csvfile = []
    csvfile.append(os.path.join(resultdir,
                                (pdbid + os.extsep + "evcovsignal" +
                                 os.extsep + "cropped" +
                                 os.extsep + "pisacov" +
                                 os.extsep + "csv")))
    csvfile.append(os.path.join(resultdir,
                                (pdbid + os.extsep + "evcovsignal" +
                                 os.extsep + "full" +
                                 os.extsep + "pisacov" +
                                 os.extsep + "csv")))

    for n in range(2):
        if scoring[n] is True:
            cpd = True if cropping else False
            pic.csvheader(csvfile[n], cropped=cpd, pisascore=True)
            if invals['OUTCSVPATH'][n] is not None:
                if os.path.isfile(invals['OUTCSVPATH'][n]) is False:
                    pic.csvheader(invals['OUTCSVPATH'][n], cropped=cpd, pisascore=True)

    logger.info('Parsing sequence files...')
    for i, fpath in fseq.items():
        # seq.imer[i].seqs['conkit'] = ckio.read(fpath, 'fasta')[0]
        seq.imer[i].seqs['conkit'] = pio.read(fpath, 'fasta', ck=True)[0]

    logger.info('Parsing contact predictions lists...')
    conpred = {}
    matches = []
    for s in seq.imer:
        if s not in conpred:
            conpred[s] = {}
        fs = fcropseq[s] if cropping else fseq[s]
        for source, attribs in sources.items():
            fc = os.path.splitext(os.path.basename(fs))[0]
            fc += os.extsep + attribs[1]
            confile = os.path.join(dmpdir, fc)
            # conpred[s][source] = ckio.read(confile, attribs[2])[0]
            conpred[s][source] = pio.read(confile, attribs[2], ck=True)[0]

    logger.info('Parsing crystal structure contacts...')
    for i in range(len(iflist)):
        logger.info(os.linesep + str(iflist[i]))
        fs = fcropstr if cropping else fstr
        fs = (os.path.splitext(os.path.basename(fs))[0] +
              os.extsep + "interface" + os.extsep + str(i+1) + os.extsep + "pdb")
        spath = os.path.join(pisadir, fs)
        # inputmap = ckio.read(spath, 'pdb')
        inputmap = pio.read(spath, 'pdb', ck=True)
        if len(inputmap) == 4:
            chnames = [iflist[i].chains[0].crystal_id,
                       iflist[i].chains[1].crystal_id]
            iflist[i].chains[0].seq_id = seq.whatseq(chnames[0])
            iflist[i].chains[1].seq_id = seq.whatseq(chnames[1])
            chseqs = [iflist[i].chains[0].seq_id,
                       iflist[i].chains[1].seq_id]

            logger.info(iflist[i].chains)
            chtypes = [iflist[i].chains[0].type,
                       iflist[i].chains[1].type]
            if (chseqs[0] != chseqs[1] or
                    (chtypes[0] != 'Protein' or chtypes[1] != 'Protein')):
                if chtypes[0] != "Protein" or chtypes[1] != "Protein":
                    logger.info('Interface ' + str(i) +
                                ' is not a Protein-Protein interface. Ignoring.')
                else:
                    logger.info('Interface ' + str(i) +
                                ' is not a homodimer. Ignoring.')
                iflist[i].structure = None
                matches.append(None)
                continue
            s = chseqs[0]

            try:
                iflist[i].structure = []
                for m in range(len(inputmap)):
                    iflist[i].structure.append(inputmap[m].as_contactmap())
                    iflist[i].structure[m].id = inputmap[m].id
            except Exception:
                logger.warning('Contact Maps obtained from a legacy ConKit ' +
                               'version with no Distograms implemented.')
                for m in range(len(inputmap)):
                    iflist[i].structure.append(inputmap[m])  # ConKit LEGACY.
            #fs = fcropstr if cropping else fstr
            #fs = (os.path.splitext(os.path.basename(fs))[0] +
            #      os.extsep + "interface" + os.extsep + str(i+1) + os.extsep + "con")
            #spath = os.path.join(pisadir, fs)
            #pio.write(spath, 'psicov', indata=iflist[i].structure[1])
            #iflist[i].contactmap = pio.read(spath, 'array')
            iflist[i].contactmap = iflist[i].structure[1].deepcopy()
            matches.append({})
            for source, attribs in sources.items():
                matches[i][source] = pcc.contact_atlas(
                                            name=pdbid+'_'+str(s),
                                            dimer_interface=iflist[i],
                                            conpredmap=conpred[s][source],
                                            conpredtype=source,
                                            sequence=seq.imer[s])
                if cropping is True:
                    matches[i][source].set_cropmap()
                matches[i][source].remove_neighbours(mindist=2)
                matches[i][source].set_conpred_seq()
                matches[i][source].remove_intra()
                matches[i][source].make_match(filterout = attribs[3])
                for cmode, cmap in matches[i][source].conkitmatch.items():
                    if (len(cmap) > 0 and
                            len(matches[i][source].interface.structure[1]) > 0):
                        for imtype in plotformats:
                            if len(matches[i][source].conkitmatch) > 1:
                                pout = (os.path.splitext(fs)[0] + os.extsep + 'match' +
                                        os.extsep + cmode + os.extsep + source +
                                        os.extsep + 'con' + os.extsep + imtype)
                            else:
                                pout = (os.path.splitext(fs)[0] + os.extsep +
                                        'match' + os.extsep + source +
                                        os.extsep + 'con' + os.extsep + imtype)
                            plotpath = os.path.join(os.path.dirname(csvfile[0]), pout)
                            matches[i][source].plot_map_alt(plotpath,
                                                            mode = cmode,
                                                            plot_type = imtype)
        else:
            iflist[i].structure = None
            iflist[i].contactmap = None
            matches.append(None)
            continue

    logger.info(os.linesep + 'Computing results and writing them to file...' +
                os.linesep)
    for i in range(len(iflist)):
        logger.info('Generating Interface ' + str(i+1) + ' data...')
        if matches[i] is None:
            continue
        results = [pdbid, str(i+1)]
        results.append(iflist[i].chains[0].crystal_id)
        results.append(iflist[i].chains[1].crystal_id)
        sid = iflist[i].chains[0].seq_id
        results.append(str(sid))
        results.append(str(seq.imer[sid].length()))
        if cropping is True:
            results.append(str(seq.imer[sid].cropmsa.meff))
        else:
            results.append(str(seq.imer[sid].msa.meff))
        results.append(str(seq.imer[sid].ncrops()))
        results.append(str(seq.imer[sid].full_length()))
        results.append(str(seq.imer[sid].msa.meff))
        for source, attribs in sources.items():
            appresults = pcs.list_scores(matches[i][source], tag=source)
            results.extend(appresults)
        results.append(str(iflist[i].stable))
        for n in range(2):
            if scoring[n] is True:
                pic.lineout(results, csvfile[n])
                pic.lineout(results, invals['OUTCSVPATH'][n])

    endmsg = pcl.ok(starttime, command=__script__)
    logger.info(endmsg)

    return

if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
