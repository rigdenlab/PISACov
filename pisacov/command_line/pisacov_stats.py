"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__
__script__ = 'PISACov Statistical Analysis script'

from pisacov import command_line as pcl
from pisacov.iomod import paths as ppaths
from pisacov.iomod import _conf_ops as pco
from pisacov.iomod import outcsv as pic
from pisacov.iomod import plots as pip
from pisacov.core import stats as pcs

import argparse
import os
import csv
import numpy as np
import copy

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in pisacov_stats."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__ + os.linesep + __doc__,
                                     epilog="Check pisacov.rtfd.io for more information.")

    parser.add_argument('scores', nargs=1, metavar=("ScoresCSVFile"),
                        help="Input scores CSV filepath.")
    parser.add_argument("-f", "--full_score_analysis", action='store_true',
                        default=False,
                        help="Produce full analysis of beta score list (beta score list required).")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input data.")
    parser.add_argument("-p", "--plot_formats", nargs='+',
                        metavar=("Plot file format(s)"),
                        help="One or more formats of 'png' and 'eps' of figures to be produced.")

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser


def main():

    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    welcomemsg, starttime = pcl.welcome(command=__script__)
    logger.info(welcomemsg)

    csvfile = ppaths.check_path(args.scores[0], 'file')
    if args.outdir is None:
        outdir = ppaths.check_path(os.path.join(os.path.dirname(csvfile), 'stats', ''))
    else:
        if os.path.basename(os.path.dirname(os.path.abspath(args.outdir[0]))) == 'stats':
            outdir = ppaths.check_path(os.path.join(args.outdir[0], ''))
        else:
            outdir = ppaths.check_path(os.path.join(args.outdir[0], 'stats', ''))
        ppaths.mdir(outdir)

    fcurves = os.path.splitext(os.path.basename(csvfile))[0] + ".TPRvFPR.rocs.csv"
    fcurves = os.path.join(outdir, fcurves)
    fcurves2 = os.path.splitext(os.path.basename(csvfile))[0] + ".HitsvTotal.tocs.csv"
    fcurves2 = os.path.join(outdir, fcurves2)
    fareas = os.path.splitext(os.path.basename(csvfile))[0] + ".TPRvFPR.roc_areas.csv"
    fareas = os.path.join(outdir, fareas)

    if args.plot_formats is None:
        plotformats = {'png'}
    else:
        plotformats = set()
        for element in args.plot_formats:
            if element.lower() in {'png', 'eps'}:
                plotformats.add(element.lower())

    # Parsing scores
    scores = []
    wholescores = {}
    names = []
    ignore = set()
    srcs = tuple(pco._sourcenames(short=True))

    if 'cropped' in os.path.basename(csvfile):
        crpd = 'cropped'
        crp = True
    elif 'full' in os.path.basename(csvfile):
        crpd = False
        crp = False
    else:
        crpd = None
        crp = None

    pic.csvheader(fcurves, cropped=crp, csvtype='rocs')
    pic.csvheader(fcurves, cropped=crp, csvtype='tocs')
    pic.csvheader(fareas, cropped=crp, csvtype='rocareas')

    with open(csvfile, newline='') as csvin:
        signals = csv.reader(csvin, delimiter=',', quotechar='|')
        col = []
        for row in signals:
            if row[0].startswith('#') is False:
                if row[-1].lower().strip() == 'true':
                    print('true')
                    wholescores['PISAscore'].append(1.0)
                elif row[-1].lower().strip() == 'false':
                    wholescores['PISAscore'].append(0.0)
                else:
                    wholescores['PISAscore'].append(None)
                c = 0
                for n in range(len(row)-1):
                    if col[n] is not None:
                        scores.append([[], []])
                        if row[n].lower().strip() != 'nan':
                            if names[c] == 'MCC_cropseq_dmp':
                                print('pares')
                            scores[c][0].append(float(row[n]))
                            wholescores[names[c]].append(float(row[n]))
                            if row[-1].lower().strip() == 'true':
                                scores[c][1].append(True)
                            elif row[-1].lower().strip() == 'false':
                                scores[c][1].append(False)
                            else:
                                scores[c][1].append(None)
                        else:
                            if names[c] == 'MCC_cropseq_dmp':
                                print('nones')
                            wholescores[names[c]].append(None)
                        c += 1
            else:
                if row[0] == '#PDB_id' :
                    for n in range(len(row)):
                        if row[n].strip().endswith(srcs):
                            col.append(row[n].strip())
                            names.append(row[n].strip())
                            wholescores[row[n].strip()] = []
                        elif row[n].strip().lower() == 'pisascore':
                            wholescores['PISAscore'] = []
                        else:
                            col.append(None)

    print('look at this:')
    print(wholescores['PISAscore'], len(wholescores['PISAscore']))
    # Calculate ROCs, TOCs areas and correlations
    # TO DO: Calcultate TOCs too! https://en.wikipedia.org/wiki/Total_operating_characteristic
    L = len(names)
    print(names)
    rates = {}
    tocs = {}
    unsrtdareas = []  # 0.5*(TPR[n]-TPR[n-1])*(FPR[n]+FPR[n-1])

    for n in range(L):
        area = 0
        rates[names[n]] = [[], []]  # FPR, TPR
        rates[names[n]][0], rates[names[n]][1], area = pcs.tpr_vs_fpr(scores[n][0], scores[n][1])
        unsrtdareas.append(area)

        tocs[names[n]] = [[], []]  # Tots, Hits
        tocs[names[n]][0], tocs[names[n]][1] = pcs.hits_vs_total(scores[n][0], scores[n][1])

    areas, names = zip(*sorted(zip(unsrtdareas, names), reverse=True))
    areas_dict = {names[i]: areas[i] for i in range(len(names))}

    # Calculate correlation matrices
    correl_matrix = np.identity(L+1)

    namex = tuple(['PISAscore'] + list(names))

    for n in range(len(namex)-1):
        for m in range(n+1, len(namex)):
            set1 = []
            set2 = []
            setr = []
            print(namex[n], len(wholescores[namex[n]]))
            print(namex[m], len(wholescores[namex[m]]))
            for dat in range(len(wholescores[namex[n]])):
                if (wholescores[namex[m]][dat] is not None and
                        wholescores[namex[n]][dat] is not None):
                    set1.append(wholescores[namex[m]][dat])
                    set2.append(wholescores[namex[n]][dat])
                    if wholescores[namex[1]][dat] is None:
                        setr.append(float('-inf'))
                    else:
                        setr.append(wholescores[namex[1]][dat])

            if len(set1) > 0:
                correlation = pcs.correl_matrix(set1, set2, setref=setr)
                correl_matrix[n][m] = correlation[0][1]
                correl_matrix[m][n] = correlation[1][0]
            else:
                correl_matrix[n][m] = 0


    # Print out results
    for n in range(L):
        f2 = fcurves2.replace("replaceme", names[n])
        pic.csvheader(f2, cropped=crp, csvtype='tocs') # MUST BE ONE BY ONE AND HAVE OWN HEADER

        pic.lineout([names[n], areas[n]], fareas)

    ignore = ()
    p = 0
    while True:
        listline = []
        for name in names:
            f2 = fcurves2.replace("replaceme", name)
            if name in ignore:
                listline.append("")
                listline.append("")
            else:
                listline.append(rates[name][0][p])
                listline.append(rates[name][1][p])
                if len(rates[name][0]) == p + 1:
                    ignore.add(name)

            pic.lineout(listline, f2)
        pic.lineout(listline, fcurves)
        p += 1
        if len(ignore) == len(names):
            break

    fname = os.path.splitext(os.path.basename(csvfile))[0]
    fout = os.path.join(outdir, fname + '.correlations.csv')
    pic.npmatrix(inmatrix = correl_matrix, outpath=fout)

    endmsg = pcl.ok(starttime, command=__script__)
    logger.info(endmsg)

    for imtype in plotformats:
        fout = os.path.join(outdir, fname + '.correlations.' + imtype )
        pip.plot_correlation_heatmap(data=correl_matrix, outpath=fout, tags=namex, plot_type=imtype)

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
