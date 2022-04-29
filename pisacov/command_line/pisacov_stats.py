"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__
__script__ = 'PISACov Statistical Analysis script'

from pisacov import command_line as pcl
from pisacov.io import paths as ppaths

import argparse
import datetime
import os
import csv

logger = None

def create_argument_parser():
    """Create a parser for the command line arguments used in pisacov_stats."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__ + os.linesep + __doc__,
                                     epilog="Check pisacov.rtfd.io for more information.")
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_rocs = subparsers.add_parser('rocs', help='Produce Receiver operating characteristic (ROC) curves for the data provided.')
    parser_rocs.add_argument('scores', nargs=1, metavar=("ScoresCSVFile"),
                        help="Input scores CSV filepath.")
    parser_rocs.add_argument("-f", "--full_score_analysis", action='store_true',
                             default=False,
                             help="Produce full analysis of beta score list (beta score list required).")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input data.")

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
    outdir = ppaths.check_path(args.outdir)
    ppaths.mdir(outdir)

    # Parsing scores
    scores = {}
    names = None
    thraw = {}
    with open(csvfile, 'r') as fin:
        scoresin = csv.reader(fin)
        for entry in scoresin:
            if entry[0][0] != "#":
                if (names is None or isinstance(names, str) or
                        (isinstance(names, list) and
                         len(names) != len(entry))):
                    names = []
                    for n in range(len(entry)):
                        names.append('sc_'+str(n+1))
                else:
                    if thraw == {}:
                        for name in names[13:-1]:
                            thraw[name] = []
                if entry[0] not in scores:
                    scores[entry[0]] = {}
                if entry[1] not in scores[entry[0]]:
                    scores[entry[0]][entry[1]] = []
                    for sc in (entry.split(sep=', ')[13:-1]):
                        scores[entry[0]][entry[1]].append(float(sc))
                    if (entry.split(sep=', ')[-1]) == 'True' or '1':
                        scores[entry[0]][entry[1]].append(True)
                    elif (entry.split(sep=', ')[-1]) == 'False' or '0':
                        scores[entry[0]][entry[1]].append(False)
                    for n in range(len(names)):
                        thraw[names[n]].append(scores[entry[0]][entry[1]][n])
                else:
                    if entry.split(sep=', ')[13:] == scores[entry[0]][entry[1]]:
                        pass
                    else:
                        raise ValueError('CSV file contains different values for same interface.')
            else:
                names = entry[1:].split(sep=', ')

    # Setting thresholds
    thr = {}
    FPR = {}
    TPR = {}
    for key, value in thraw.items():
        thr[key] = list(set(thraw)).sort()
        FPR[key] = []
        TPR[key] = []
        for t in thr[key]:
            FP = 0
            TP = 0
            FN = 0
            TN = 0
            for pdbid in scores:
                for iface in scores[pdbid]:
                    stable = scores[pdbid][iface][-1]
                    for n in range(len(names)):
                        if scores[pdbid][iface][n] < t:
                            if stable is True:
                                FN += 1
                            else:
                                TN += 1
                        else:
                            if stable is True:
                                TP += 1
                            else:
                                FP += 1
            FPR[key].append(FP/(FP+TN))
            TPR[key].append(TP/(TP+FN))
        fnameout = os.path.join(outdir, (key +
                                         os.path.splitext(os.path.basename(csvfile))[0] +
                                         'roc.dat'))
        with open(fnameout, 'w') as fout:
            for n in range(len(FPR[key])):
                fout.write(str(FPR[key][n]) + ' ' + str(TPR[key][n]))

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
