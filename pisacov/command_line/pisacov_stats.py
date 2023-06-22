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
    parser.add_argument("-n", "--none_is_false_pisa", action='store_true',
                        default=False,
                        help="Include additional data with PISA None results considered False.")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input data.")
    parser.add_argument("-p", "--plot_formats", nargs='+',
                        metavar=("Plot file format(s)"),
                        help="One or more formats of 'png' and 'eps' and 'svg' of figures to be produced.")

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
    #fcurves2 = os.path.splitext(os.path.basename(csvfile))[0] + ".HitsvTotal.tocs.csv"
    #fcurves2 = os.path.join(outdir, fcurves2)
    fcurvesB = "replaceme." + os.path.splitext(os.path.basename(csvfile))[0] + ".TPRvFPR.rocs_bezier.csv"
    fcurvesB = os.path.join(outdir, fcurvesB)
    fareas = os.path.splitext(os.path.basename(csvfile))[0] + ".TPRvFPR.roc_areas.csv"
    fareas = os.path.join(outdir, fareas)
    #fareas2 = os.path.splitext(os.path.basename(csvfile))[0] + ".TPRvFPR.roc_areas.unsorted.csv"
    #fareas2 = os.path.join(outdir, fareas2)

    if args.plot_formats is None:
        plotformats = {'png'}
    else:
        plotformats = set()
        for element in args.plot_formats:
            if element.lower() in {'png', 'eps', 'svg'}:
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
    pic.csvheader(fcurvesB, cropped=crp, csvtype='rocs_bezier')
    #pic.csvheader(fcurves, cropped=crp, csvtype='tocs')
    pic.csvheader(fareas, cropped=crp, csvtype='rocareas')
    #pic.csvheader(fareas2, cropped=crp, csvtype='rocareas')

    with open(csvfile, newline='') as csvin:
        signals = csv.reader(csvin, delimiter=',', quotechar='|')
        col = []
        for row in signals:
            if row[0].startswith('#') is False:
                if row[-1].lower().strip() == 'true':
                    #print('true')
                    wholescores['PISAscore_NN'].append(1.0)
                    if args.none_is_false_pisa is True:
                        wholescores['PISAscore_NF'].append(1.0)
                elif row[-1].lower().strip() == 'false':
                    wholescores['PISAscore_NN'].append(0.0)
                    if args.none_is_false_pisa is True:
                        wholescores['PISAscore_NF'].append(0.0)
                else:
                    wholescores['PISAscore_NN'].append(None)
                    if args.none_is_false_pisa is True:
                        wholescores['PISAscore_NF'].append(0.0)
                c = 0
                for n in range(len(row)-1):
                    if col[n] is not None:
                        scores.append([[], []])
                        if row[n].lower().strip() != 'nan':
                            #if names[c] == 'MCC_cropseq_dmp':
                                #print('pares')
                            scores[c][0].append(float(row[n]))
                            wholescores[names[c]].append(float(row[n]))
                            if row[-1].lower().strip() == 'true':
                                scores[c][1].append(True)
                            elif row[-1].lower().strip() == 'false':
                                scores[c][1].append(False)
                            else:
                                scores[c][1].append(None)
                        else:
                            #if names[c] == 'MCC_cropseq_dmp':
                                #print('nones')
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
                            wholescores['PISAscore_NN'] = []
                            if args.none_is_false_pisa is True:
                                wholescores['PISAscore_NF'] = []
                        else:
                            col.append(None)

    #print('look at this:')
    #print(wholescores['PISAscore'], len(wholescores['PISAscore']))
    # Calculate ROCs, TOCs areas and correlations
    L = len(names)
    scores_dict = {names[i]: scores[i][0] for i in range(len(names))}
    #print(names)
    rocs = {}
    rocs_scores = {}
    rocs_bezier = {}
    tocs = {}
    unsrtdareas = []  # 0.5*(TPR[n]-TPR[n-1])*(FPR[n]+FPR[n-1])

    names2 = []
    pindex = ("NN", "NF")
    pbool = (False, True)
    for p in range(len(pindex)):
        if args.none_is_false_pisa is False and pindex[p] == "NF":
            break
        for n in range(L):
            names2.append(names[n] + '_' + pindex[p])
            area = 0
            rocs[names2[L*p+n]] = [[], []]  # FPR, TPR
            rocs[names2[L*p+n]][0], rocs[names2[L*p+n]][1], rocs_scores[names2[L*p+n]], area = (
                pcs.tpr_vs_fpr(scores[n][0], scores[n][1], noneisfalse=pbool[p]))
            unsrtdareas.append(area)

            #tocs[names2[L*p+n]] = [[], []]  # Tots, Hits
            #tocs[names2[L*p+n]][0], tocs[names2[L*p+n]][1] = (
            #    pcs.hits_vs_total(scores[n][0], scores[n][1], noneisfalse=pbool[p]))

    unsrtdnames = copy.deepcopy(names2)
    areas, names2 = zip(*sorted(zip(unsrtdareas, names2), reverse=True))
    areas_dict = {names2[i]: areas[i] for i in range(len(names2))}
    areas_dict_best = {}
    for name, auc in areas_dict.items():
        if auc >= 0.7 or auc <= 0.3:
            areas_dict_best[name] = auc
            isconvex = True if auc < 0.5 else False
            rocs_bezier[name] = pcs.bezier_parametrization([rocs[name][0],
                                                            rocs[name][1]],
                                                           scores=rocs_scores[name],
                                                           convex=isconvex)

    # Calculate correlation matrices
    names_best = []
    for key, auc in areas_dict.items():
        for name in names:
            if key.startswith(name):
                if auc >= 0.7 or auc <= 0.3:
                    if name not in names_best:
                        names_best.append(name)

    Lbest = len(names_best)
    if args.none_is_false_pisa is True:
        correl_matrix = np.identity(L+2)
        correl_matrix[0][1] = 1.
        correl_matrix[1][0] = 1.
        namex = tuple(['PISAscore_NN', 'PISAscore_NF'] + list(names))
        correl_matrix_best = np.identity(Lbest+2)
        correl_matrix_best[0][1] = 1.
        correl_matrix_best[1][0] = 1.
        namex_best = tuple(['PISAscore_NN', 'PISAscore_NF'] + list(names_best))
    else:
        correl_matrix = np.identity(L+1)
        namex = tuple(['PISAscore_NN'] + list(names))
        correl_matrix_best = np.identity(Lbest+1)
        namex_best = tuple(['PISAscore_NN'] + list(names_best))

    for n in range(len(namex)-1):
        for m in range(n+1, len(namex)):
            if args.none_is_false_pisa is True and n == 0 and m == 1:
                continue
            set1 = []
            set2 = []
            setr = []
            for dat in range(len(wholescores[namex[n]])):
                if (wholescores[namex[m]][dat] is not None and
                        wholescores[namex[n]][dat] is not None):
                    set1.append(wholescores[namex[m]][dat])
                    set2.append(wholescores[namex[n]][dat])
                    if wholescores[namex[0]][dat] is None:
                        setr.append(float('-inf'))
                    else:
                        setr.append(wholescores[namex[0]][dat])

            if len(set1) > 0:
                correlation = pcs.correl_matrix(set1, set2, setref=setr)
                correl_matrix[n][m] = correlation[0][1]
                correl_matrix[m][n] = correlation[1][0]
                if namex[m] in namex_best and namex[n] in namex_best:
                    nb = namex_best.index(namex[n])
                    mb = namex_best.index(namex[m])
                    correl_matrix_best[nb][mb] = correlation[0][1]
                    correl_matrix_best[mb][nb] = correlation[1][0]
            else:
                correl_matrix[n][m] = 0.
                if namex[m] in namex_best and namex[n] in namex_best:
                    nb = namex_best.index(namex[n])
                    mb = namex_best.index(namex[m])
                    correl_matrix_best[nb][mb] = 0.

    # Print out results
    for n in range(len(names2)):
        pic.lineout([names2[n], areas[n]], fareas)
        # pic.lineout([unsrtdnames[n], unsrtdareas[n]], fareas2)
        # f2 = fcurves2.replace("replaceme", names2[n])
        # pic.csvheader(f2, cropped=crp, csvtype='tocs') # MUST BE ONE BY ONE AND HAVE OWN HEADER

    ignore = []
    p = 0
    while True:
        listline = []
        for name in names2:
            if name in ignore:
                listline.append("")
                listline.append("")
                listline.append("")
            else:
                listline.append(rocs_scores[name][p])
                listline.append(rocs[name][0][p])
                listline.append(rocs[name][1][p])
                if len(rocs[name][0]) == p + 1:
                    ignore.append(name)

        pic.lineout(listline, fcurves)
        p += 1
        if len(ignore) == len(names2):
            break

    var = ["param", "bezier", "bezier_der1", "bezier_der2", "bezier_der3",
           "curvature", "LR", "scores", "probability", "lambda", "Youden",
           "area"]
    for name in areas_dict_best:
        froc = fcurvesB.replace("replaceme", name)
        for n in range(len(rocs_bezier[name]["param"])):
            listline = []
            for v in var:
                if v == "area" and n > 0:
                    listline.append("")
                elif v.startswith("bezier"):
                    listline.append(rocs_bezier[name][v][0][n])
                    listline.append(rocs_bezier[name][v][1][n])
                else:
                    listline.append(rocs_bezier[name][v][n])
                    pic.lineout(listline, froc)
        pic.lineout(listline, froc)


    fname = os.path.splitext(os.path.basename(csvfile))[0]
    fout = os.path.join(outdir, fname + '.correlations.csv')
    # pic.npmatrix(inmatrix = correl_matrix, outpath=fout)

    for imtype in plotformats:
        pdir = os.path.join(outdir, imtype, '')
        ppaths.mdir(pdir)
        # ROC AREAS - HISTOGRAM
        fout = os.path.join(pdir, fname + '.roc_areas.' + imtype)
        pip.plot_area_histogram(data=areas_dict, outpath=fout, plot_type=imtype)
        # ROC CURVES
        fout = os.path.join(pdir, fname + '.rocs.' + imtype)
        pip.plot_rocs(data=rocs, outpath=fout, areas_for_color=areas_dict, plot_type=imtype)
        fout = os.path.join(pdir, fname + '.norandom.rocs.' + imtype)
        pip.plot_rocs(data=rocs, outpath=fout, areas_for_color=areas_dict, plot_type=imtype, norand=True)
        # ROC CURVES with BÉZIER APPROXIMATION
        for name in areas_dict_best.items():
            fout = os.path.join(pdir, fname + '.' + name + '.bezier.roc.' + imtype)
            pip.plot_roc_parametrization(data=rocs[name], bezier=rocs_bezier[name],
                                         outpath=fout, datatag=name,
                                         area_for_color=areas_dict_best[name],
                                         plot_type=imtype)
            fout = os.path.join(pdir, fname + '.' + name + '.bezier.probvsscore.' + imtype)
            pip.plot_roc_parametrization(data=rocs[name], bezier=rocs_bezier[name],
                                         outpath=fout, datatag=name,
                                         area_for_color=areas_dict_best[name],
                                         plot_type=imtype,
                                         roc_type='probvsscore',
                                         scores=rocs_scores[name])
        # TOC CURVES
        #for name in names2:
        #    fout = os.path.join(pdir, fname + '.' + name + '.toc.' + imtype)
        #    pip.plot_toc(data=tocs[name], datatag=name, outpath=fout,
        #                 area_for_color=areas_dict[name], plot_type=imtype)
        # CORRELATION HEATMAP
        fout = os.path.join(pdir, fname + '.correlations.' + imtype)
        pip.plot_correlation_heatmap(data=correl_matrix, outpath=fout,
                                     labels=namex, plot_type=imtype, light0=True)
        fout = os.path.join(pdir, fname + '.norandom.correlations.' + imtype)
        pip.plot_correlation_heatmap(data=correl_matrix_best, outpath=fout,
                                     labels=namex_best, plot_type=imtype, light0=True)
        #fout = os.path.join(pdir, fname + '.correlations2.' + imtype)
        #pip.plot_correlation_sns(data=correl_matrix, outpath=fout,
        #                         labels=namex, plot_type=imtype, light0=True,
        #                         clustered=False)
        # CLUSTERED CORRELATIONS - HEATMAP + DENDROGRAM
        fout = os.path.join(pdir, fname + '.clusters.' + imtype)
        pip.plot_correlation_sns(data=correl_matrix, outpath=fout,
                                 labels=namex, plot_type=imtype, light0=True,
                                 clustered=True, areas_for_color=areas_dict)

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
