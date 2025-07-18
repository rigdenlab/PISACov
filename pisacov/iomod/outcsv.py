"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.core import scores as psc
from pisacov.iomod import _conf_ops as pco

import csv
import os
import datetime
import logging


def csvheader(outpath, csvtype='scores', cropped=False, pisascore=False,
              upth=False, scoreth=None, singularID=None):
    """
    Create a new CSV file with only the header.

    :param outpath: CSV filepath.
    :type outpath: str
    :param csvtype: 'scores' or 'rocs' or 'tocs' or 'rocareas' or 'rocs_bezier', defaults to 'scores'.
    :type csvtype: str, optional
    :param cropped: True when results have been obtained from crops_cropstr inputs, defaults to False.
    :type cropped: bool, optional
    :param pisascore: If using CCP4 PISA interfaces add to csv header, defaults to False.
    :type pisascore: bool, optional
    :param upth: UniProt single proportion threshold, defaults to False.
    :type upth: float, optional
    :param scoreth: Low Threshold for predicted contacts' scores [psicov, ccmpred, dmp], defaults to False.
    :type scoreth: list [float], optional
    :param singularID: If the file contains information for a single system, include information of the ID in the header, defaults to None.
    :type singularID: str
    """
    begin = datetime.datetime.now()
    tzbegin = begin.astimezone().strftime("%d %B %Y - %H:%M:%S %Z")
    if singularID is None:
        sID = ''
    else:
        sID = singularID + '. '
    if csvtype == 'scores':
        csvline = '# PISACov run. ' + sID + tzbegin + ". Using "
        if cropped is True:
            csvline += "cropped "
        elif cropped is False:
            csvline += "full "
        elif cropped is None:
            csvline += "unknown "
        csvline += "version of sequences. "
        if upth is False:
            upth = 0.0
        if upth > 0.0:
            csvline += "UniProt single sequence proportion threshold = " + str(upth)
        csvline += ". "
        if scoreth is None:
            csvline += "No predicted contacts filtered out according to score."
        else:
            csvline += "Predicted contacts filtered out when score: "
            csvline += "PSICOV < " + str(scoreth[0])
            csvline += "; CCMPred < " + str(scoreth[1])
            csvline += "; DeepMetaPSICOV < " + str(scoreth[2]) +"."
    elif csvtype == 'rocs':
        csvline = "# PISACov stats: TPR vs FPR Receiver operating characteristic curves (ROCs)." + tzbegin + "."
    elif csvtype == 'tocs':
        csvline = "# PISACov stats: Hits (TPs) vs total (TPs+FPs) Total operating characteristic curves (TOCs)." + tzbegin + "."
    elif csvtype == 'rocs_bezier':
        csvline = "# PISACov stats: TPR vs FPR Receiver operating characteristic curves (ROCs) data obtained from Bézier curve (B(t)). " + sID + tzbegin + "."
        csvline += os.linesep
        csvline += "# t_param, Bx(t), By(t), Bx'(t), By'(t), Bx''(t), By''(t), Bx'''(t), By'''(t), "
        csvline += "K(t), LR(t), Scores(t), P(t), lambda(t), Youden(t), AUC"
    elif csvtype == 'rocareas':
        csvline = "# PISACov stats: Areas under TPR vs FPR Receiver operating characteristic curves (ROC areas). Sorted by area." + tzbegin + "."
        csvline += os.linesep
        csvline += "# Score, Area under ROC"
    else:
        logging.critical("        pisacov.iomod.csvheader input 'csvtype' must be one of scores' or rocs' or 'tocs' or 'rocareas'.")
        raise ValueError

    csvline += os.linesep

    if csvtype == 'scores':
        scnames = psc._scorenames(crop=cropped)
        croptag = 'cropseq' if cropped is True else 'fullseq'
        names = pco._sourcenames(short=True)

        csvline += ('#PDB_id, Interface, Chain1, Chain2, Sequence, ' +
                   'L' + croptag + ', ' + 'Neff_' + croptag + ', ' +
                   'Ncrops, ' + 'Lfullseq, ' + 'Neff_fullseq, ' +
                   'N' + croptag + '_IFcontacts, ' +
                   'N' + croptag + '_IFcontacts_used, ')

        for source in names:
            for score in scnames[source]:
                csvline += score + ', '

        csvline = csvline[:-2]

        if pisascore is True and csvtype == 'scores':
            csvline += ', PISAscore'

        csvline += os.linesep

    with open(outpath, "w") as outcsv:
        outcsv.write(csvline)

    return


def lineout(alist, outpath):
    """
    Print sets of results to output csv file.

    :param alist: List containing all the data.
    :type alist: list [str]
    :param outpath: Output filepath.
    :type outpath: str

    """
    with open(outpath, "a", newline='') as f:
        csv_writer = csv.writer(f, delimiter=',', quotechar='#',
                                quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(alist)

    return
