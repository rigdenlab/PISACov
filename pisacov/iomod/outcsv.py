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


def csvheader(outpath, cropped=False, pisascore=False,
              upth=False, scoreth=False):
    """
    Create a new CSV file with only the header.

    :param outpath: CSV filepath.
    :type outpath: str
    :param cropped: True when results have been obtained from crops_cropstr inputs, defaults to False.
    :type cropped: bool, optional
    :param pisascore: If using CCP4 PISA interfaces add to csv header, defaults to False.
    :type pisascore: bool, optional
    :param upth: UniProt single proportion threshold, defaults to False.
    :type upth: float, optional
    :param scoreth: Low Threshold for predicted contacts' scores [psicov, ccmpred, dmp], defaults to False.
    :type scoreth: list [float], optional
    """
    csvline = '# PISACov run. '
    begin = datetime.datetime.now()
    tzbegin = begin.astimezone().strftime("%d %B %Y - %H:%M:%S %Z")
    csvline += tzbegin
    csvline += ". Using "
    if cropped:
        csvline += "cropped "
    else:
        csvline += "full "
    csvline += "version of sequences. "
    if upth is False:
        upth = 0.0
    if upth > 0.0:
        csvline += "UniProt single sequence proportion threshold = " + str(upth)
    csvline += ". "
    if scoreth is False:
        csvline += "No predicted contacts filtered out according to score."
    else:
        csvline += "Predicted contacts filtered out when score below: "
        csvline += "PSICOV = " + str(scoreth[0])
        csvline += "; CCMPred = " + str(scoreth[1])
        csvline += "; DeepMetaPSICOV = " + str(scoreth[2]) +"."

    csvline += os.linesep

    scnames = psc._scorenames(crop=cropped)
    names = pco._sourcenames()
    croptag = 'cropseq' if cropped is True else 'fullseq'
    names = pco._sourcenames()
    csvline = ('#PDB_id, Interface, Chain1, Chain2, Sequence, ' +
               'L' + croptag + ', ' + 'Neff_' + croptag + ', ' +
               'Ncrops, ' + 'Lfullseq, ' + 'Neff_fullseq, ' +
               'N' + croptag + '_IFcontacts, ' +
               'N' + croptag + '_IFcontacts_used, ')
    for source in names:
        for score in scnames[source]:
            csvline += score + ', '

    csvline = csvline[:-2]

    if pisascore is True:
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
