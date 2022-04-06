"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import _conf_ops as pco
from pisacov.core import scores as psc

import csv


def csvheader(outpath, cropped=False):
    """
    Create a new CSV file with only the header.

    :param outpath: CSV filepath.
    :type outpath: str
    :param cropped: True when results have been obtained from crops-cropstr inputs, defaults to False.
    :type cropped: bool, optional
    """
    scnames = psc._scorenames(crop=cropped)
    croptag = 'cropseq' if cropped is True else 'fullseq'
    csvline = ('#PDB_id, Interface, Chain1, Chain2, Sequence, ' +
               'L' + croptag + ', ' + 'Neff' + croptag + ', ' +
               'Ncrops, ' + 'Lfullseq, ' + 'Nefffullseq, ', +
               croptag + 'Depth, ',
               'N' + croptag + '_contacts, ' +
               'N' + croptag + '_usedcontacts, ')
    for score in scnames:
        csvline += score + ', '

    csvline = csvline[:-2]

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
