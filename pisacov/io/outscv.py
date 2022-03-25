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
    scnames = psc._scorenames(crop=cropped)
    croptag = 'cropseq' if cropped is True else 'fullseq'
    csvline = ('#PDB_id, Interface, Chain1, Chain2, Sequence, ' +
               'L' + croptag + ', ' +
               'Ncrops' + 'Lfullseq' +
               croptag + 'Depth, ',
               'N' + croptag + '_contacts, ' +
               'N' + croptag + '_usedcontacts, ')
    for score in scnames:
        csvline += score + ', '

    csvline = csvline[:-2]

    with open(outpath, "w") as outcsv:
        outcsv.write(csvline)

    return
