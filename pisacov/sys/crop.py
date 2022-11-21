"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import os
import importlib
import sys
import logging

def runcrops(seqin, strin, dbin, thin=None, upin=None, outdirin = None):#, loggingfile):
    """Run CROPS to produce clean sequence and structure files.

    :param seqin: Sequence filepath.
    :type seqin: str
    :param strin: Structure filepath.
    :type strin: str
    :param dbin: SIFTS database filepath.
    :type dbin: str
    :param thin: Uniprot threshold, defaults to None.
    :type thin: str, int, optional
    :param upin: Uniclust filepath (or server-only url), defaults to None.
    :type upin: str, int, optional
    :param outdirin: Output directory's path. If not given, seqin dir will use instead, defaults to None.
    :type outdirin: str, optional

    """

    #cropsdir = os.path.dirname(crops.__file__)
    pythonexec = '"'+sys.executable+'"'

    #CROP SEQUENCE AND STRUCTURE

    cst = importlib.import_module('crops.command_line.crops_cropstr')
    cropspy = (cst.__file__)

    if outdirin is None:
        outdirin = os.dirname(seqin)

    command = pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + strin + ' ' + dbin
    if thin is not None:
        command += ' -u ' + thin + ' ' + upin
    command += ' -o ' + outdirin + ' -i'
    try:
        os.system(command)  # + ' > ' + loggingfile)
    except Exception:
        logging.critical('        An error occurred while executing crops_cropstr')
        raise OSError

    return


def renumcrops(seqin, strin, outdirin = None):#, loggingfile):
    """Run CROPS to produce clean sequence and structure files.

    :param seqin: Sequence filepath.
    :type seqin: str
    :param strin: Structure filepath.
    :type strin: str
    :param outdirin: Output directory's path. If not given, seqin dir will use instead, defaults to None.
    :type outdirin: str, optional

    """

    #cropsdir = os.path.dirname(crops.__file__)
    pythonexec = '"'+sys.executable+'"'

    #Renumber STRUCTURE
    logging.info('    Running crops_renumber...')

    cst = importlib.import_module('crops.command_line.crops_renumber')
    cropspy = (cst.__file__)

    if outdirin is None:
        outdirin = os.dirname(seqin)

    command = (pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + strin + ' ' +
               ' -o ' + outdirin)
    try:
        os.system(command)  # + ' > ' + loggingfile)
    except Exception:
        logging.critical('        An error occurred while executing crops_renumber')
        raise OSError

    return


def splitseqs(seqin, outdirin = None):#, loggingfile):
    """Run CROPS to split sequence files.

    :param seqin: Sequence filepath.
    :type seqin: str
    :param outdirin: Output directory's path. If not given, seqin dir will use instead, defaults to None.
    :type outdirin: str, optional

    """

    #cropsdir = os.path.dirname(crops.__file__)
    pythonexec = '"'+sys.executable+'"'

    #Renumber STRUCTURE
    logging.info('    Running crops_renumber...')

    cst = importlib.import_module('crops.command_line.crops_splitseqs')
    cropspy = (cst.__file__)

    if outdirin is None:
        outdirin = os.dirname(seqin)

    command = (pythonexec + ' ' + cropspy + ' ' + seqin + ' '  + ' -i'
               ' -o ' + outdirin)
    try:
        os.system(command)  # + ' > ' + loggingfile)
    except Exception:
        logging.critical('        An error occurred while executing crops_renumber')
        raise OSError

    return
