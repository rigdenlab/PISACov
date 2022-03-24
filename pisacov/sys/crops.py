"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import os
# import crops
import importlib
import sys
import logging

def runcrops(seqin, strin, dbin, thin=None, upin=None, outdirin = None):#, loggingfile):
    """
    Run CROPS to produce clean sequence and structure files.

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
    logging.info('    Running crops-cropstr...')

    cst = importlib.import_module('crops.command_line.crops-cropstr')
    cropspy = (cst.__file__)

    command = pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + strin + ' ' + dbin
    if thin is not None:
        command += ' -u ' + thin + ' ' + upin
    command += ' -o ' + outdirin + ' -i'
    try:
        os.system(command)  # + ' > ' + loggingfile)
    except:
        logging.critical('        An error occurred while executing Crops-cropstr')

    logging.info('    Done' + os.linesep)

    return
