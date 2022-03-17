"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import os
import crops
import sys
import logging

def runcrops(seqin, strin, dbin, thin, upin, outdirin):#, loggingfile):
    """
    Run CROPS to produce clean sequence and structure files.

    :param seqin: Sequence filepath.
    :type seqin: str
    :param strin: Structure filepath.
    :type strin: str
    :param dbin: SIFTS database filepath.
    :type dbin: str
    :param outdirin: Output directory's path.
    :type outdirin: str

    """

    cropsdir = os.path.dirname(crops.__file__)
    pythonexec = '"'+sys.executable+'"'

    # CROP SEQUENCE
    logging.info('    Running crops-cropseq...')
    cropspy = os.path.join(cropsdir, 'crops', 'command_line', 'crops-cropseq.py')
    try:
        os.system(pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + dbin + ' ' +
                  '-u ' + thin + ' ' + upin + ' -o ' + outdirin)  # + ' > ' + loggingfile)
    except:
        logging.critical('        An error occurred while executing Crops-cropseq')

    logging.info('    Done\n')

    #CROP STRUCTURE
    logging.info('    Running crops-cropstr...')
    cropspy = os.path.join(cropsdir, 'crops', 'command_line', 'crops-cropstr.py')
    try:
        os.system(pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + strin + ' ' +
                  dbin + ' ' + '-u ' + thin + ' ' + upin + ' -o ' + outdirin)  # + ' > ' + loggingfile)
    except:
        logging.critical('        An error occurred while executing Crops-cropstr')

    logging.info('    Done\n')

    return
