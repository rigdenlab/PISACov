"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io.conf import DMP_PATH
from pisacov.io import paths as ppaths

import os
import crops
import sys
import logging

def rundmp(spath, msapath, outdir):
    """
    Run DeepMetaPSICOV to produce contact prediction lists.

    :param seqpath: Input sequence filepath.
    :type seqpath: str
    :param msapath: Input MSA filepath.
    :type msapath: str
    """
    dmp_exec = '"'+ppaths.check_path(DMP_PATH, 'file')+'"'
    try:
        os.system(dmp_exec + ' -i '+ spath + ' -a ' + msapath)
    except:
        logging.critical('        An error occurred while executing DeepMetaPSICOV.')
        raise OSError
