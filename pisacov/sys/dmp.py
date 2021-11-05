"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov.about import __prog__, __description__, __version__
from pisacov.about import  __author__, __date__, __copyright__

from pisacov.io.conf import DMP_PATH
from pisacov import io as pio

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
    dmp_exec = '"'+pio.check_path(DMP_PATH, 'file')+'"'
    try:
        os.system(dmp_exec + ' -i '+ spath + ' -a ' + msapath)
    except:
        logging.critical('        An error occurred while executing DeepMetaPSICOV.')