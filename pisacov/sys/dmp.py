"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.iomod.conf import DMP_PATH
from pisacov.iomod import paths as ppaths

import os
import logging


def rundmp(spath, msapath, outdir):
    """Run DeepMetaPSICOV to produce contact prediction lists.

    :param seqpath: Input sequence filepath.
    :type seqpath: str
    :param msapath: Input MSA filepath.
    :type msapath: str
    """
    dmp_exec = '"'+ppaths.check_path(DMP_PATH, 'file')+'"'
    try:
        oc = os.system(dmp_exec + ' -i ' + spath + ' -a ' + msapath)
        if oc != 0:
            raise Exception
    except Exception:
        logging.critical('        An error occurred while executing DeepMetaPSICOV.')
        raise OSError
