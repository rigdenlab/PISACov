"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import logging

def _ftypelist():
    lout = ['psicov', 'fasta', 'pdb', 'a3m']
    return lout

def read(infile, ftype):
    """

    :param infile: Path to input file.
    :type infile: str
    :param ftype: One of ['psicov', 'fasta', 'pdb', 'a3m']
    :type ftype: str
    :return: Parsed file.
    :rtype: ....

    """
    if (ftype.lower() not in _ftypelist() or
            isinstance(ftype, str) is not True):
        logging.critical('Specified type not valid.')
        raise ValueError

    if ftype.lower() == 'psicov':
        pass
    elif ftype.lower() == 'fasta':
        pass
    elif ftype.lower() == 'pdb':
        pass
    elif ftype.lower() == 'a3m':
        pass
