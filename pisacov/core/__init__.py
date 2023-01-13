"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__


def _psicov_modes(raw=False):
    """Return a list with PSICOV modes used.

    :param raw: Include 'raw' as the first mode in list, defaults to False.
    :type raw: bool, optional
    :type raw: Include 'raw' as the first mode in list, defaults to False.

    :return: Psicov Modes
    :rtype: list [str]

    """
    if raw:
        retlist = ['raw', 'norm', 'abs', 'shifted']
    else:
        retlist = ['norm', 'abs', 'shifted']

    return retlist
