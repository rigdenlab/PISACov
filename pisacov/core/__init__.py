"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__


def _psicov_modes():
    """Return a list with PSICOV modes used.

    :return: Psicov Modes
    :rtype: list [str]

    """
    return ['norm', 'abs', 'shifted']
