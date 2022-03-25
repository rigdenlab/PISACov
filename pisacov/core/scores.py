"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import _conf_ops as pco


def _scorenames(crop=False):
    """
    Return scorenames.

    :param crop: Include tag for cropped sequence, defaults to False
    :type crop: bool, optional
    :return: scorenames
    :rtype: dict [str: list [str]]

    """
    croptag = 'fullseq' if crop is False else 'cropseq'
    names = pco._sourcenames()
    shortnames = pco._sourcenames(short=True)
    scorenames = {}
    mainname = ['AVScoreRaw', 'AVScoreNorm', 'ACCScoreRaw', 'ACCScoreNorm',
                'TP', 'PREC', 'COVER', 'MCC', 'JAC']
    for i in range(len(names)):
        if names[i] not in scorenames:
            scorenames[names[i]] = []
        for mn in mainname:
            if mn[-4:] == 'Norm' and names[i] != 'psicov':
                pass
            else:
                scorenames[names[i]].append(mainname + '_' +
                                            croptag + '_' +
                                            shortnames[i])

    return scorenames


def avraw(inmap):

    return outval


def avnorm(inmap):

    return outval


def accraw(inmap):

    return outval


def accnorm(inmap):

    return outval


def n_tps(inmap):

    return outval


def precision(inmap):

    return outval


def coverage(inmap):

    return outval


def mcc(inmap):

    return outval


def jaccard(inmap):

    return outval
