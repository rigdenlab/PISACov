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
    mainname = ['AVScoreRaw', 'AVScoreNorm', 'AVScoreAbs',
                'ACCScoreRaw', 'ACCScoreNorm', 'ACCScoreAbs',
                'TP', 'PREC', 'COVER', 'MCC', 'JAC']
    for i in range(len(names)):
        if names[i] not in scorenames:
            scorenames[names[i]] = []
        for mn in mainname:
            if mn[-4:] == 'Norm' and names[i] != 'psicov':
                pass
            elif mn[-3:] == 'Abs' and names[i] != 'psicov':
                pass
            else:
                scorenames[names[i]].append(mainname + '_' +
                                            croptag + '_' +
                                            shortnames[i])

    return scorenames


def avscore(inmap, mode=None):
    """
    Return the average score of the True Positive contacts.

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :param mode: None (raw), or 'absolute', or 'normal', or 'shifted', defaults to None.
    :type mode: TYPE, optional
    :return: Cumulative value
    :rtype: float

    """

    av = 0.0

    nc = float(len(inmap))
    if mode.lower() == 'shifted' or mode.lower() == 'normal':
        # Introduce minimum and maximum value calculation

    for contact in inmap:
        if contact.true_positive is True:
            if mode is None or mode.lower() == 'raw':
                av += contact.raw_score
            elif mode.lower() == 'absolute':
                av += abs(contact.raw_score)
            elif mode.lower() == 'shifted':
                av +=

    av = av / nc

    return av


def accscore(inmap, mode=None):
    """
    Return the cumulative score of the True Positive contacts.

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :param mode: None (raw), or 'absolute', or 'normal', defaults to None.
    :type mode: TYPE, optional
    :return: Cumulative value
    :rtype: float

    """

    acc = 0.0

    nc = float(len(inmap))
    for contact in inmap:
        if contact.true_positive is True:
            acc += contact.raw_score

    return acc


def n_tps(inmap):
    """
    Return the number of True Positive contacts.

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Number of true positives
    :rtype: int

    """

    ntp = 0

    for contact in inmap:
        if contact.true_positive is True:
            ntp += 1

    return ntp


def n_fps(inmap):
    """
    Return the number of False Positive contacts.

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Number of false positives
    :rtype: int

    """

    nfp = 0

    for contact in inmap:
        if contact.false_positive is True:
            nfp += 1

    return nfp


def n_tns(inmap):
    """
    Return the number of True Negative contacts.

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Number of true negatives
    :rtype: int

    """

    ntn = 0

    for contact in inmap:
        if contact.true_negative is True:
            ntn += 1

    return ntn


def n_fns(inmap):
    """
    Return the number of False Negative contacts.

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Number of false negatives
    :rtype: int

    """

    nfn = 0

    for contact in inmap:
        if contact.false_negative is True:
            nfn += 1

    return nfn


def precision(inmap):
    """
    Return the precision of a given matched map. Prec = TP / (TP+FP).

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Precision
    :rtype: float

    """
    ntp = float(n_tps(inmap))
    nfp = float(n_fps(inmap))

    prec = ntp / (ntp + nfp)

    return prec


def coverage(inmap):
    """
    Return the coverage of a given matched map. Cover = TP / (TP+FN).

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Coverage
    :rtype: float

    """
    ntp = float(n_tps(inmap))
    nfn = float(n_fns(inmap))

    cover = ntp / (ntp + nfn)

    return cover


def mcc(inmap):
    """
    Return the Matthew’s Correlation Coefficient (MCC) of a given matched map. MCC = TP / (TP+FN).

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: MCC
    :rtype: float

    """
    ntp = float(n_tps(inmap))
    nfp = float(n_fps(inmap))
    nfn = float(n_fns(inmap))
    ntn = float(n_tns(inmap))

    mcc = (ntp*ntn - nfp*nfn) / ((ntp + nfp)*(ntp + nfn)*(ntn + nfp)*(ntn + nfn))

    return mcc


def jaccard(inmap):
    """
    Return the Jaccard Index of a given matched map. Jaccard = TP / (TP+FP+FN).

    :param inmap: Conkit matched map.
    :type inmap: conkit map
    :return: Jaccard Index
    :rtype: float

    """
    ntp = float(n_tps(inmap))
    nfp = float(n_fps(inmap))
    nfn = float(n_fns(inmap))

    jacc = ntp / (ntp + nfp + nfn)

    return jacc
