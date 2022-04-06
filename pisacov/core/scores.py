"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import _conf_ops as pco
from pisacov.core import contacts as pcc

import logging
from numpy import sqrt


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
    mainname = ['Nconpred', 'Nconused',
                'ACCScoreRaw', 'AVScoreRaw',
                'ACCScoreNorm', 'AVScoreNorm',
                'ACCScoreAbs', 'AVScoreAbs',
                'ACCScoreShift', 'AVScoreShift',
                'TP', 'PREC', 'COVER', 'MCC', 'JACCARD']
    for i in range(len(names)):
        if names[i] not in scorenames:
            scorenames[names[i]] = []
        for mn in mainname:
            if mn[-4:] == 'Norm' and names[i] != 'psicov':
                pass
            elif mn[-3:] == 'Abs' and names[i] != 'psicov':
                pass
            elif mn[-5:] == 'Shift' and names[i] != 'psicov':
                pass
            else:
                scorenames[names[i]].append(mainname + '_' +
                                            croptag + '_' +
                                            shortnames[i])

    return scorenames


def accscore(inatlas, mode=None):
    """
    Return the cumulative score of the True Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param mode: None (raw), or 'absolute', or 'normal', or 'shifted', defaults to None.
    :type mode: str, optional
    :return: Cumulative value
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('First argument of accscore must be a Contact Atlas.')
        raise TypeError
    if (mode is not None and mode.lower() != 'raw' and
            mode.lower() != 'absolute' and mode.lower() != 'normal' and
            mode.lower() != 'shifted'):
        logging.critical('Mode must be either None or a string.')
        raise TypeError

    inmap = inatlas.conkitmatch
    acc = 0.0
    minval = 0.0
    maxval = 1.0
    ntp = float(inatlas.tp)

    for contact in inmap:
        if contact.true_positive is True:
            if mode is None or mode.lower() == 'raw':
                acc += contact.raw_score
            elif mode.lower() == 'absolute':
                if contact.raw_score < 0.0:
                    acc -= contact.raw_score
                else:
                    acc += contact.raw_score
            elif mode.lower() == 'shifted' or mode.lower() == 'normal':
                acc += contact.raw_score
                if contact.raw_score < minval:
                    minval = contact.raw_score
                if contact.raw_score > maxval and mode.lower() == 'normal':
                    maxval = contact.raw_score

    if mode.lower() == 'shifted' or mode.lower() == 'normal':
        acc -= minval*ntp
        if mode.lower() == 'normal':
            acc /= (maxval - minval)

    return acc


def avscore(inatlas, mode=None):
    """
    Return the average score of the True Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param mode: None (raw), or 'absolute', or 'normal', or 'shifted', defaults to None.
    :type mode: str, optional
    :return: Average value
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('First argument of avscore must be a Contact Atlas.')
        raise TypeError
    if (mode is not None and mode.lower() != 'raw' and
            mode.lower() != 'absolute' and mode.lower() != 'normal' and
            mode.lower() != 'shifted'):
        logging.critical('Mode must be either None or a string.')
        raise TypeError

    ntp = float(inatlas.tp)

    av = accscore(inatlas, mode=mode) / ntp

    return av


def n_tps(inatlas):
    """
    Return the number of True Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Number of true positives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    return inatlas.tp


def n_fps(inatlas):
    """
    Return the number of False Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Number of false positives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    return inatlas.fp


def n_tns(inatlas):
    """
    Return the number of True Negative contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Number of true negatives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    return inatlas.tn


def n_fns(inatlas):
    """
    Return the number of False Negative contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Number of false negatives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    return inatlas.fn


def mcc(inatlas):
    """
    Return the number of Matthew's Correlation Coefficient.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: MCC
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    ntp = float(inatlas.tp)
    nfp = float(inatlas.fp)
    nfn = float(inatlas.fn)
    ntn = float(inatlas.tn)

    mcc = ntp*ntn - nfp*nfn
    denom = (ntp+nfp)
    denom *= (ntp+nfn)
    denom *= (ntn+nfp)
    denom *= (ntn+nfn)
    mcc /= sqrt(denom)

    return mcc


def precision(inatlas):
    """
    Return the precision of the contact map.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Precision, TP/(TP+FP)
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    ntp = float(inatlas.tp)
    nfp = float(inatlas.fp)

    p = ntp / (ntp + nfp)

    return p


def coverage(inatlas):
    """
    Return the coverage, also known as recall, of the contact map.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Coverage, TP/(TP+FN)
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    ntp = float(inatlas.tp)
    nfn = float(inatlas.fn)

    c = ntp / (ntp + nfn)

    return c


def jaccard(inatlas):
    """
    Return the Jaccard Index of a given matched map. Jaccard = TP / (TP+FP+FN).

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :return: Jaccard Index
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    ntp = float(inatlas.tp)
    nfp = float(inatlas.fp)
    nfn = float(inatlas.fn)

    jacc = ntp / (ntp + nfp + nfn)

    return jacc


def list_scores(inatlas, tag=None):
    """
    Return a list containing all the results in order.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param tag: None ('psicov'), or 'ccmpred', or 'deepmetapsicov', defaults to None.
    :type tag: str, optional
    :return: Results list.
    :rtype: list

    """
    values = []
    psicovmodes = ['normal', 'absolute', 'shifted']

    values.append(str(inatlas.conpred_raw.ncontacts))
    values.append(str(inatlas.conpred.ncontacts))

    acc = accscore(inatlas)
    values.append(str(acc))
    values.append(str(acc/inatlas.tp))
    if tag == 'psicov':
        for m in psicovmodes:
            acc = accscore(inatlas, mode=m)
            values.append(str(acc))
            values.append(str(acc/inatlas.tp))

    values.append(str(inatlas.tp))
    values.append(str(precision(inatlas)))
    values.append(str(coverage(inatlas)))
    values.append(str(mcc(inatlas)))
    values.append(str(jaccard(inatlas)))

    return values
