"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import _conf_ops as pco
from pisacov.core import contacts as pcc
from pisacov.core import _psicov_modes as PSICOV_modes

import logging
from numpy import sqrt


def _scorenames(crop=False):
    """Return scorenames.

    :param crop: Include tag for cropped sequence, defaults to False
    :type crop: bool, optional
    :return: scorenames
    :rtype: dict [str: list [str]]

    """
    croptag = 'fullseq' if crop is False else 'cropseq'
    names = pco._sourcenames()
    shortnames = pco._sourcenames(short=True)
    scorenames = {}
    mainnameraw = ['Nconpred', 'Nconused',
                   'ACCScoreRaw', 'AVScoreRaw',
                   'TP', 'PREC', 'COVER', 'MCC', 'JACCARD']
    mainnamepsi = ['Nconpred', 'Nconused',
                   'ACCScoreRaw', 'AVScoreRaw',
                   'ACCScoreNorm', 'AVScoreNorm',
                   'ACCScoreAbs', 'AVScoreAbs',
                   'ACCScoreShift', 'AVScoreShift',
                   'TP', 'PREC', 'COVER', 'MCC', 'JACCARD']
    for i in range(len(names)):
        if names[i] not in scorenames:
            scorenames[names[i]] = []
        lnames = mainnameraw if names[i] != 'psicov' else mainnamepsi
        for mn in lnames:
            scorenames[names[i]].append(mn + '_' +
                                        croptag + '_' +
                                        shortnames[i])

    return scorenames


def accscore(inatlas, alt=None):
    """Return the cumulative score of the True Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Cumulative value
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('First argument of accscore must be a Contact Atlas.')
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    inmap = inatlas.conkitmatch[alt]
    acc = 0.0
    ntp = float(inatlas.tp[alt])

    for contact in inmap:
        if contact.true_positive is True:
            acc += contact.raw_score

    return acc


def avscore(inatlas, alt=None):
    """Return the average score of the True Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Average value
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('First argument of avscore must be a Contact Atlas.')
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    ntp = float(inatlas.tp[alt])

    av = accscore(inatlas, alt) / ntp

    return av


def n_tps(inatlas, alt=None):
    """Return the number of True Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Number of true positives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical("Argument 'inatlas' must be a Contact Atlas.")
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    return inatlas.tp[alt]


def n_fps(inatlas, alt=None):
    """Return the number of False Positive contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Number of false positives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical("Argument 'inatlas' must be a Contact Atlas.")
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    return inatlas.fp[alt]


def n_tns(inatlas, alt=None):
    """Return the number of True Negative contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Number of true negatives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical("Argument 'inatlas' must be a Contact Atlas.")
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    return inatlas.tn[alt]


def n_fns(inatlas, alt=None):
    """
    Return the number of False Negative contacts.

    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Number of false negatives.
    :rtype: int

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical("Argument 'inatlas' must be a Contact Atlas.")
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    return inatlas.fn[alt]


def mcc(inatlas, alt=None):
    """
    Return the number of Matthew's Correlation Coefficient.

    MCC = (TP*TN-FP*FN) / ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: MCC
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    ntp = float(inatlas.tp[alt])
    nfp = float(inatlas.fp[alt])
    nfn = float(inatlas.fn[alt])
    ntn = float(inatlas.tn[alt])

    mcc = ntp*ntn - nfp*nfn
    denom = (ntp+nfp)
    denom *= (ntp+nfn)
    denom *= (ntn+nfp)
    denom *= (ntn+nfn)
    if denom == 0:
        mcc = 'NaN'
    else:
        mcc /= sqrt(denom)

    return mcc


def precision(inatlas, alt=None):
    """
    Return the precision of the contact map.

    Prec = TP / (TP+FP)
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Precision.
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    ntp = float(inatlas.tp[alt])
    nfp = float(inatlas.fp[alt])

    if (ntp + nfp) == 0:
        p = 'NaN'
    else:
        p = ntp / (ntp + nfp)

    return p


def coverage(inatlas, alt=None):
    """
    Return the coverage of the contact map.

    Also known as True Positive Rate, Recall and Sensitivity.
    Cover = TP / (TP+FN)
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Coverage.
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError
    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    ntp = float(inatlas.tp[alt])
    nfn = float(inatlas.fn[alt])

    if (ntp + nfn) == 0:
        c = 'NaN'
    else:
        c = ntp / (ntp + nfn)

    return c


def jaccard(inatlas, alt=None):
    """
    Return the Jaccard Index of a given matched map.

    Jaccard = TP / (TP+FP+FN).
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Jaccard Index
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    ntp = float(inatlas.tp[alt])
    nfp = float(inatlas.fp[alt])
    nfn = float(inatlas.fn[alt])

    if (ntp + nfp + nfn) == 0:
        jacc = 'NaN'
    else:
        jacc = ntp / (ntp + nfp + nfn)

    return jacc


def accuracy(inatlas, alt=None):
    """
    Return the Accuracy of a given matched map.

    Accuracy = (TP+TN) / (TP+TN+FP+FN).
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Accuracy Index
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    ntp = float(inatlas.tp[alt])
    nfp = float(inatlas.fp[alt])
    nfn = float(inatlas.fn[alt])
    ntn = float(inatlas.tn[alt])

    if (ntp + ntn + nfp + nfn) == 0:
        accindex = 'NaN'
    else:
        accindex = (ntp + ntn) / (ntp + nfp + nfn)

    return accindex


def fn_rate(inatlas, alt=None):
    """
    Return the False Negative Rate of a given matched map.

    FNR = TN / (TN+FN).
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: False Negative Rate.
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    nfn = float(inatlas.fn[alt])
    ntn = float(inatlas.tn[alt])

    if (ntn + nfn) == 0:
        fnr = 'NaN'
    else:
        fnr = nfn / (ntn + nfn)

    return fnr


def specificity(inatlas, alt=None):
    """
    Return the False Negative Rate of a given matched map.

    Also known as True Negative Rate.
    SPC = TN / (TN+FP).
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: Specificity.
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    nfp = float(inatlas.fp[alt])
    ntn = float(inatlas.tn[alt])

    if (ntn + nfp) == 0:
        spc = 'NaN'
    else:
        spc = ntn / (ntn + nfp)

    return spc


def fp_rate(inatlas, alt=None):
    """
    Return the False Positive Rate of a given matched map.

    FPR = 1 - SPC = FP / (TN+FP).
    :param inatlas: Contact atlas.
    :type inatlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param alt: Alternative scoring (abs, shifted, norm), defaults to None.
    :type alt: str, optional
    :return: False Positive Rate.
    :rtype: float

    """
    if isinstance(inatlas, pcc.contact_atlas) is False:
        logging.critical('Argument must be a Contact Atlas.')
        raise TypeError

    if (alt is not None and alt != 'raw' and alt != 'abs'
            and alt != 'shifted' and alt != 'norm'):
        logging.critical("Argument must be one of 'abs', 'shifted' or 'norm'.")
        raise TypeError

    if alt is None:
        alt = 'raw'

    nfp = float(inatlas.fp[alt])
    ntn = float(inatlas.tn[alt])

    if (ntn + nfp) == 0:
        fpr = 'NaN'
    else:
        fpr = nfp / (ntn + nfp)

    return fpr


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
    psicovmodes = PSICOV_modes()

    values.append(str(inatlas.conpred_raw.ncontacts))
    values.append(str(inatlas.conkitmatch['raw'].ncontacts))

    acc = accscore(inatlas)
    values.append(str(acc))
    if inatlas.tp['raw'] == 0:
        values.append(str(0.0))
    else:
        values.append(str(acc/inatlas.tp['raw']))
    if tag == 'psicov':
        for m in psicovmodes:
            acc = accscore(inatlas, alt=m)
            values.append(str(acc))
            if inatlas.tp[m] == 0:
                values.append(str(0.0))
            else:
                values.append(str(acc/inatlas.tp[m]))

    values.append(str(inatlas.tp['raw']))
    values.append(str(precision(inatlas)))
    values.append(str(coverage(inatlas)))
    values.append(str(mcc(inatlas)))
    values.append(str(jaccard(inatlas)))

    return values
