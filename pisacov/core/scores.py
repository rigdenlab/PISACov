"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.iomod import _conf_ops as pco
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
    shortnames = pco._sourcenames(short=True)
    scorenames = {}
    psicovmodes = PSICOV_modes(raw=True)


    nconpredname = 'Nconpred_' + croptag
    mainname = ['Nconused_'+croptag,
                'AccScore_'+croptag,
                'AvScore_'+croptag,
                'TP'+croptag,
                'PREC'+croptag,
                'COVER'+croptag,
                'MCC'+croptag,
                'JACCARD'+croptag]

    for s in shortnames:
        if s == 'psicov':
            modes = psicovmodes
        else:
            modes = ['raw']

        for m in modes:
            if m =='raw':
                scorenames[s] = []
                scorenames[s].append(nconpredname+'_'+s)
            for n in mainname:
                scorenames[s].append(n+'_'+m+s)

    return scorenames

def _score_properties(input_score, descriptor):
    """Return score markers according to descriptor.

    :param input_score: Score(s) from which markers are generated.
    :type input_score: str, list[str]
    :param descriptor: Either 'scores' or 'sources'.
    :type descriptor: str

    :raises TypeError: If input_score is not a string or a list of strings, or descriptor is not a string.
    :raises ValueError: If descriptor is not one of the two valid strings.

    :return: The marker(s) for the given score(s) and descriptor.
    :rtype: str, list[str]

    """
    if input_score.isinstance(str):
        i = [input_score]
    elif input_score.isinstance(list):
        for n in range(len(input_score)):
            if input_score[n].isinstance(str):
                pass
            else:
                msg = 'Input must be either a string or a list of strings.'
                logging.critical(msg)
                raise TypeError
        i = input_score
    else:
        msg = 'Input must be either a string or a list of strings.'
        logging.critical(msg)
        raise TypeError

    if descriptor.isinstance(str):
        if descriptor.lower() != 'source' and descriptor.lower() != 'score':
            msg = "Descriptor must be either 'score' or 'source'."
            logging.critical(msg)
            raise ValueError
    else:
        msg = 'Descriptor must a string.'
        logging.critical(msg)
        raise TypeError

    r =[]
    d = descriptor.lower()
    if d == 'scores':
        s = ['PISA', 'Nconpred', 'Nconused', 'AccScore', 'AvScore',
             'TP', 'PREC', 'COVER', 'MCC', 'JACCARD']
    elif d == 'sources':
        s = pco._sourcenames(short=True)

    for n in len(i):
        for element in s:
            if element in i[n]:
                r.append(element)

    if input_score.isinstance(str):
        return r[0]
    else:
        return r


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
    psicovmodes = PSICOV_modes(raw=True)

    values.append(str(inatlas.conpred_raw.ncontacts))

    if tag == 'psicov':
        modes = psicovmodes
    else:
        modes = ['raw']

    for m in modes:
        values.append(str(inatlas.conkitmatch[m].ncontacts))
        acc = accscore(inatlas, alt=m)
        values.append(str(acc))
        if inatlas.tp[m] == 0:
            values.append(str(0.0))
        else:
            values.append(str(acc/inatlas.tp[m]))

        values.append(str(inatlas.tp[m]))
        values.append(str(precision(inatlas, alt=m)))
        values.append(str(coverage(inatlas, alt=m)))
        values.append(str(mcc(inatlas, alt=m)))
        values.append(str(jaccard(inatlas, alt=m)))

    return values
