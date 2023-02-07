"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import numpy as np
import copy

def tpr_vs_fpr(scores, against):
    """
    Produce True Positive rate and False Positive rate data by comparing results against standard (ROCs).

    :param scores: Numeric scores.
    :type scores: list [float]
    :param against: Standard for scores to be compared against.
    :type against: list [bool]

    :return fpr: False Positive Rate data.
    :rtype fpr: set [float]
    :return tpr: True Positive Rate data.
    :rtype tpr: set [float]
    :return area: Area under the ROC curve.
    :rtype area: float

    """
    if len(scores)==0:
        fpr = None
        tpr = None
        area = 0
        return fpr, tpr, area

    scores, against = zip(*sorted(zip(scores, against), reverse=True))
    thr = (list(set(scores)))

    tlist = sorted(thr, reverse=True)
    tpr = []
    fpr = []
    fpr.append(0)
    tpr.append(0)
    area = 0
    for t in tlist:
        FP = 0
        TP = 0
        FN = 0
        TN = 0
        for s in range(len(scores)):
            if scores[s] < t:
                if against[s]:
                    FN += 1
                else:
                    TN += 1
            else:
                if against[s]:
                    TP += 1
                else:
                    FP += 1
        if (FP+TN) == 0 or (TP+FN) == 0:
            fpr.append(None)
            tpr.append(None)
        else:
            fpr.append(FP/(FP+TN))
            tpr.append(TP/(TP+FN))
            p = -1
            while True:
                if fpr[-1] is None:
                    p -= 1
                else:
                    break
            area += 0.5*(tpr[-1]-tpr[p-1])*(fpr[-1]+fpr[p-1])

    return fpr, tpr, area


def hits_vs_total(scores, against):
    """
    Produce True Positive hits and Total hits and misses data by comparing results against standard (TOCs).

    :param scores: Numeric scores.
    :type scores: list [float]
    :param against: Standard for scores to be compared against.
    :type against: list [bool]

    :return ntot: Hits (TPs) + False Alarms (FPs) data.
    :rtype ntot: set [float]
    :return hits: Hits (TPs) data.
    :rtype hits: set [float]
    :return area: Area under the ROC curve.
    :rtype area: float

    """
    if len(scores)==0:
        ntot = None
        hits = None
        # area = 0
        return ntot, hits  # , area

    scores, against = zip(*sorted(zip(scores, against), reverse=True))
    thr = (list(set(scores)))

    tlist = sorted(thr, reverse=True)
    hits = []
    ntot = []
    ntot.append(0)
    hits.append(0)
    # area = 0
    for t in tlist:
        FP = 0
        TP = 0
        FN = 0
        TN = 0
        for s in range(len(scores)):
            if scores[s] < t:
                if against[s]:
                    FN += 1
                else:
                    TN += 1
            else:
                if against[s]:
                    TP += 1
                else:
                    FP += 1
        if (FP+TN) == 0 or (TP+FN) == 0:
            ntot.append(None)
            hits.append(None)
        else:
            ntot.append(FP+TP)
            hits.append(TP)
            # p = -1
            # while True:
            #     if ntot[-1] is None:
            #         p -= 1
            #     else:
            #         break
            # area += 0.5*(tpr[-1]-tpr[p-1])*(fpr[-1]+fpr[p-1])

    return ntot, hits  # , area


def correl_matrix(set1, set2, setref=None):
    """
    Calculate various correlation matrices for the given data vectors.

    :param set1: Data set 1, sort according to this set unless setref given.
    :type set1: list [float or None]
    :param set2: Data set 2.
    :type set2: list [float or None]
    :param setref: Reference data set, sort the others according to it, defaults to None
    :type setref: list [float or None], optional

    :return: Pearson correlation matrix.
    :rtype: :class:`~numpy.matrix`

    """
    if setref is None:
        setref = copy.deepcopy(set1)

    wref, set1 = zip(*sorted(zip(setref, set1), reverse=True))
    wref, set1 = zip(*sorted(zip(setref, set2), reverse=True))

    set1 = np.asarray(set1)
    set2 = np.asarray(set2)

    # corr = np.cov(array1, array2)
    correl = np.corrcoef(set1, set2)

    return correl
