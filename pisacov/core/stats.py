"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import numpy as np
import copy

def tpr_vs_fpr(scores, against, noneisfalse=True):
    """
    Produce True Positive rate and False Positive rate data by comparing results against standard (ROCs).

    :param scores: Numeric scores.
    :type scores: list [float]
    :param against: Standard for scores to be compared against.
    :type against: list [bool]
    :param noneisfalse: If True, when PISA score given is None, it is considered False, otherwise it is ignored, defaults to True.
    :param noneisfalse: bool, optional

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
        area = 0.5
        return fpr, tpr, area

    numagainst = []
    for b in against:
        if b is True:
            numagainst.append(1)
        elif b is False:
            numagainst.append(0)
        elif b is None:
            numagainst.append(-1)

    scores, against = zip(*sorted(zip(scores, numagainst), reverse=True))

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
                if against[s] == 1:
                    FN += 1
                elif against[s] == 0:
                    TN += 1
                elif against[s] == -1 and noneisfalse is True:
                    TN += 1
                else:
                    pass
            else:
                if against[s] == 1:
                    TP += 1
                elif against[s] == 0:
                    FP += 1
                elif against[s] == -1 and noneisfalse is True:
                    FP += 1
                else:
                    pass
        if (FP+TN) == 0 or (TP+FN) == 0:
            fpr.append(None)
            tpr.append(None)
        else:
            fpr.append(FP/(FP+TN))
            tpr.append(TP/(TP+FN))
            p = -1
            cont = False
            while True:
                if fpr[p-1] is None:
                    p -= 1
                else:
                    break
                if p-1 == -len(tpr):
                    cont = True
                    break
            if cont:
                continue
            else:
                area += 0.5*(tpr[-1]+tpr[p-1])*(fpr[-1]-fpr[p-1])

    return fpr, tpr, area


def hits_vs_total(scores, against, noneisfalse=True):
    """
    Produce True Positive hits and Total hits and misses data by comparing results against standard (TOCs).

    :param scores: Numeric scores.
    :type scores: list [float]
    :param against: Standard for scores to be compared against.
    :type against: list [bool]
    :param noneisfalse: If True, when PISA score given is None, it is considered False, otherwise it is ignored, defaults to True.
    :param noneisfalse: bool, optional

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

    numagainst = []
    for b in against:
        if b is True:
            numagainst.append(1)
        elif b is False:
            numagainst.append(0)
        elif b is None:
            numagainst.append(-1)

    scores, against = zip(*sorted(zip(scores, numagainst), reverse=True))

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
                if against[s] == 1:
                    FN += 1
                elif against[s] == 0:
                    TN += 1
                elif against[s] == -1 and noneisfalse is True:
                    TN += 1
                else:
                    pass
            else:
                if against[s] == 1:
                    TP += 1
                elif against[s] == 0:
                    FP += 1
                elif against[s] == -1 and noneisfalse is True:
                    FP += 1
                else:
                    pass
        if (FP+TN) == 0 or (TP+FN) == 0:
            ntot.append(None)
            hits.append(None)
        else:
            ntot.append(FP+TP)
            hits.append(TP)
            # p = -1
            # cont = False
            # while True:
            #     if fpr[p-1] is None:
            #         p -= 1
            #     else:
            #         break
            #     if p-1 == -len(tpr):
            #         cont = True
            #         break
            # if cont:
            #     continue
            # else:
            #     area += 0.5*(tpr[-1]+tpr[p-1])*(fpr[-1]-fpr[p-1])

    return ntot, hits  # , area


def bezier_parametrization(data, scores, npoints=101, convex=True, emp_tangent=False):
    """
    Calculate a Bézier curve that fits the ROC curve input data.

    :param data: False and True Positive Rate data.
    :type data: list [list [float], list [float]] or list [set [float], set [float]]
    :param scores: Empirical scores.
    :type scores: list [float]
    :param npoints: Number of regular interval points of brezier curve returned, defaults to 101.
    :type npoints: int, optional
    :param convex: Is the curve convex? (i.e. is the slope decreasing for all the curve?), defaults to True.
    :type convex: bool, optional
    :param emp_tangent: Use empirical tangent at (1,1), instead of ROC ideal, defaults to False.
    :type emp_tangent: bool, optional

    :return results: A dictionary with: t parameter, Bézier curve, first and second derivatives, likelihood ratio, probability, scores, lambda value, Younden index, curvature and area.
    :rtype results: dict [list [float], list [list [float]], float]
    """
    # Fierz W (2018) PLoS ONE 13(2): e0192420. https://doi.org/10.1371/journal.pone.0192420
    # Fierz W. MethodsX 7 (2020) 100915. https://doi.org/10.1016/j.mex.2020.100915
    # Hermes, (2017). Journal of Open Source Software, 2(16), 267, https://doi.org/10.21105/joss.00267
    # Yang SN, Huang ML. A New Shape Control and Classification for Cubic Bezier Curves. N. M. Thalmann et al. (eds.), Communicating with Virtual Worlds. © Springer-Verlag Tokyo 1993
    # Zhang, Z, Chen, M, Zhang, X, Wang, Z. IEEE (2009), pp 1-4, https://doi.org/10.1109/CISE.2009.5366218.
    # Cubic Bézier curves: B(t) = (1-t)^3*P_0 + 3*t*(1-t)^2*P_1 + 3*t^2*(1-t))P_2 + t^3*P_3
    # Endpoints: P_0 = (0.0, 0.0), P_3 = (1.0, 1.0)
    # Shape control point: B(t) = S(u) is the data point closest to the top-left vertex (0.0, 1.0) (if convex). Equivalent: Max Youden index (J = Se+Sp-1 = TPR - FPR)
    # Parameter t of the Bézier curve: t = (x+y)/2, with x, y the experimental data. For S(u): u = (x(u) + y(u))/2.
    # End tangents: T_0 def= T_1 = (0.0, 1.0), T_3 = (1.0, 0.0) if convex, else T_0 = (1.0, 0.0), T_3 = (0.0, 1.0).


    P = [(data[0][0], data[1][0]), "P1", "P2", (data[0][-1], data[1][-1])]
    txy = []
    j = 0.0
    for n in range(len(data[0])):
        if P[0] == (0.0, 0.0) and P[3] == (1.0, 1.0):
            x = data[0][n]
            y = data[1][n]
        else:
            x = (data[0][n]-P[0][0])/(P[3][0]-P[0][0])
            y = (data[1][n]-P[0][1])/(P[3][1]-P[0][1])
        txy.append((x+y)/2.0)
        jn = np.abs(y - x)
        if jn > j and convex:
            u = n
            j = jn
        elif jn > j and not convex:
            u = n
            j = jn

    if convex:
        T = [(0.0, 1.0), (0.0, 1.0), "T2", (1.0, 0.0), "T4", "T5" ]
    else:
        T = [(1.0, 0.0), (1.0, 1.0), "T2", (0.0, 1.0), "T4", "T5"]
    if emp_tangent:
        x = data[0][-1]-data[0][-2]
        y = data[1][-1]-data[1][-2]
        norm = np.sqrt(x*x+y*y)
        T[3] = [x/norm, y/norm]

    a = (1.0-txy[u])*(1.0-txy[u])*(1.0+2.0*txy[u])
    b = txy[u]*txy[u]*(3.0-2.0*txy[u])
    D = (T[0][0]*T[3][1])-(T[3][0]*T[0][1])

    tu = txy[u]
    xu = data[0][u]
    yu = data[1][u]

    alpha = (xu-a*P[0][0]-b*P[3][0])*T[3][1] - (yu-a*P[0][1]-b*P[3][1])*T[3][0]
    alpha = alpha / (tu*(1.0-tu)*(1.0-tu)*D)

    beta = (xu-a*P[0][0]-b*P[3][0])*T[0][1] - (yu-a*P[0][1]-b*P[3][1])*T[0][0]
    beta = beta / (tu*tu*(1.0-tu)*D)

    P[1] = ((alpha/3.0)*T[0][0]+P[0][0], (alpha/3.0)*T[0][1]+P[0][1])
    P[2] = ((-beta/3.0)*T[3][0]+P[3][0], (-beta/3.0)*T[3][1]+P[3][1])

    B = [[], []] # BÉZIER CURVE'S POINTS
    V = [[], []] # BÉZIER CURVE'S VELOCITY (1st DERIVATIVE)
    A = [[], []] # BÉZIER CURVE'S ACCELERATION (2nd DERIVATIVE)
    J = [[], []] # BÉZIER CURVE'S JERK (3rd DERIVATIVE)
    K = []       # BÉZIER CURVE'S CURVATURE
    c = [0, 0, 0, 0]
    LR = [] # SLOPE OF TANGENT VECTORS (LIKELIHOOD RATE)
    l = [] # 1 / 1 + LR
    Y = [] # YOUDEN INDICES
    P = [] # PROBABILITY
    area = 0.0
    tparam=[]
    jx = -6.0*P[0][0] + 18.0*P[1][0] - 18.0*P[2][0] + 6.0*P[3][0]
    jy = -6.0*P[0][1] + 18.0*P[1][1] - 18.0*P[2][1] + 6.0*P[3][1]
    for tpoints in range(npoints):
        t = float(tpoints) / float(npoints-1)
        tparam.append(t)
        c[0] = (1.0-t)*(1.0-t)*(1.0-t)
        c[1] = 3.0*t*(1.0-t)*(1.0-t)
        c[2] = 3.0*t*t*(1.0-t)
        c[3] = t*t*t
        Bx = 0.0
        By = 0.0
        for n in range(len(c)):
            Bx += c[n]*P[n][0]
            By += c[n]*P[n][1]

        B[0].append(Bx)
        B[1].append(By)

        c[0] = -3.0*(1.0-t)*(1.0-t)
        c[1] = 3.0*(1.0-t)*(1.0-3.0*t)
        c[2] = 3.0*t*(2.0-3.0*t)
        c[3] = 3.0*t*t

        Bx = 0.0
        By = 0.0
        for n in range(len(c)):
            Bx += c[n]*P[n][0]
            By += c[n]*P[n][1]
        V[0].append(Bx)
        V[1].append(By)

        c[0] = 6.0*(1.0-t)
        c[1] = -12.0 + 18.0*t
        c[2] = 6.0 - 18.0*t
        c[3] = 6.0*t

        Bx = 0.0
        By = 0.0
        for n in range(len(c)):
            Bx += c[n]*P[n][0]
            By += c[n]*P[n][1]
        A[0].append(Bx)
        A[1].append(By)

        J[0].append(jx)
        J[1].append(jy)

        num = V[0][-1]*A[1][-1] - V[1][-1]*A[0][-1]
        den = (V[0][-1]*V[0][-1] + V[1][-1]*V[1][-1])**1.5

        if den == 0.0:
            K.append(float('inf'))
        else:
            K.append(num/den)

        if V[0][-1] == 0.0:
            LR.append(float('inf'))
            l.append(0.0)
            P.append(1.0)
        else:
            LR.append(V[1][-1]/V[0][-1])
            l.append(1.0/(1.0+LR[-1]))
            P.append(l[-1]*LR[-1])
        Y.append(2.0*(l[-1]*By[-1]+(1.0-l)*(1.0-Bx))-1.0)

        if tpoints > 0:
            area += 0.5*(B[1][-1]+B[1][-2])*(B[0][-1]-B[0][-2])

    SCt = list(np.interp(tparam, txy, scores))

    results = {}
    results["param"] = tparam
    results["bezier"] = [B[0], B[1]]
    results["bezier_der1"] = [[V[0], V[1]]]
    results["bezier_der2"] = [[A[0], A[1]]]
    results["bezier_der3"] = [[J[0], J[1]]]
    results["LR"] = LR
    results["curvature"] = K
    results["Youden"] = Y
    results["lambda"] = l
    results["area"] = area
    results["probability"] = P
    results["scores"] = SCt

    return results


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
    wref, set2 = zip(*sorted(zip(setref, set2), reverse=True))

    set1 = np.asarray(set1)
    set2 = np.asarray(set2)

    # corr = np.cov(array1, array2)
    correl = np.corrcoef(set1, set2)

    return correl
