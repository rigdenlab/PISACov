"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.core import scores as psc
from pisacov.iomod import _conf_ops as pco
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import numpy as np

_cschemes = {
    'roc_grad':
        {'red':
         ((0.0, 0.8392156862745098, 0.8392156862745098),
          (0.2, 0.8392156862745098, 0.9882352941176471),
          (0.3, 0.9882352941176471, 1.0),
          (0.7, 1.0, 0.30196078431372547),
          (0.8, 0.30196078431372547, 0.21568627450980393),
          (1.0, 0.21568627450980393, 0.21568627450980393)),
         'green':
         ((0.0, 0.15294117647058825, 0.15294117647058825),
          (0.2, 0.15294117647058825, 0.5529411764705883),
          (0.3, 0.5529411764705883, 0.8509803921568627),
          (0.7, 0.8509803921568627, 0.6862745098039216),
          (0.8, 0.6862745098039216, 0.49411764705882355),
          (1.0, 0.49411764705882355, 0.49411764705882355)),
         'blue':
         ((0.0, 0.1568627450980392, 0.1568627450980392),
          (0.2, 0.1568627450980392, 0.3843137254901961),
          (0.3, 0.3843137254901961, 0.1843137254901961),
          (0.7, 0.1843137254901961, 0.2901960784313726),
          (0.8, 0.2901960784313726, 0.7215686274509804),
          (1.0, 0.7215686274509804, 0.7215686274509804))},
    'roc_grad_alt':
        {'red':
         ((0.0, 0.5, 0.5),
          (0.5, 1.0, 1.0),
          (0.7, 1.0, 0.0),
          (0.8, 0.0, 0.0),
          (1.0, 0.0, 0.0)),
         'green':
         ((0.0, 0.0, 0.0),
          (0.5, 0.0, 0.0),
          (0.7, 0.91, 0.5),
          (0.8, 1.0, 0.0),
          (1.0, 0.0, 0.0)),
         'blue':
         ((0.0, 0.0, 0.0),
          (0.5, 0.0, 0.0),
          (0.7, 0.0, 0.3),
          (0.8, 0.0, 0.5),
          (1.0, 1.0, 1.0))},
    'correl_grad':
        {'red':
         ((0.0, 1.0, 1.0),
          (0.25, 0.5, 0.5),
          (0.5, 0.0, 0.0),
          (0.75, 0.0, 0.0),
          (1.0, 0.0, 0.0)),
         'green':
         ((0.0, 1.0, 1.0),
          (0.25, 0.0, 0.0),
          (0.5, 0.0, 0.0),
          (0.75, 0.0, 0.0),
          (1.0, 1.0, 1.0)),
         'blue':
         ((0.0, 0.0, 0.0),
          (0.25, 0.0, 0.0),
          (0.5, 0.0, 0.0),
          (0.75, 0.5, 0.5),
          (1.0, 1.0, 1.0))}
    }

def _gen_cmap(scheme='rocs'):
    # https://github.com/frankligy/scTriangulate/blob/main/image/colors_module/README.md
    if scheme == 'rocs':
        name = 'roc_grad'
    elif scheme == 'rocs2':
        name = 'roc_grad_alt'
    elif scheme == 'correlations':
        name = 'correl_grad'
    elif scheme == 'RdBu' or scheme == 'Pastel1' or scheme == 'Dark2':
        return scheme

    cdict = _cschemes[name]

    return mpl.colors.LinearSegmentedColormap(name, segmentdata=cdict)

class colour_scheme:
    # https://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
    def __init__(self, scheme):
        self.name = scheme
        self.cmap = _gen_cmap(scheme)
        if scheme == 'rocs':
            self.norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
        elif scheme == 'correlations' or scheme == 'RdBu':
            self.norm = mpl.colors.Normalize(vmin=-1.0, vmax=1.0)
        self.scalarMap = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val, pos=None):
        """
        Get RGBA colour coordinates from colour map.

        :param val: Input value.
        :type val: float
        :param pos: If two values correspond to coordinate in the colour scheme, use 'h' or 'l' to retrieve higher or lower value, respectively, defaults to 'l'.
        :type pos: str, optional

        :return: RGBA coordinates.
        :rtype: set[float]

        """
        if pos == 'h':
            plus = 0.05
        elif pos == 'l':
            plus = -0.05
        else:
            plus = 0.0

        return self.scalarMap.to_rgba(val+plus)

color_scheme=colour_scheme

def _get_dpi(ptype):
    """
    Return the resolution of the figure in dpi units.

    :param ptype: Figure's file format. 'png' or 'eps' or 'svg', defaults to 'png'.
    :type ptype: str, optional

    :return: Resolution in dpi units.
    :rtype: int

    """
    if ptype == 'png':
        r = 141
    elif ptype == 'eps' or ptype == 'svg':
        r = 1200
    else:
        logging.warning('Unrecognised plot_type in _get_dpi. Using default PNG.')
        r = 141

    return r

def _set_dpi(ptype='png', seaborn=False):
    """
    Creates figure and axes according to plot type.

    :param ptype: Figure's file format. 'png' or 'eps' or 'svg', defaults to 'png'.
    :type ptype: str, optional

    :return: Figure and axis (Matplotlib)
    :rtype: :class:`~matplotlib.figure.Figure`, :class:`~matplotlib.axes.Axes`

    """
    if ptype == 'png' or ptype == 'eps' or ptype == 'svg':
        fig, ax = plt.subplots(dpi=_get_dpi(ptype))

    else:
        logging.warning('Unrecognised plot_type in plot_rocs. Using default PNG.')
        fig, ax = plt.subplots(dpi=_get_dpi('png'))

    return fig, ax

def plot_matched_map(input_atlas, outpath, mode='raw', plot_type='png', xL=None):
    """Plot matched contact map.

    :param input_atlas: Atlas containing conctact matched contact maps.
    :type input_atlas: :class:`~pisacov.core.contacts.contact_atlas`
    :param outpath: Path to output file.
    :type outpath: str
    :param mode: Mode, if any, defaults to 'raw'.
    :type mode: str, optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image or 'svg' vector image, defaults to 'png'.
    :type plot_type: str, optional
    :param xL: Number of contacts plotted as a function of L, defaults to None (all contacts).
    :type xL: int or float, optional

    """
    if xL is not None:
        nc = round(input_atlas.sequence.length()*xL)
    else:
        nc = len(input_atlas.conkitmatch[mode])

    fpx = []
    fpy = []
    tpx = []
    tpy = []
    fnx = []
    fny = []
    # TO DO: Adapt title to non-crystal dimer inputs
    title = (input_atlas.name + ', ' + 'Interface ' + input_atlas.interface.name +
             ', Chains ' + input_atlas.interface.chains[0].crystal_id
             + input_atlas.interface.chains[1].crystal_id +
             ', ' + input_atlas.conpred_source)
    if input_atlas.conpred_source == 'psicov':
        title += ' (' + mode + ')'

    n = 0
    for contact in input_atlas.conkitmatch[mode]:
        c1 = contact.id[0]
        c2 = contact.id[1]
        if contact.true_positive is True and n < nc:
            n += 1
            tpx.append(c1)
            tpy.append(c2)
        elif contact.false_positive is True and n < nc:
            n += 1
            fpx.append(c1)
            fpy.append(c2)
        else:
            fnx.append(c1)
            fny.append(c2)

    if plot_type == 'agr':
        if len(tpx) > 0:
            tpx, tpy = (list(e) for e in zip(*sorted(zip(tpx, tpy))))
        if len(fpx) > 0:
            fpx, fpy = (list(e) for e in zip(*sorted(zip(fpx, fpy))))
        if len(fnx) > 0:
            fnx, fny = (list(e) for e in zip(*sorted(zip(fnx, fny))))
        xdat = [tpx, fpx, fnx]
        ydat = [tpy, fpy, fny]
        with open(outpath, 'w') as fout:
            for p in range(3):
                fout.write('@target G0.S' + str(p) + '\n@type xy\n')
                for n in range(len(xdat[p])):
                    fout.write(str(xdat[p][n]) + '  ' + str(ydat[p][n]) + '\n')
                fout.write('&\n')
        return
    else:
        fig, ax = _set_dpi(ptype=plot_type)

    ax.set_title(title, y=1.08)

    vmin = 1
    vmax = input_atlas.sequence.length()
    ax.axis([vmin, vmax, vmin, vmax])
    ax.set_xlim(vmin - 0.5, vmax + 0.5)
    ax.set_ylim(vmin - 0.5, vmax + 0.5)

    ax.set_xlabel('Residues from Chain 1')
    ax.set_ylabel('Residues from Chain 2')

    s = ((ax.get_window_extent().width / (vmax-vmin+1.) * 50./fig.dpi) ** 2)
    ax.scatter(tpx, tpy, s=s, marker='o', linewidth=0, c='k', label='Matched (TP)')
    ax.scatter(fpx, fpy, s=s, marker='o', linewidth=0, c='r', label='Unmatched (FP)')
    ax.scatter(fnx, fny, s=s, marker='o', linewidth=0, c='lightgrey', label='Structure (FN)')

    ax.legend(numpoints=1,
              fontsize=10,
              bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
              loc=3,
              ncol=3,
              mode="expand",
              borderaxespad=0.0)
    for n in range(len(ax.legend_.legendHandles)):
        ax.legend_.legendHandles[n]._sizes = [30]
    fig.tight_layout()
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, facecolor='white', transparent=False)
        plt.close(fig)

def plot_rocs(data, outpath, areas_for_color=None, plot_type='png', roc_type='tprvsfpr', norand=False):
    """
    Plot ROC curves.

    :param data: X, Y values of ROCs.
    :type data: dict[list[float], list[float]]
    :param outpath: Output filepath.
    :type outpath: str
    :param areas_for_color: Colour-index assigned to each label, defaults to None.
    :type areas_for_color: dict[float], optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
    :type plot_type: str, optional
    :param roc_type: Type of ROC curve, tprvsfpr' only, defaults to 'tprvsfpr'.
    :type roc_type: str, optional
    :param norand: If True, it omits those curves with 0.3 < area < 0.7, defaults to False.
    :type norand: bool, optional

    """
    if not isinstance(areas_for_color, dict):
        logging.critical('Input must be a dictionary.')
        raise TypeError

    if roc_type == 'tprvsfpr':
        xaxis = 'False Positive Rate (FPR)'
        yaxis = 'True Positive Rate (TPR)'
        title = ('Receiver operating characteristic (ROC) curves: TPR vs. FPR.')
    else:
        logging.warning('Unrecognised roc_type in plot_rocs. Using default TPR vs. FPR.')
        xaxis = 'False Positive Rate (FPR)'
        yaxis = 'True Positive Rate (TPR)'
        title = ('Receiver operating characteristic (ROC) curves: TPR vs. FPR.')

    if plot_type == 'agr':
        # TO DO
        return
    else:
        fig, ax = _set_dpi(ptype=plot_type)
        cmap = colour_scheme('rocs')
        ax.set_box_aspect(1)

    if norand:
        data2 = dict()
        for key in data:
            if areas_for_color[key] <= 0.3 or areas_for_color[key] >= 0.7:
                data2[key] = data[key]
    else:
        data2 = data

    ax.plot([0, 1], [0, 1], linestyle=":", color="dimgray")
    for key in data2:
        if areas_for_color is None:
            clr = 'k'
        else:
            clr = cmap.get_rgb(areas_for_color[key])

        ax.plot(data2[key][0], data2[key][1], linestyle="-", label=key, color=clr)

    ax.set_title(title, y=1.08)

    vmin = 0
    vmax = 1
    ax.axis([vmin, vmax, vmin, vmax])
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)

    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)

    #ax.legend();
    fig.tight_layout()
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, facecolor='white', transparent=False)
        plt.close(fig)


def plot_roc_parametrization(data, bezier, outpath, datatag,
                             area_for_color=None, plot_type='png',
                             roc_type='tprvsfpr', scores=None):
    """
    Plot ROC parametrized curve and probability values.

    :param data: X, Y values of ROCs (FPR, TPR).
    :type data: list[list[float], list[float]]
    :param bezier: Dictionary containing data from Bézier curve parametrization as generated by :func:`core.stats.bezier_parametrization`.
    :type bezier: dict [list [float], list [list [float]], float]
    :param outpath: Output filepath.
    :type outpath: str
    :param datatag: Score's name.
    :type datatag: str
    :param area_for_color: Colour-index assigned to each label, defaults to None.
    :type area_for_color: float, optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
    :type plot_type: str, optional
    :param roc_type: Type of curve, 'tprvsfpr' or 'probvsscore', defaults to 'tprvsfpr'.
    :type roc_type: str, optional
    :param scores: Empirical scores, must have the same length as FPR and TPR (only needed if roc_type=='probvsscore'), defaults to None.
    :type scores: list[float], optional

    """
    if not isinstance(area_for_color, float):
        logging.critical('Input must be a float.')
        raise TypeError

    if roc_type not in ['tprvsfpr', 'probvsscore']:
        logging.warning('Unrecognised roc_type in plot_roc_parametrization. Using default TPR vs. FPR.')
        roc_type="tprvsfpr"

    if roc_type == 'probvsscore' and scores is None:
        logging.critical("The roc_type='probvsscore' requires a list of empirical scores.")
        raise ValueError

    if plot_type == 'agr':
        # TO DO
        return
    else:
        fig, ax = _set_dpi(ptype=plot_type)
        cmap = colour_scheme('rocs')

    if roc_type=="tprvsfpr":
        ax.plot([0, 1], [0, 1], linestyle=":", color="dimgray")

    if area_for_color is None:
        clr = 'k'
    else:
        clr = cmap.get_rgb(area_for_color)

    p_emp = []
    for n in range(len(data[0])):
        t_emp = (data[0][n]+data[1][n])/2.0
        p_emp.append(np.interp(t_emp, bezier["param"], bezier["probability"]))

    ac = bezier["anticorrelated"]
    if roc_type == 'probvsscore':
        xaxis = datatag
        yaxis = 'True Positive Rate (TPR), Probability (P)'
        if not ac:
            title = ('Probability and TPR vs score for ' + datatag + ".")
        else:
            title = ('Probability (anticorrelated function) and TPR vs score for ' + datatag + ".")
        xmin = scores[-1]
        xmax = scores[0]
        xcont = bezier["scores"]
        xdisc = scores
    elif roc_type == 'tprvsfpr':
        ax.set_box_aspect(1)
        xaxis = 'False Positive Rate (FPR)'
        yaxis = 'True Positive Rate (TPR), Probability (P)'
        if not ac:
            title = ('Bézier approximation to ROC curve. ' + datatag + ".")
        else:
            title = ('Bézier approximation to ROC curve (anticorrelated prob.). ' + datatag + ".")
        xmin = -0.03
        xmax = 1.03
        xcont = bezier["bezier"][0]
        xdisc = data[0]
    ax.plot(xdisc, data[1], linestyle=" ", marker=".", label="TPR (emp)", color=clr)
    ax.plot(xcont, bezier['bezier'][1], linestyle="--", label="TPR (t)", color=clr)
    ax.plot(xdisc, p_emp, linestyle=" ", marker = ".", markerfacecolor='none', label="P (emp)", markeredgecolor="dimgray")
    ax.plot(xcont, bezier['probability'], linestyle="-", label="P (t)", color="dimgray")

    ax.set_title(title, y=1.08)

    ymin = -0.03
    ymax = 1.03
    ax.axis([xmin, xmax, ymin, ymax])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)

    ax.legend();
    fig.tight_layout()
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, facecolor='white', transparent=False)
        plt.close(fig)


def plot_toc(data, datatag, outpath, area_for_color=None, plot_type='png'):
    """
    Plot TOC curve.

    :param data: X, Y values of TOC.
    :type data: list[list[float], list[float]]
    :param datatag: Data description of TOC (score's name).
    :type datatag: str
    :param outpath: Output filepath.
    :type outpath: str
    :param area_for_color: Colour-index assigned to curve, defaults to None.
    :type area_for_color: float, optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
    :type plot_type: str, optional

    """
    if not isinstance(area_for_color, float):
        logging.critical('Input must be a float.')
        raise TypeError

    xaxis = 'Hits + False Alarms (FP + TP)'
    yaxis = 'Hits (TP)'
    title = ('Total operating characteristic (TOC) curve.')

    if plot_type == 'agr':
        # TO DO
        return
    else:
        fig, ax = _set_dpi(ptype=plot_type)
        cmap = colour_scheme('rocs')

    xmax = max(data[0])
    ymax = max(data[1])
    ax.plot([0, xmax], [0, ymax], linestyle=":", color="dimgray")           # RANDOM
    ax.plot([0, ymax], [0, ymax], linestyle=":", color="seagreen")          # MAX
    ax.plot([xmax-ymax, xmax], [0, ymax], linestyle=":", color="indianred") # MIN

    if area_for_color is None:
        clr = 'k'
    else:
        clr = cmap.get_rgb(area_for_color)

    ax.plot(data[0], data[1], linestyle="-", label=datatag, color=clr)

    ax.set_title(title, y=1.08)

    ax.axis([0, xmax, 0, ymax])
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)

    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)

    ax.legend();
    fig.tight_layout()
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, facecolor='white', transparent=False)
        plt.close(fig)

def plot_correlation_sns(data, outpath, labels=None, plot_type='png',
                     light0=False, clustered=False, show_values=False,
                     areas_for_color=None):
    """Plot correlation heatmap with Seaborn. Add dendrograms optionally.

    :param data: Correlation matrix.
    :type data: :class:`~np.array`
    :param outpath: Output filepath.
    :type outpath: str
    :param labels: Label list, defaults to None.
    :type labels: list[str], optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, or 'svg' vector image, defaults to 'png'.
    :type plot_type: str, optional
    :param light0: Use 'RdBu' colour map if True, custom yellow-black-cyan map if False, defaults to False.
    :type light0: bool
    :param clustered: Cluster scores by similarity and produce dendogram on heatmap, defaults to False.
    :type clustered: bool
    :param show_values: Show correlation values on each cell, defaults to False.
    :type show_values: bool, optional
    :param areas_for_color: Area-related colour-index assigned to each label (only when clustered is True), defaults to None.
    :type areas_for_color: dict[float], optional

    """
    # CHECK https://seaborn.pydata.org/generated/seaborn.clustermap.html
    # CHECK https://seaborn.pydata.org/generated/seaborn.heatmap.html
    import seaborn as sns; sns.set(color_codes=True)
    import pandas as pd
    import copy

    colours = 'RdBu' if light0 is True else 'correlations'
    cmap = colour_scheme(colours)

    title = 'Correlation matrix for PISACov scores.'

    if labels is None:
        labels = list(range(len(data)))
        labels = [str(a) for a in labels]

    # labels_orig = copy.deepcopy(labels)
    # labels_rev = list(labels).reverse()

    if clustered:
        row_coloring = []
        markers1 = pco._sourcenames(short=True)
        #lut1 = dict(zip(set(markers1),
        #                sns.hls_palette(len(set(markers1)), l=0.3, s=0.8)))
        lut1 = dict(zip(set(markers1),
                        sns.color_palette(palette='Dark2',
                                          n_colors=len(set(markers1)))))
        lut1_long = {}
        for label in labels:
            if label == 'PISAscore':
                lut1_long[label] = (1.0, 1.0, 1.0)
            else:
                for m in markers1:
                    if m in label:
                        lut1_long[label] = lut1[m]
        row_coloring.append(pd.DataFrame(labels)[0].map(lut1_long))

        markers2 = ['Nconpred', 'Nconused', 'AccScore', 'AvScore',
                    'TP', 'PREC', 'COVER', 'MCC', 'JACCARD']
        #lut2 = dict(zip(set(markers2),
        #                sns.hls_palette(len(set(markers2)), l=0.5, s=0.8)))
        lut2 = dict(zip(set(markers2),
                        sns.color_palette(palette='Pastel1',
                                          n_colors=len(set(markers2)))))
        lut2_long = {}

        for label in labels:
            if label == 'PISAscore':
                lut2_long[label] = (1.0, 1.0, 1.0)
            else:
                for m in markers2:
                    if m in label:
                        lut2_long[label] = lut2[m]
        row_coloring.append(pd.DataFrame(labels)[0].map(lut2_long))

        if areas_for_color is not None:
            cmap2 = colour_scheme('rocs')
            lut3_long = {}
            lut3_long['PISAscore'] = (1.0, 1.0, 1.0)
            for key in areas_for_color:
                lut3_long[key] = cmap2.get_rgb(val=areas_for_color[key])
            row_coloring.append(pd.DataFrame(labels)[0].map(lut3_long))

        title += '\nClustered by similarity.'

        df = pd.DataFrame(data=data, index=labels, columns=labels)
        cg = sns.clustermap(df, vmin=-1, vmax=1, cmap=cmap.cmap, annot=False,
                            cbar_pos=(1.05, 0, 0.08, 0.8),
                            xticklabels=True, yticklabels=True,
                            linewidths=0.1) #, row_colors=row_coloring)

        #cg.ax_row_dendrogram.set_visible(False)
        #cg.ax_row_dendrogram.set_xlim([0,0])
    else:
        title += '\nSorted by Area under the ROC.'

        df = pd.DataFrame(data=data, index=labels, columns=labels)
        plt.subplots(figsize=(25,25))
        cg = sns.heatmap(df, vmin=-1, vmax=1, annot=False, cmap=cmap.cmap,
                         xticklabels=True, yticklabels=True, linewidths=0.1)

        plt.tight_layout()
    plt.title(title, y=1.40, x=-8)
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        plt.savefig(outpath, format=plot_type, facecolor='white', transparent=False)
        plt.close()

def plot_correlation_heatmap(data, outpath, labels=None, plot_type='png',
                             light0=False, show_values=False):
    """
    Plot correlation heatmap.

    :param data: Correlation matrix.
    :type data: :class:`~np.array`
    :param outpath: Output filepath.
    :type outpath: str
    :param labels: Label list, defaults to None.
    :type labels: list[str], optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, or 'svg' vector image, defaults to 'png'.
    :type plot_type: str, optional
    :param light0: Use 'RdBu' colour map if True, custom yellow-black-cyan map if False, defaults to False.
    :type light0: bool
    :param show_values: Show correlation values in each cell, defaults to False.
    :type show_values: bool

    """
    # CHECK https://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-on-top-of-a-matrix-of-data
    # CHECK https://www.python-graph-gallery.com/405-dendrogram-with-heatmap-and-coloured-leaves
    colours = 'RdBu' if light0 is True else 'correlations'

    title = 'Correlation matrix for PISACov scores.\nSorted by Area under the ROC.'

    if labels is None:
        labels = list(range(len(data)))
        labels = [str(a) for a in labels]

    if plot_type == 'agr':
        # TO DO
        return
    else:
        fig, ax = _set_dpi(ptype=plot_type)
        cmap = colour_scheme(colours)

    ax.imshow(np.fliplr(data), cmap=cmap.cmap, norm=cmap.norm)
    # ax.imshow(data, cmap=cmap.cmap)

    # ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_xticks(np.flip(np.arange(len(labels))),
                  labels=labels, fontsize=4)
    ax.set_yticks(np.arange(len(labels)),
                  labels=labels, fontsize=4)

    # Minor ticks
    ax.set_xticks(np.arange(-.5, len(labels), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)

    #ax.axhline(0.5, color="white", lw=2)
    #ax.axvline(len(labels)-1.5, color="white", lw=2)
    plt.grid(False)
    plt.grid(which='minor', linewidth=0.1, color='w')

    # Remove minor ticks
    ax.tick_params(which='minor', bottom=False, left=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    if show_values is True:
        fdatar = np.fliplr(np.around(data,decimals=2))
        for n in range(len(labels)):
            for m in range(len(labels)):
                v = fdatar[n,m]
                if light0:
                    c = "k" if abs(v) < 0.3 else "w"
                else:
                    c = "k" if abs(v) > 0.8 else "w"
                ax.text(m, n, v,
                        ha="center", va="center", color=c)

    ax.set_title(title, y=1.08)

    cb=fig.colorbar(cmap.scalarMap, ax=ax)
    cb.set_ticks(np.arange(-1.0, 1.01, 0.25))
    cb.outline.set_visible(False)

    fig.tight_layout()
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type)
        plt.close(fig)


def plot_area_histogram(data, outpath, plot_type='png'):
    """
    Plot histogram with scores' areas.

    :param data: Data (areas) and labels (dictionary), added by area.
    :type data: dict[float]
    :param outpath: Output filepath.
    :type outpath: str, optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
    :type plot_type: str, optional

    """
    # https://www.tutorialspoint.com/matplotlib/matplotlib_bar_plot.htm
    # https://stackoverflow.com/questions/28129606/how-to-create-a-matplotlib-bar-chart-with-a-threshold-line

    title = 'Histogram of ROC areas per score.'
    xaxis = "Score"
    yaxis = "ROC area"

    if plot_type == 'agr':
        # TO DO
        return
    else:
        fig, ax = _set_dpi(ptype=plot_type)
        cmap = colour_scheme('rocs')

    clist = [cmap.get_rgb(data[i]) for i in data]

    # Plot histogram data
    ax.bar(list(data.keys()), list(data.values()), color=clist)

    ax.set_xticks(np.arange(len(data)), labels=list(data.keys()), fontsize=4)
    ax.set_ylim([0, 1])

    cb = fig.colorbar(cmap.scalarMap, ax=ax)
    cb.set_ticks(np.arange(0.0, 1.01, 0.1))

    # Add horizontal lines
    # ax.axhline(y=0.5, color = cmap.get_rgb(0.5), linestyle="--")
    ax.axhline(y=0.7, color = cmap.get_rgb(0.75), linestyle="--")
    ax.axhline(y=0.8, color = cmap.get_rgb(1.00), linestyle="--")

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    plt.setp(ax, yticks=np.arange(0.0, 1.01, 0.1))

    ax.set_title(title, y=1.08)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    plt.grid(False)

    fig.tight_layout()
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, facecolor='white', transparent=False)
        plt.close(fig)

####################################################################################################
#"""
#=====================
#3D surface (colormap)
#=====================

#Demonstrates plotting a 3D surface colored with the coolwarm colormap.
#The surface is made opaque by using ``antialiased=False``.

#Also demonstrates using the `.LinearLocator` and custom formatting for the
#z axis tick labels.
#"""

#import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator
#import numpy as np

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
#X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)

# Plot the surface.
#surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
#ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.show()


#############################################################################
#
# .. admonition:: References
#
#    The use of the following functions, methods, classes and modules is shown
#    in this example:
#
#    - `matplotlib.pyplot.subplots`
#    - `matplotlib.axis.Axis.set_major_formatter`
#    - `matplotlib.axis.Axis.set_major_locator`
#    - `matplotlib.ticker.LinearLocator`
#    - `matplotlib.ticker.StrMethodFormatter`
