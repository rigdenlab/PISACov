"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
import numpy as np

_cschemes = {
    'roc_grad':
        {'red':
         ((0.0, 0.5, 0.5),
          (0.5, 1.0, 1.0),
          (0.6, 1.0, 1.0),
          (0.7, 0.6, 0.6),
          (0.8, 0.0, 0.0),
          (1.0, 0.0, 0.0)),
         'green':
         ((0.0, 0.0, 0.0),
          (0.5, 0.0, 0.0),
          (0.6, 0.91, 0.91),
          (0.7, 0.6, 0.6),
          (0.8, 0.5, 0.5),
          (1.0, 0.0, 0.0)),
         'blue':
         ((0.0, 0.0, 0.0),
          (0.5, 0.0, 0.0),
          (0.6, 0.0, 0.0),
          (0.7, 0.0, 0.0),
          (0.8, 0.5, 0.5),
          (1.0, 0.7, 0.7))},
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

def gen_cmap(scheme='rocs'):
    # https://github.com/frankligy/scTriangulate/blob/main/image/colors_module/README.md
    if scheme == 'rocs':
        name = 'roc_grad'
    elif scheme == 'correlations':
        name = 'correl_grad'

    cdict = _cschemes[name]

    return mpl.colors.LinearSegmentedColormap(name, segmentdata=cdict)


class colour_scheme:
    # https://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
    def __init__(self, scheme):
        self.name = scheme
        self.cmap = gen_cmap(scheme)
        if scheme == 'rocs':
            self.norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
        elif scheme == 'correlations':
            self.norm = mpl.colors.Normalize(vmin=-1.0, vmax=1.0)
        self.scalarMap = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val):
        return self.scalarMap.to_rgba(val)


color_scheme=colour_scheme


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

    if plot_type == 'png':
        fig, ax = plt.subplots(dpi=141)
    elif plot_type == 'eps':
        fig, ax = plt.subplots(dpi=1200)
    elif plot_type == 'svg':
        fig, ax = plt.subplots(dpi=1200)
    elif plot_type == 'agr':
        tpx, tpy = (list(e) for e in zip(*sorted(zip(tpx, tpy))))
        fpx, fpy = (list(e) for e in zip(*sorted(zip(fpx, fpy))))
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

    ax.set_title(title, y=1.08)

    vmin = 1
    vmax = input_atlas.sequence.length()
    ax.axis([vmin, vmax, vmin, vmax])
    ax.set_xlim(vmin - 0.5, vmax + 0.5)
    ax.set_ylim(vmin - 0.5, vmax + 0.5)

    ax.set_xlabel('Residues from Chain 1')
    ax.set_ylabel('Residues from Chain 2')

    s = ((ax.get_window_extent().width  / (vmax-vmin+1.) * 50./fig.dpi) ** 2)
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
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, overwrite=True)
        plt.close(fig)


def plot_rocs(data, outpath, areas_for_color=None, plot_type='png', roc_type='tprvsfpr'):
    """
    Plot ROC curves.

    :param data: X, Y values of ROCs.
    :type data: dict[list[float], list[float]]
    :param outpath: Output filepath.
    :type outpath: str, optional
    :param areas_for_color: Colour-index assigned to each label, defaults to None.
    :type areas_for_color: dict[float], optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
    :type plot_type: str, optional
    :param roc_type: Type of ROC curve, tprvsfpr' only, defaults to 'tprvsfpr'
    :type roc_type: str, optional

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

    if plot_type == 'png':
        fig, ax = plt.subplots(dpi=141)
    elif plot_type == 'eps':
        fig, ax = plt.subplots(dpi=1200)
    elif plot_type == 'svg':
        fig, ax = plt.subplots(dpi=1200)
    else:
        logging.warning('Unrecognised plot_type in plot_rocs. Using default PNG.')
        fig, ax = plt.subplots(dpi=141)

    if plot_type == 'png' or plot_type == 'eps':
        cmap = colour_scheme('rocs')

    ax.plot([0, 1], [0, 1], linestyle=":", color="dimgray")
    for key in data:
        if areas_for_color is None:
            clr = 'k'
        else:
            clr = cmap.get_rgb(areas_for_color[key])

        ax.plot(data[key][0], data[key][1], linestyle="-", label=key, color=clr)

    ax.set_title(title, y=1.08)

    vmin = 0
    vmax = 1
    ax.axis([vmin, vmax, vmin, vmax])
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)

    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)

    ax.legend();
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, overwrite=True)
        plt.close(fig)


def plot_toc(data, datatag, outpath, area_for_color=None, plot_type='png'):
    """
    Plot TOC curve.

    :param data: X, Y values of TOC.
    :type data: list[list[float], list[float]]
    :param datatag: Data description of TOC.
    :type datatag: str
    :param outpath: Output filepath.
    :type outpath: str, optional
    :param areas_for_color: Colour-index assigned to each label, defaults to None.
    :type areas_for_color: dict [float], optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, 'svg' vector image, or 'agr' in raw grace format, defaults to 'png'.
    :type plot_type: str, optional

    """
    if not isinstance(area_for_color, float):
        logging.critical('Input must be a float.')
        raise TypeError

    xaxis = 'Hits + False Alarms (FP + TP)'
    yaxis = 'Hits (TP)'
    title = ('Total operating characteristic (TOC) curve.')

    if plot_type == 'png':
        fig, ax = plt.subplots(dpi=141)
    elif plot_type == 'eps':
        fig, ax = plt.subplots(dpi=1200)
    elif plot_type == 'svg':
        fig, ax = plt.subplots(dpi=1200)
    else:
        logging.warning('Unrecognised plot_type in plot_toc. Using default PNG.')
        fig, ax = plt.subplots(dpi=141)

    if plot_type == 'png' or plot_type == 'eps':
        cmap = colour_scheme('rocs')


    xmax = max(data[0])
    ymax = max(data[1])
    ax.plot([0, xmax], [0, ymax], linestyle=":", color="dimgray")
    ax.plot([0, xmax], [0, xmax], linestyle=":", color="seagreen")
    ax.plot([0, ymax], [0, ymax], linestyle=":", color="indianred")


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
    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, overwrite=True)
        plt.close(fig)


def plot_correlation_heatmap(data, outpath, labels=None, plot_type='png'):
    """
    Plot correlation heatmap.

    :param data: Correlation matrix.
    :type data: :class:`~np.array`
    :param outpath: Output filepath.
    :type outpath: str, optional
    :param labels: Label list, defaults to None.
    :type labels: list[str], optional
    :param plot_type: Plot either as a 'png' image, 'eps' vector image, or 'svg' vector image, defaults to 'png'.
    :type plot_type: str, optional

    """


    title = ('Correlation matrix for PISACov scores.')

    if labels is None:
        labels = list(range(len(data)))
        labels = [str(a) for a in labels]

    if plot_type == 'png':
        fig, ax = plt.subplots(dpi=141)
    elif plot_type == 'eps':
        fig, ax = plt.subplots(dpi=1200)
    elif plot_type == 'svg':
        fig, ax = plt.subplots(dpi=1200)
    else:
        logging.warning('Unrecognised plot_type in plot_rocs. Using default PNG.')
        fig, ax = plt.subplots(dpi=141)

    im = ax.imshow(data)

    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for n in range(len(labels)):
        for m in range(len(labels)):
            text = ax.text(m, n, data[n, m],
                           ha="center", va="center", color="w")

    ax.set_title(title, y=1.08)

    if plot_type == 'png' or plot_type == 'eps' or plot_type == 'svg':
        fig.savefig(outpath, format=plot_type, overwrite=True)
        plt.close(fig)


#def area_histogram(data, outpath, labels=None, plot_type='png'):
#https://www.tutorialspoint.com/matplotlib/matplotlib_bar_plot.htm
#https://stackoverflow.com/questions/28129606/how-to-create-a-matplotlib-bar-chart-with-a-threshold-line
