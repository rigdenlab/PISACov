"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import os
import logging

from conkit import io as ckio
from conkit import plot as ckplot

from pisacov.io.conf import HHBLITS_PATH
from pisacov.io.conf import HHBLITS_DATABASE_DIR
from pisacov.io.conf import HHBLITS_DATABASE_NAME
from pisacov.io import paths as ppaths

from matplotlib import pyplot as plt

def runhhblits(spath, hhparam, outdirmsa):
    """
    Run HHSUITE to produce a Multiple Sequence Alignment from a single fasta sequence.

    :param spath: Sequence path.
    :type spath: str
    :param hhparam: HHBLITS sequence alignment parameters.
    :type hhparam: list [str], dim=5
    :param outdirmsa: Output path
    :type outdirmsa: str
    :return: MSA
    :rtype: str

    """
    hhsuite_exec = '"' + ppaths.check_path(HHBLITS_PATH, 'file') + '"'

    msaa3mfile = (os.path.splitext(os.path.basename(spath))[0] +
                  os.extsep + "msa" + os.extsep + "a3m")
    msaa3mpath = os.path.join(outdirmsa, msaa3mfile)

    hhrfile = (os.path.splitext(os.path.basename(spath))[0] +
                  os.extsep + "hhr")
    hhrpath = os.path.join(outdirmsa, hhrfile)

    execstring = (hhsuite_exec + ' -i ' + spath +
                  ' -d ' + os.path.join(ppaths.check_path(HHBLITS_DATABASE_DIR, 'dir'),
                                        HHBLITS_DATABASE_NAME) +
                  ' -n ' + str(hhparam[0]) +
                  ' -e ' + str(hhparam[1]) + ' -diff ' + str(hhparam[2]) +
                  ' -cov ' + str(hhparam[3]) + ' -id ' + str(hhparam[4]) +
                  ' -oa3m ' + msaa3mpath + ' -o ' + hhrpath)
    try:
        os.system(execstring)
    except Exception:
        logging.critical('        An error occurred while executing HHBLITS.')
        raise OSError

    # Convert A3M MSA file to Jones format (DMP standard input format)
    parsedmsa = msafilesgen(msaa3mpath)

    return parsedmsa


def msafilesgen(inpath_a3m):
    """Convert a3m format alignment to "jones" format (.aln) and print msa coverage graph.

    :param inpath_a3m: Multiple Sequence Alignment (.a3m) path
    :type inpath_a3m: str
    :return: ConKit parsed MSA
    :rtype: :obj:`~conkit.core.sequencefile.SequenceFile`

    """
    parsedmsa = ckio.read(inpath_a3m, 'a3m')

    # Convert to 'jones' format
    outpath_aln = os.path.splitext(inpath_a3m)[0] + os.extsep + "aln"
    ckio.write(outpath_aln, 'jones', parsedmsa)

    # Plot Coverage
    try:
        msacoveragepath = (os.path.splitext(inpath_a3m)[0] + os.extsep +
                           "coverage" + os.extsep + "png")
        fig = ckplot.SequenceCoverageFigure(parsedmsa)
        fig.savefig(msacoveragepath, overwrite=True)
        del fig
    except Exception:
        logging.warning('Something went wrong with ConKit ' +
                        'and Coverage Plot was not produced.')

    return parsedmsa
