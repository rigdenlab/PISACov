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

from pisacov.io.conf import HHSUITE_PATH
from pisacov.io.conf import HHBLITS_DATABASE_DIR
from pisacov.io.conf import HHBLITS_DATABASE_NAME
from pisacov import io as pio

def runhhblits(spath, hhparam, outdirmsa):
    """
    Run HHSUITE to produce a Multiple Sequence Alignment from a single fasta sequence.

    :param spath: Sequence path.
    :type spath: str
    :param hhparam: HHBLITS sequence alignment parameters.
    :type hhparam: list [str], dim=5
    :param outdirmsa: Output path
    :type outdirmsa: str
    :return: MSA filepath
    :rtype: str

    """
    hhsuite_exec='"'+pio.check_path(HHSUITE_PATH, 'file')+'"'

    msaa3mfile= os.path.splitext(os.path.basename(spath))[0] +".msa.a3m"
    msaa3mpath = os.path.join(outdirmsa, msaa3mfile)

    try:
        os.system(hhsuite_exec + ' -i '+ spath +
                  ' -d ' + pio.check_path(HHBLITS_DATABASE_DIR, 'dir') +
                  HHBLITS_DATABASE_NAME + ' -n ' + hhparam[0] +
                  ' -e ' + hhparam[1] + ' -diff ' + hhparam[2] +
                  ' -cov ' + hhparam[3] + ' -id ' + hhparam[4] +
                  ' -oa3m ' + msaa3mpath)
    except:
        try:
            os.system(hhsuite_exec + ' -i '+ spath +
                      ' -d ' + pio.check_path(HHBLITS_DATABASE_DIR, 'dir') +
                      HHBLITS_DATABASE_NAME + ' -n ' + hhparam[0] +
                      ' -e ' + hhparam[1] + ' -diff ' + hhparam[2] +
                      ' -cov ' + hhparam[3] + ' -id ' + hhparam[4] +
                      ' -oa3m ' + msaa3mpath)
        except:
            logging.critical('        An error occurred while executing HHBLITS.')


    # Convert A3M MSA file to Jones format (DMP standard input format)
    parsedmsa, msajonespath = msafilesgen(msaa3mpath)
    logging.info('    Done\n')

    return parsedmsa, msajonespath

def msafilesgen(inpath_a3m):
    """Convert a3m format alignment to "jones" format (.aln) and print msa coverage graph.

    :param inpath_a3m: Multiple Sequence Alignment (.a3m) path
    :type inpath_a3m: str
    :return: ConKit parsed MSA
    :rtype: :obj:`~conkit.core.sequencefile.SequenceFile`

    """
    parsedmsa=ckio.read(inpath_a3m,'a3m')

    # Convert to 'jones' format
    outpath_aln = os.path.join(os.path.splitext(inpath_a3m)[0],".aln")
    ckio.write(outpath_aln,'jones',parsedmsa)

    # Plot Coverage
    msacoveragepath = os.path.join(os.path.splitext(inpath_a3m)[0],".coverage.png")
    fig = ckplot.SequenceCoverageFigure(parsedmsa)
    fig.savefig(msacoveragepath, overwrite=True)

    neff=parsedmsa.meff

    return neff
