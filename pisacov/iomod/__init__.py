"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import logging
from crops.iomod import parsers as cps
from conkit import io as ckio
import xml.etree.ElementTree as ET
import numpy as np

def _ftypelist():
    lout = ['psicov', 'fasta', 'pdb', 'a3m', 'jones', 'xml']
    return lout

def read(infile, ftype, ck=False):
    """

    :param infile: Path to input file.
    :type infile: str
    :param ftype: One of 'psicov', 'ccmpred', 'fasta', 'pdb', 'a3m', 'jones', 'xml', 'array'.
    :type ftype: str
    :param ck: Open alternative conkit version instead of default, defaults to False.
    :type ck: bool, optional
    :return: Parsed file (and, for 'pdb', list of filenames).
    :rtype: One or two of list[str], :class:`~crops.elements.sequences.sequence`, :class:`~conkit.core.sequence.Sequence`,

    """
    if (ftype.lower() not in _ftypelist() or
            isinstance(ftype, str) is not True):
        logging.critical('Specified type not valid.')
        raise ValueError

    if ck is True and ftype.lower() != 'xml' and ftype.lower() != 'array':
        output = ckio.read(infile, ftype.lower())
    else:
        if ftype.lower() == 'psicov':
            pass
        if ftype.lower() == 'ccmpred':
            pass
        elif ftype.lower() == 'fasta':
            output = cps.parseseqfile(infile)
        elif ftype.lower() == 'pdb':
            output1, output2 = cps.parsestrfile(infile)
            return output1, output2
        elif ftype.lower() == 'a3m' or 'jones':
            output = ckio.read(infile, ftype.lower())
        elif ftype.lower() == 'xml':
            output = ET.parse(infile)
        elif ftype.lower() == 'array':
            output = np.loadtxt(infile)
    return output


def write(infile, ftype, indata, ck=False):
    """

    :param infile: Path to input file.
    :type infile: str
    :param ftype: One of 'psicov', 'ccmpred', 'fasta', 'pdb', 'a3m', 'jones', 'xml'.
    :type ftype: str
    :param ck: Open alternative conkit version instead of default, defaults to False.
    :type ck: bool, optional
    :return: Parsed file (and, for 'pdb', list of filenames).
    :rtype: One or two of list[str], :class:`~crops.elements.sequences.sequence`, :class:`~conkit.core.sequence.Sequence`,

    """
    if (ftype.lower() not in _ftypelist() or
            isinstance(ftype, str) is not True):
        logging.critical('Specified type not valid.')
        raise ValueError

    if ck is True and ftype.lower() != 'xml':
        output = ckio.write(infile, ftype.lower(), hyerarchy=indata)
    else:
        if ftype.lower() == 'psicov':
            pass
        if ftype.lower() == 'ccmpred':
            pass
        elif ftype.lower() == 'fasta':
            output = cps.parseseqfile(infile)
        elif ftype.lower() == 'pdb':
            output1, output2 = cps.parsestrfile(infile)
            return output1, output2
        elif ftype.lower() == 'a3m' or 'jones':
            output = ckio.read(infile, ftype.lower())
        elif ftype.lower() == 'xml':
            output = ET.parse(infile)
    return output
