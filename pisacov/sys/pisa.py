"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import os
import logging
import xml.etree.ElementTree as ET

from pisacov.io.conf import PISA_PATH
from pisacov import io as pio
from pisacov.core import interfaces as pci

def runpisa(instr, outdir):
    """
    Run PISA and generate interface PDB files for given PDB

    :param instr: Input structure filepath.
    :type instr: str
    :param outdir: Directory where results of PISA will be printed out.
    :type outdir: str

    :return: Interfaces.
    :rtype: :class:`~pisacov.core.contacts.contact_atlas`

    """
    pisa_exec = '"' + PISA_PATH + '"'
    # ANALYSE PDB STRUCTURE
    try:
        os.system(pisa_exec + ' pisacov -analyse ' + instr)
    except Exception:
        logging.critical('        An error occurred while executing PISA.')
        raise OSError

    # CREATE XML FILES
    xmlfile = (os.path.splitext(os.path.basename(instr))[0] +
               os.extsep + "interface" + os.extsep + "xml")
    ixmlpath = os.path.join(outdir, xmlfile)
    try:
        os.system(pisa_exec + ' pisacov -xml interface > ' + ixmlpath)
    except Exception:
        logging.critical('        An error occurred while executing PISA.')
        raise OSError

    xmlfile = (os.path.splitext(os.path.basename(instr))[0] +
               os.extsep + "assembly" + os.extsep + "xml")
    axmlpath = os.path.join(outdir, xmlfile)
    try:
        os.system(pisa_exec + ' pisacov -xml assembly > ' + axmlpath)
    except Exception:
        logging.critical('        An error occurred while executing PISA.')
        raise OSError

    # Obtain interface list
    ilist = pci.parse_interface_xml(ixmlpath, axmlpath)

    # Generate PDB files for each interface
    logging.info("        Generating interface PDB files...")
    for nif in range(1, len(ilist)+1):
        # Write pdb files
        pdbfile = (os.path.splitext(os.path.basename(instr))[0] + os.extsep +
                   "interface" + os.extsep + str(nif) + os.extsep + "pdb")
        try:
            os.system(pisa_exec + ' pisacov -pdb interface ' + str(nif) +
                      ' > ' + os.path.join(outdir, pdbfile))
        except Exception:
            logging.critical("ERROR: An error occurred during the execution of PISA (interface pdb files production).")
            raise OSError

    logging.info('    Done\n')

    return ilist
