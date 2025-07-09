"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import os
import logging

from pisacov.iomod.conf import PISA_PATH
from pisacov.core import interfaces as pci


def runpisa(instr, outdir, sessionid=None):
    """Run PISA and generate interface PDB files for a given PDB file.

    :param instr: Input structure filepath.
    :type instr: str
    :param outdir: Directory where results of PISA will be printed out.
    :type outdir: str
    :param sessionid: Local name for PISA session (whole name will read 'pisacov_sessionid'), defaults to None.
    :type sessionid: str, optional

    :return: Interfaces.
    :rtype: :class:`~pisacov.core.contacts.contact_atlas`

    """
    pisa_exec = '"' + PISA_PATH + '"'
    sname = 'pisacov'
    if sessionid is not None:
        sname += '_' + sessionid
    # CLEAR SESSION
    try:
        oc = os.system(pisa_exec + ' ' + sname + ' -erase')
        if oc != 0:
            raise Exception
    except Exception:
        logging.critical('        An error occurred while executing PISA -erase.')
        raise OSError

    # ANALYSE PDB STRUCTURE
    try:
        oc = os.system(pisa_exec + ' ' + sname + ' -analyse ' + instr)
        if oc != 0:
            raise Exception
    except Exception:
        logging.critical('        An error occurred while executing PISA -analyse.')
        raise OSError

    # CREATE XML FILES
    xmlfile = (os.path.splitext(os.path.basename(instr))[0] +
               os.extsep + "interface" + os.extsep + "xml")
    ixmlpath = os.path.join(outdir, xmlfile)
    try:
        oc = os.system(pisa_exec + ' ' + sname + ' -xml interface > ' + ixmlpath)
        if oc != 0:
            raise Exception
    except Exception:
        logging.critical('        An error occurred while executing PISA -xml interface.')
        raise OSError

    xmlfile = (os.path.splitext(os.path.basename(instr))[0] +
               os.extsep + "assembly" + os.extsep + "xml")
    axmlpath = os.path.join(outdir, xmlfile)
    try:
        oc = os.system(pisa_exec + ' ' + sname + ' -xml assembly > ' + axmlpath)
        if oc != 0:
            raise Exception
    except Exception:
        logging.critical('        An error occurred while executing PISA -xml assembly.')
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
            oc = os.system(pisa_exec + ' ' + sname +' -pdb interface ' +
                           str(nif) + ' > ' + os.path.join(outdir, pdbfile))
            if oc != 0:
                raise Exception
        except Exception:
            logging.critical("ERROR: An error occurred during the execution of PISA (interface pdb files production).")
            raise OSError

    return ilist
