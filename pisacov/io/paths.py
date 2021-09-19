"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov.io import conf as pcnf

import os
import argparse

import logging

def check_path(path, typeofpath=None):
    """Return full path. If typeofpath is given, path is checked.

    :param path: Input (local) path.
    :type path: str
    :param typeofpath: The type of path, 'dir', 'file' or 'either', defaults to None.
    :type typeofpath: str, optional
    :raises ValueError: When given typeofpath is none of 'dir', 'file' or 'either'.
    :raises argparse: If wrong path or pathtype given.
    :return: Absolute path.
    :rtype: str
    """
    pathok = False
    if typeofpath == 'dir' or typeofpath == 'either':
        path = os.path.abspath(path)
        if os.path.isdir(os.path.join(path, '')):
            path = os.path.abspath(os.path.join(path, ''))
            pathok = True
    elif typeofpath == 'file' or typeofpath == 'either':
        path = os.path.abspath(path)
        if os.path.isfile(path):
            pathok = True
    elif typeofpath is None:
        path = os.path.abspath(path)
        pathok = True
    else:
        raise ValueError("Input string 'typeofpath' should be either 'dir' or 'file'.")
    if pathok:
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def check_hhparams(paramlist):
    """Return a list with the validated HHBLITS input arguments.

    :param paramlist: User-provided list of HHBLITS arguments.
    :type paramlist: list of str
    :raises ValueError: Arguments are not valid.
    :return: Complete and checked list of HHBLITS parameters
    :rtype: list of (int, float, str)
    """
    try:
        int(float(paramlist[0]))
        float(paramlist[1])
        if paramlist[2] != 'inf':
            int(float(paramlist[2]))
        float(paramlist[4])
        float(paramlist[5])
    except:
        raise ValueError('One or more of HHblits arguments given are not valid')

    outlist=[str(int(float(paramlist[0]))), str(float(paramlist[1])),
             paramlist[2], str(float(paramlist[4])), str(float(paramlist[5]))]
    if paramlist[2] != 'inf':
        outlist[2] = str(int(float(paramlist[2])))

    return [outlist]

def check_uniprot(inuniprot):
    """Return Uniprot segment threshold value and Uniprot database path.

    :param inuniprot: Initial argument for Uniprot threshold
    :type inuniprot: str
    :raises ValueError: Argument is not valid.
    :return: Threshold value and database path.
    :rtype: float, str

    """
    try:
        float(inuniprot)
    except:
        raise ValueError('Uniprot threshold given not valid.')

    try:
        dbpath = check_path(pcnf.UNICLUST_FASTA_PATH, 'file')
    except:
        dbpath = 'https://www.uniprot.org/uniprot/'

    return float(inuniprot), dbpath

def output_dir(root,baseid):
    """Return path to output directory.

    :param baseid: PDB ID or other id to identify the molecule.
    :type baseid: str
    :return: Path to output directory.
    :rtype: str

    """
    if not isinstance(baseid, str):
        try:
            str(baseid)
        except:
            raise ValueError('Not a valid argument. String expected.')
    outpath = os.path.join(root, str(baseid).lower(),"")

    return outpath

def mkdirout(root,baseid=None):
    """Check or create output directories.

    :param root: Path to output directory.
    :type root: str
    :param baseid: PDB ID of molecule been analised, defaults to None.
    :type baseid: str
    :raises OSError: If the system does not allow to create the directories.
    :return: Dictionary containing absolute paths to output directories.
    :rtype: dict of str

    """
    outpaths={}

    try:
        outpaths['root'] = check_path(root, 'dir')
    except:
        try:
            os.mkdir(os.path.join(root,""))
        except:
            errormsg = "Unable to create output directory " + str(os.path.join(root,""))
            raise OSError(errormsg)
        outpaths['root'] = check_path(root, 'dir')

    if baseid is not None:
        try:
            outpaths['pdbid'] = check_path(os.path.join(outpaths['root'], baseid, ""), 'dir')
        except:
            try:
                os.mkdir(os.path.join(outpaths['root'], baseid, ""))
            except:
                errormsg = "Unable to create output directory " + str(os.path.join(outpaths['root'], baseid, ""))
                raise OSError(errormsg)
            outpaths['pdbid'] = check_path(os.path.join(outpaths['root'], baseid, ""), 'dir')

        try:
            outpaths['pisa'] = check_path(os.path.join(outpaths['pdbid'], 'pisa_interfaces', ""), 'dir')
        except:
            try:
                os.mkdir(os.path.join(outpaths['pdbid'], 'pisa_interfaces', ""))
            except:
                errormsg = "Unable to create output directory " + str(os.path.join(outpaths['pdbid'], 'pisa_interfaces', ""))
                raise OSError(errormsg)
            outpaths['pisa'] = check_path(os.path.join(outpaths['pdbid'], 'pisa_interfaces', ""), 'dir')

        try:
            outpaths['pisacov'] = check_path(os.path.join(outpaths['pdbid'], 'pisacov', ""), 'dir')
        except:
            try:
                os.mkdir(os.path.join(outpaths['pdbid'], 'pisacov', ""))
            except:
                errormsg = "Unable to create output directory " + str(os.path.join(outpaths['pdbid'], 'pisacov', ""))
                raise OSError(errormsg)
            outpaths['pisacov'] = check_path(os.path.join(outpaths['pdbid'], 'pisacov', ""), 'dir')

        try:
            outpaths['dmp'] = check_path(os.path.join(outpaths['pdbid'], 'dmp_predictions', ""), 'dir')
        except:
            try:
                os.mkdir(os.path.join(outpaths['pdbid'], 'dmp_predictions', ""))
            except:
                errormsg = "Unable to create output directory " + str(os.path.join(outpaths['pdbid'], 'dmp_predictions', ""))
                raise OSError(errormsg)
            outpaths['dmp'] = check_path(os.path.join(outpaths['pdbid'], 'dmp_predictions', ""), 'dir')

    return outpaths
