"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""
from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import conf as pcnf

import os
import argparse
import glob

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
        if os.path.isdir(os.path.join(path, '')) is True:
            path = os.path.abspath(os.path.join(path, ''))
            pathok = True
        else:
            if os.path.isfile(path) is True:
                pathok = True
    elif typeofpath == 'file':
        path = os.path.abspath(path)
        if os.path.isfile(path) is True:
            pathok = True
    elif typeofpath is None:
        path = os.path.abspath(path)
        pathok = True
    else:
        raise ValueError("Input string 'typeofpath' should be either 'dir' or 'file'.")
    if pathok is True:
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def check_wildcard(path):
    """Return list of several files full path.

    :param path: Path containing wildcard '*'.
    :type path: str
    :return: List of absolute paths.
    :rtype: list [str]
    """

    abspath = check_path(path)

    return glob.glob(abspath)


def mdir(dirpath):
    """Create directory recursively if not existent.

    :param dirpath: Path to directory to be created.
    :type dirpath: str
    """
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    return

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
