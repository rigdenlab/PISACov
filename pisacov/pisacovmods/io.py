"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

import os
import sys
import argparse


def check_path(path, typeofpath=None):
    """Returns full path if correct.
    :param path: Input (local) path.
    :type path: str
    :param typeofpath: The type of path, 'dir' or 'file', defaults to None.
    :type typeofpath: str, optional
    :raises ValueError: When given typeofpath is neither 'dir' nor 'file'.
    :raises argparse: If wrong path given.
    :return: Complete checked path.
    :rtype: str
    """
    pathok = False
    if typeofpath == 'dir':
        path = os.path.abspath(path)
        if os.path.isdir(os.path.join(path, '')) is True:
            path = os.path.abspath(os.path.join(path, ''))
            pathok = True
    elif typeofpath == 'file':
        path = os.path.abspath(path)
        if os.path.isfile(path) is True:
            pathok = True
    elif typeofpath is None:
        if os.path.isdir(os.path.abspath(os.path.join(path, ''))) is True:
            path = os.path.abspath(os.path.join(path, ''))
            pathok = True
        else:
            path = os.path.abspath(path)
            if os.path.isfile(path) is True:
                pathok = True
    else:
        raise ValueError("Input string 'typeofpath' should be either 'dir' or 'file'.")
    if pathok:
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def output_dir(rootlevel, extralevel=None):
    """
    Returns path for output directory

    :param extralevel: additional directory level, defaults to None
    :type extralevel: str, optional
    :return: outpath
    :rtype: Output directory

    """
    outpath = os.path.join()
    outpath = os.path.join(OUTPUT_DIR, pdbid(), "")

    return outpath


def output_tmpdir(extrapath=None):
    """
    Returns path for output directory

    Returns
    -------
    outpath : str
        Output directory

    """
    if extrapath == None:
        outpath = os.path.join(OUTPUT_DIR, pdbid(),"tmp","")   # OUTPUT DIRECTORY WHERE OUTPUT FILES WILL GO
    else:
        outpath = os.path.join(OUTPUT_DIR, pdbid(),"tmp",extrapath,"")   # OUTPUT DIRECTORY WHERE OUTPUT FILES WILL GO

    return outpath

def mkdirout():
    """
    Create output directory

    Returns
    -------
    None.

    """
    #pdbid=os.path.splitext(os.path.basename(PDB_PATH))[0]
    #outdir = os.path.join(OUTPUT_DIR, pdbid(),"")   # OUTPUT DIRECTORY WHERE OUTPUT FILES WILL GO

    if os.path.exists(output_dir()):
        sys.exit("ERROR. Unable to create output directory. %s already exists. Please, make sure you choose an output path not containing former results." % output_dir() ) # LOGGING?
    else:
        try:
            os.mkdir(output_dir())
        except OSError:
            sys.exit("ERROR. Unable to create output directory %s." % output_dir() )
    os.mkdir(output_tmpdir())
    os.mkdir(output_tmpdir("pisacov"))
    os.mkdir(output_tmpdir("pisa"))
    os.mkdir(output_tmpdir("deepmetapsicov"))
#    if LOG_STDOUT:
#        with open(stdoutpath(), 'w') as out:
#            pass
