# -*- coding: utf-8 -*-

import os
import sys
#from StringIO import StringIO  # Python2
from io import StringIO  # Python3

from .inputvalues import OUTPUT_DIR
from .inputvalues import LOG_STDOUT

from .init2 import pdbid
from .init2 import gettime

#####################################################################
##### OUTPUT FUNCTIONS ##############################################

def output_dir():
    """
    Returns path for output directory

    Returns
    -------
    outpath : str
        Output directory

    """
    #pdbid=os.path.splitext(os.path.basename(PDB_PATH))[0]
    outpath = os.path.join(OUTPUT_DIR, pdbid(),"")   # OUTPUT DIRECTORY WHERE OUTPUT FILES WILL GO

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

def stdoutpath():
    """
    Path to output logfile.

    Returns
    -------
    stdout : str
        Path to ouput file.

    """
    stdoutfile=pdbid()+".stdout.log"
    stdout = os.path.join(output_dir(), stdoutfile)

    return stdout

def errorpath():
    """
    Path to error logfile.

    Returns
    -------
    stdout : str
        Path to ouput file.

    """
    stdoutfile=pdbid()+".error.log"
    stdout = os.path.join(output_dir(), stdoutfile)

    return stdout

def printout(string,extraline=False,errorlog=False):
    """
    Writes string into the log file(s)

    Parameters
    ----------
    string : str (if not a string, try to keep stdout)
        output string
    extraline : bool, optional
        Add extra blank line after string. The default is False.
    errorlog : bool, optional
        Print this information to the error log too. The default is False.

    Returns
    -------
    None.

    """
    leap='\n' if not extraline else '\n\n'

    try:
        string = string.decode()
    except (UnicodeDecodeError, AttributeError):
        pass

    if isinstance(string,list):
        string="\n".join(str(row) for row in string)

    elif not isinstance(string,(str,list)):
        old_stdout = sys.stdout
        output = StringIO()
        sys.stdout=output
        print(string)
        sys.stdout = old_stdout
        string = output.getvalue()

    if LOG_STDOUT:
        with open(stdoutpath(), 'a') as out:
            out.write(string+leap)
    else:
        print(string+leap)

    if errorlog:
        try:
            with open(errorpath(), 'a') as out:
                out.write(string+leap)
        except:
            with open(errorpath(), 'w') as out:
                out.write(string+leap)

def sysfileout():
    """
    If systems allow, print stdout of external processes to stdout file

    Returns
    -------
    fileout : str
        Append stdout to file command.

    """

    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        fileout=' >> '+stdoutpath()
    else:
        fileout=''

    return fileout

def interrupt(errstring):
    """
    Interrupt execution.

    Parameters
    ----------
    errstring : str
        Error message to be returned.

    Returns
    -------
    None.

    """
    inttime = gettime()
    printout("***********************************",errorlog=True)
    printout("***   I N T E R R U P T I O N   ***",errorlog=True)
    printout("***********************************",errorlog=True, extraline=True)
    printout(errstring,errorlog=True)
    printout('Interruption occurred at ' + inttime[1].strftime("%-d %B %Y, %X") +' UTC',errorlog=True)
    sys.exit(errstring)

