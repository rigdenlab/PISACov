# -*- coding: utf-8 -*-

import os
import sys
import datetime
import time

from .inputvalues import PDB_PATH

#####################################################################
##### ININITIALISATION FUNCTIONS ####################################

def pdbid():
    """
    PDB ID from PDB file.

    Returns
    -------
    locpdbid : str
        PDB ID (lower case).

    """
    locpdbid=os.path.splitext(os.path.basename(PDB_PATH))[0]

    return locpdbid

def gettime():
    """
    Returns date, time and wallclock time

    Returns
    -------
    list
        [Wallclock time, date and time].

    """
    if sys.platform == 'win32':
        # On Windows, the best timer is time.clock
        wallclock=time.clock()
    else:
        # On most other platforms the best timer is time.time
        wallclock= time.time()

    nowdt = datetime.datetime.now(datetime.timezone.utc)
    return [wallclock, nowdt]

def readabletime(inputtime):
    if inputtime >= 60:
        inputtime /= 60
        mins = int(inputtime)
        secs=(inputtime - mins)*60
        outstring=str(round(secs,2)) + ' seconds'
        inputtime=mins
        if inputtime >= 60:
            inputtime /= 60
            hours = int(inputtime)
            mins = (inputtime - hours)*60
            outstring = str(round(mins,0))+ ' minutes '+ outstring
            inputtime=hours
            if inputtime >= 24:
                inputtime /= 24
                days = int(inputtime)
                hours = (inputtime - days)*24
                outstring = str(days)+ ' days '+str(round(hours,0))+ ' hours '+ outstring
            else:
                outstring = str(round(hours,0))+ ' hours '+ outstring
        else:
            outstring = str(round(mins,0))+ ' minutes '+ outstring
    else:
        secs = inputtime
        outstring=str(secs) + ' seconds'

    return outstring

def ressymbol(name):
    """
    Dictionary containing conversion from residue 3-letter symbol to 1-letter symbol

    Parameters
    ----------
    name : str
        Residue symbol (3-letter convention)

    Returns
    -------
    oneletter : str
        Residue symbol (1-letter convention)

    """

    reslist =	{
      "ALA": "A",
      "ARG": "R",
      "ASN": "N",
      "ASP": "D",
      "CYS": "C",
      "GLN": "Q",
      "GLU": "E",
      "GLY": "G",
      "HIS": "H",
      "ILE": "I",
      "LEU": "L",
      "LYS": "K",
      "MSE": "M",
      "MET": "M",
      "PHE": "F",
      "PRO": "P",
      "SER": "S",
      "THR": "T",
      "TRP": "W",
      "TYR": "Y",
      "VAL": "V",
      "SEC": "U",
      "PYL": "O",
      "XAA": "X",
      "ASX": "B",
      "GLX": "Z",
      "XLE": "J"
    }

    oneletter=reslist[name]
    return oneletter