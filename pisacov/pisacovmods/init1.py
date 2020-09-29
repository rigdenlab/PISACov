# -*- coding: utf-8 -*-

from .inputvalues import NEIGHBOURS_MINDISTANCE
from .inputvalues import HHBLITS_PARAMETERS
from .inputvalues import HHBLITS_VIA_DMP
from .inputvalues import DMPSCORE_THRESHOLD
from .inputvalues import PSICOVSCORE_THRESHOLD
from .inputvalues import REMOVE_INTRA_CONTACTS

from .output import printout

#####################################################################
##### ININITIALISATION FUNCTIONS ####################################

def minneigh(logit=None):
    """
    Check that this variable has the appropriate format and reformat otherwise
    ----------
    logit : str, optional
        If not none, it activates the logging. The default is None.

    Returns
    -------
    newvalue : int
        Corrected value for NEIGHBOURS_MINDISTANCE.

    """
    if REMOVE_INTRA_CONTACTS:
        newvalue=2
    else:
        newvalue=NEIGHBOURS_MINDISTANCE
        if isinstance(NEIGHBOURS_MINDISTANCE,int):
            if (NEIGHBOURS_MINDISTANCE<2):
                if logit is not None:
                    printout('WARNING: NEIGHBOURS_MINDISTANCE can only assume values of 2 or higher. This value has been now changed to 2.',errorlog=True)
                newvalue=2
        elif isinstance(NEIGHBOURS_MINDISTANCE,float):
            if logit is not None:
                printout('WARNING: NEIGHBOURS_MINDISTANCE can only assume integer values. To avoid crashing, its value has been updated from', NEIGHBOURS_MINDISTANCE, ' to ', round(NEIGHBOURS_MINDISTANCE),errorlog=True)
            newvalue=round(NEIGHBOURS_MINDISTANCE)
            if (newvalue<2):
                if logit is not None:
                    printout('WARNING: NEIGHBOURS_MINDISTANCE can only assume values of 2 or higher. This value has been now changed to 2.',errorlog=True)
                newvalue=2
        else:
            if logit is not None:
                printout("ERROR: NEIGHBOURS_MINDISTANCE is not an integer. Assigning default value = 2.", errorlog=True)

    return newvalue

def scorethreshold(source,logit=None):
    """
    Check that this variable has the appropriate format and reformat otherwise
    ----------
    source : str
        Source of the contact list ("deepmetapsicov" or "psicov")
    logit : str, optional
        If not none, it activates the logging. The default is None.

    Returns
    -------
    newvalue : float
        Corrected value for DMPSCORE_THRESHOLD.

    """
    if source == "deepmetapsicov":
        if not DMPSCORE_THRESHOLD:
            if isinstance(DMPSCORE_THRESHOLD,int) or isinstance(DMPSCORE_THRESHOLD,float):
                newvalue=0.0
            else:
                if logit is not None:
                    printout('WARNING: A DMPSCORE_THRESHOLD value was not given. Setting it to default value 0.2.', errorlog=True)
                newvalue=0.2
        else:
            if isinstance(DMPSCORE_THRESHOLD,int) or isinstance(DMPSCORE_THRESHOLD,float):
                if DMPSCORE_THRESHOLD<0 or DMPSCORE_THRESHOLD>=1:
                    if logit is not None:
                        printout('DMPSCORE_THRESHOLD = '+DMPSCORE_THRESHOLD+' but it should be within the interval [0,1). Setting it to default value 0.2.', errorlog=True)
                    newvalue=0.2
                else:
                    newvalue=DMPSCORE_THRESHOLD
            else:
                if logit is not None:
                    printout('WARNING: DMPSCORE_THRESHOLD must be a number within the interval [0,1). Setting it to default value 0.2.', errorlog=True)
                newvalue=0.2
    elif source == "psicov":
        if not PSICOVSCORE_THRESHOLD:
            if isinstance(PSICOVSCORE_THRESHOLD,int) or isinstance(PSICOVSCORE_THRESHOLD,float):
                newvalue=0.0
            else:
                if logit is not None:
                    printout('WARNING: A PSICOVSCORE_THRESHOLD value was not given. Setting it to default value 0.2.', errorlog=True)
                newvalue=0.2
        else:
            if isinstance(PSICOVSCORE_THRESHOLD,int) or isinstance(PSICOVSCORE_THRESHOLD,float):
                if PSICOVSCORE_THRESHOLD<0:
                    if logit is not None:
                        printout('PSICOVSCORE_THRESHOLD = '+PSICOVSCORE_THRESHOLD+' but it should be within the interval [0,inf). Setting it to default value 0.2.', errorlog=True)
                    newvalue=0.2
                else:
                    newvalue=PSICOVSCORE_THRESHOLD
            else:
                if logit is not None:
                    printout('WARNING: PSICOVSCORE_THRESHOLD must be a number within the interval [0,inf). Setting it to default value 0.2.', errorlog=True)
                newvalue=0.2
    return newvalue

def hhparam(logit=None):
    """
    Default HHBLITS PARAMETERS if not defined

    Parameters
    ----------
    logit : str, optional
        If not none, it activates the logging. The default is None.

    Returns
    -------
    newvalue : list [int,float]
        Input parameters for HHSUITE run.

    """
    if HHBLITS_VIA_DMP:
        if logit is not None:
            printout('Using default HHBLITS_PARAMETERS given by the DeepMetaPSICOV script: [ 3 , 0.001, inf, 50, 99 ] ')
        newvalue=[ 3 , 0.001, "inf", 50, 99 ]
    else:
        if not HHBLITS_PARAMETERS:
            if logit is not None:
                printout('WARNING: HHBLITS_PARAMETERS values not given. Setting parameters to default [ 2 , 0.001, 1000, 0, 90 ] ', errorlog=True)
            newvalue = [ 2 , 0.001, 1000, 0, 90 ]
        else:
            newvalue = HHBLITS_PARAMETERS

    return newvalue