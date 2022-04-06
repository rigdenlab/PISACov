"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.io import conf as pconf
from pisacov.io.paths import check_path
import logging

_surl = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz'
_uniurl = 'https://www.uniprot.org/uniprot/'

def _default_values(key):
    defaultvals = {}
    defaultvals['HHBLITS_PARAMETERS'] = [3, 0.001, 'inf', 50, 99]
    defaultvals['UNICLUST_FASTA_PATH'] = None
    defaultvals['NEIGHBOURS_MINDISTANCE'] = 2
    defaultvals['REMOVE_INTRA_CONTACTS'] = True

    return defaultvals[key]

def _default_keys():

    defaultkeys = ['SIFTS_PATH', 'PISA_PATH', 'DMP_PATH', 'HHBLITS_PATH',
                 'HHBLITS_DATABASE_NAME', 'HHBLITS_DATABASE_DIR', 'HHBLITS_PARAMETERS',
                 'UNICLUST_FASTA_PATH', 'NEIGHBOURS_MINDISTANCE', 'REMOVE_INTRA_CONTACTS']

    return defaultkeys

def _check_input(val, key):

    compmsg = {'SIFTS_PATH': "SIFTS_PATH file not found.",
               'PISA_PATH': "'PISA_PATH file not found.",
               'HHBLITS_PATH': "HHBLITS_PATH file not found.",
               'HHBLITS_DATABASE_NAME': "HHBLITS_DATABASE_DIR should be a string.",
               'HHBLITS_DATABASE_DIR': "HHBLITS_DATABASE_DIR directory not found.",
               'DMP_PATH': "DMP_PATH file not found.",
               'HHBLITS_PARAMETERS': "Format should be [int, float, int or 'inf', int, int].",
               'UNICLUST_FASTA_PATH': "UNICLUST_FASTA_PATH file not found.",
               'NEIGHBOURS_MINDISTANCE': "NEIGHBOURS_MINDISTANCE should be an integer.",
               'REMOVE_INTRA_CONTACTS': "REMOVE_INTRA_CONTACTS should be a boolean."}
    errormsg = compmsg[key] + ' Please, update the configuration file using pisaconf.'

    if (key == 'SIFTS_PATH' or key == 'PISA_PATH' or
            key == 'HHBLITS_PATH' or key == 'DMP_PATH'):
        try:
            val = check_path(val, 'file')
        except NameError:
            logging.critical(errormsg)
    if (key == 'UNICLUST_FASTA_PATH'):
        if val == "" or val is None:
            val = None
        else:
            try:
                val = check_path(val, 'file')
            except NameError:
                logging.critical(errormsg)
    elif (key == 'HHBLITS_DATABASE_DIR'):
        try:
            val = check_path(val, 'dir')
        except NameError:
            logging.critical(errormsg)
    elif (key == 'HHBLITS_DATABASE_NAME'):
        if isinstance(val, str) is False:
            try:
                val = str(val)
            except Exception:
                logging.critical(errormsg)
        else:
            pass
    elif (key == 'HHBLITS_PARAMETERS'):
        if val is None or val == "":
            val = _default_values(key)
        else:
            if isinstance(val, list) is False:
                logging.critical('HHBLITS_PARAMETERS is not a python list.')
            else:
                for n in range(len(val)):
                    if n == 1:
                        if isinstance(val[n], float) is False:
                            try:
                                val[n] = float(val[n])
                            except TypeError:
                                logging.critical(errormsg)
                    else:
                        if n == 2 and val[n] == 'inf':
                            pass
                        else:
                            if isinstance(val[n], int) is False:
                                try:
                                    val[n] = int(val[n])
                                except TypeError:
                                    logging.critical(errormsg)
    elif (key == 'NEIGHBOURS_MINDISTANCE'):
        if val is None or val == "":
            val = _default_values(key)
        else:
            if isinstance(val[n], int) is False:
                try:
                    val = int(val)
                except TypeError:
                    logging.critical(errormsg)
    elif (key == 'REMOVE_INTRA_CONTACTS'):
        if val is None:
            val = _default_values(key)
        elif isinstance(val, str) is True:
            if val == "":
                val = _default_values(key)
            elif val.lower() == "true":
                val = True
            elif val.lower() == "false":
                val = False
            else:
                logging.critical(errormsg)
                raise TypeError
        elif isinstance(val, bool) is True:
            pass
        else:
            logging.critical(errormsg)
            raise TypeError

    return val

def _parse_conf():

    trymsg = ' Please, update the configuration file using pisaconf.'
    compmsg = 'One or more of the requested parameters not found. '
    try:
        configfile = {'SIFTS_PATH': pconf.SIFTS_PATH,
                      'PISA_PATH': pconf.PISA_PATH,
                      'HHBLITS_PATH': pconf.HHBLITS_PATH,
                      'HHBLITS_DATABASE_NAME': pconf.HHBLITS_DATABASE_NAME,
                      'HHBLITS_DATABASE_DIR': pconf.HHBLITS_DATABASE_DIR,
                      'DMP_PATH': pconf.DMP_PATH}
    except NameError:
        logging.critical(compmsg + trymsg)

    try:
        configfile['UNICLUST_FASTA_PATH'] = pconf.UNICLUST_FASTA_PATH
    except Exception:
        configfile['UNICLUST_FASTA_PATH'] = None

    try:
        configfile['HHBLITS_PARAMETERS'] = pconf.HHBLITS_PARAMETERS
    except Exception:
        configfile['HHBLITS_PARAMETERS'] = None

    try:
        configfile['NEIGHBOURS_MINDISTANCE'] = pconf.NEIGHBOURS_MINDISTANCE
        configfile['REMOVE_INTRA_CONTACTS'] = False
    except Exception:
        configfile['NEIGHBOURS_MINDISTANCE'] = None
        configfile['REMOVE_INTRA_CONTACTS'] = True

    config = {}

    for keystr, value in configfile.items():
        config[keystr] = _check_input(value, keystr)

    return config


def _initialise_inputs():

    outvalues = _parse_conf()

    outvalues['INSEQ'] = None
    outvalues['INSTR'] = None
    outvalues['ALTDB'] = None
    outvalues['OUTROOT'] = None
    outvalues['OUTCSVPATH'] = None
    outvalues['UPTHRESHOLD'] = None

    return outvalues


def _check_hhparams(paramlist):
    """Return a list with the validated HHBLITS input arguments.

    :param paramlist: User-provided list of HHBLITS arguments.
    :type paramlist: list of str
    :raises ValueError: Arguments are not valid.
    :return: Complete and checked list of HHBLITS parameters
    :rtype: list of (int, float, str)
    """
    if (paramlist == 'dmp' or paramlist == [3, 0.001, 'inf', 50, 99] or
        paramlist == ['3', '0.001', 'inf', '50', '99']):
        outparams = ['3', '0.001', 'inf', '50', '99']
    elif (paramlist == 'hhblits' or paramlist == [2, 0.001, 1000, 0, 90] or
        paramlist == ['2', '0.001', '1000', '0', '90']):
        outparams = ['2', '0.001', '1000', '0', '90']
    else:
        try:
            int(float(paramlist[0]))
            float(paramlist[1])
            if paramlist[2] != 'inf':
                int(float(paramlist[2]))
            float(paramlist[4])
            float(paramlist[5])
        except:
            raise ValueError('One or more of HHblits arguments given are not valid')

        outparams=[str(int(float(paramlist[0]))), str(float(paramlist[1])),
            paramlist[2], str(float(paramlist[4])), str(float(paramlist[5]))]
        if paramlist[2] != 'inf':
            outparams = str(int(float(paramlist[2])))

    return outparams

def _check_uniprot(inuniprot):
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
        dbpath = check_path(pconf.UNICLUST_FASTA_PATH, 'file')
    except:
        dbpath = 'server-only'

    return float(inuniprot), dbpath

def _sourcenames(short=False):

    """Return a list with the source names.

    :param short: True for shorter names, defaults to False
    :type short: bool, optional
    :return: Source names.
    :rtype: dict [list [str]]

    """
    if short is False:
        sources = ["psicov", "ccmpred", "deepmetapsicov"]
    else:
        sources =["psicov", "ccmpred", "dmp"]

    return sources

def _sources():
    """Return the subdir name and extension of each of the contact prediction types.

    :return: Contact prediction types and location.
    :rtype: dict [list [str]]

    """
    sources = _sourcenames()
    confiledir = ["deepmetapsicov", "deepmetapsicov", "deepmetapsicov"]
    confilesuffix = ["psicov", "ccmpred", "deepmetapsicov.con"]
    conkittype = ["psicov", "ccmpred", "psicov"]

    outsinfo = {}
    for n in range(len(sources)):
        outsinfo[sources[n]] = [confiledir[n], confilesuffix[n], conkittype[n]]

    return outsinfo
