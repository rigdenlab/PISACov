"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

from pisacov.iomod import conf as pconf
from pisacov.iomod.paths import check_path
import logging

_surl = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz'
_uniurl = 'https://www.uniprot.org/uniprot/'

def _default_values(key):
    defaultvals = {}
    defaultvals['HHBLITS_PARAMETERS'] = [3, 0.001, 'inf', 50, 99]
    defaultvals['LOWER_THRESHOLD'] = [0.2, 0.7, 0.2]
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
               'LOWTHRESHOLD': "Format should be [float or 'inf' or '-inf', " +
                   "float or 'inf' or '-inf', float or 'inf' or '-inf'].",
               'UNICLUST_FASTA_PATH': "UNICLUST_FASTA_PATH file not found.",
               'NEIGHBOURS_MINDISTANCE': "NEIGHBOURS_MINDISTANCE should be an integer.",
               'REMOVE_INTRA_CONTACTS': "REMOVE_INTRA_CONTACTS should be a boolean."}
    errormsg = compmsg[key] + ' Please, update the configuration file using pisaconf.'

    if (key == 'SIFTS_PATH' or key == 'PISA_PATH' or
            key == 'HHBLITS_PATH' or key == 'DMP_PATH'):
        try:
            val = check_path(val, 'file')
        except Exception:
            logging.critical(errormsg)
            raise NameError()
    if (key == 'UNICLUST_FASTA_PATH'):
        if val == "" or val is None:
            val = None
        else:
            try:
                val = check_path(val, 'file')
            except Exception:
                logging.critical(errormsg)
                raise NameError()
    elif (key == 'HHBLITS_DATABASE_DIR'):
        try:
            val = check_path(val, 'dir')
        except Exception:
            logging.critical(errormsg)
            raise NameError()
    elif (key == 'HHBLITS_DATABASE_NAME'):
        if isinstance(val, str) is False:
            try:
                val = str(val)
            except Exception:
                logging.critical(errormsg)
                raise ValueError()
        else:
            pass
    elif (key == 'HHBLITS_PARAMETERS'):
        if val is None or val == "":
            val = _default_values(key)
        else:
            if isinstance(val, list) is False:
                logging.critical('HHBLITS_PARAMETERS is not a python list.')
                raise ValueError()
            else:
                for n in range(len(val)):
                    if n == 1:
                        if isinstance(val[n], float) is False:
                            try:
                                val[n] = float(val[n])
                            except Exception:
                                logging.critical(errormsg)
                                raise TypeError
                    else:
                        if n == 2 and val[n] == 'inf':
                            pass
                        else:
                            if isinstance(val[n], int) is False:
                                try:
                                    val[n] = int(val[n])
                                except Exception:
                                    logging.critical(errormsg)
                                    raise TypeError
    elif (key == 'LOWTHRESHOLD'):
        if val is None or val == "":
            val = _default_values(key)
        else:
            if isinstance(val, list) is False:
                logging.critical('LOWTHRESHOLD is not a python list.')
                raise ValueError()
            else:
                for n in range(len(val)):
                    if isinstance(val[n], float) is False:
                        if val[n] == 'inf' or val[n] == '-inf':
                            pass
                        else:
                            try:
                                val[n] = float(val[n])
                            except Exception:
                                logging.critical(errormsg)
                                raise TypeError
    elif (key == 'NEIGHBOURS_MINDISTANCE'):
        if val is None or val == "":
            val = _default_values(key)
        else:
            if isinstance(val[n], int) is False:
                try:
                    val = int(val)
                except Exception:
                    logging.critical(errormsg)
                    raise TypeError
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
    except Exception:
        logging.critical(compmsg + trymsg)
        raise ValueError

    try:
        configfile['UNICLUST_FASTA_PATH'] = pconf.UNICLUST_FASTA_PATH
    except Exception:
        configfile['UNICLUST_FASTA_PATH'] = None

    try:
        configfile['HHBLITS_PARAMETERS'] = pconf.HHBLITS_PARAMETERS
    except Exception:
        configfile['HHBLITS_PARAMETERS'] = None

    try:
        configfile['LOWTHRESHOLD'] = pconf.HHBLITS_PARAMETERS
    except Exception:
        configfile['LOWTHRESHOLD'] = None

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
        except Exception:
            logging.critical('One or more of HHblits arguments given are not valid')
            raise ValueError

        outparams = [str(int(float(paramlist[0]))), str(float(paramlist[1])),
                     paramlist[2], str(float(paramlist[4])),
                     str(float(paramlist[5]))]
        if paramlist[2] != 'inf':
            outparams = str(int(float(paramlist[2])))

    return outparams


def _check_lowth(paramlist):
    """Return a list with the validated scores' lower thresholds.

    :param paramlist: User-provided list of thresholds.
    :type paramlist: list of str
    :raises ValueError: Arguments are not valid.
    :return: Complete and checked list of thresholds.
    :rtype: list of (float, str)
    """
    if not isinstance(paramlist, list):
        logging.critical('A list is required as input.')
    else:
        if len(paramlist) != 3:
            logging.critical('3 parameters are required as lower thresholds.')

    outparams = []
    for n in range(len(paramlist)):
        try:
            float(paramlist[n])
        except Exception:
            logging.critical('One or more of thresholds given are not valid')
            raise

        if paramlist[n] != "inf" and paramlist[n] != "-inf":
            outparams.append(float(paramlist[n]))
        else:
            outparams.append(paramlist[n])

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
    except Exception:
        logging.critical('Uniprot threshold given not valid.')
        raise ValueError

    try:
        dbpath = check_path(pconf.UNICLUST_FASTA_PATH, 'file')
    except Exception:
        logging.warning('Uniprot database file does not exist. Switching to server-only.')
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
        sources = ["psicov", "ccmpred", "dmp"]

    return sources


def _sources(lowth=[0.2, 0.7, 0.2]):
    """Return the subdir name and extension of each of the contact prediction types.

    :param lowth: Low threshold for each data source, defaults to [0.2, 0.7, 0.2].
    :type lowth: list [float], optional

    :return: Contact prediction types and location.
    :rtype: dict [list [str]]

    """
    sources = _sourcenames()
    confiledir = ["deepmetapsicov", "deepmetapsicov", "deepmetapsicov"]
    confilesuffix = ["psicov", "ccmpred", "deepmetapsicov.con"]
    conkittype = ["psicov", "ccmpred", "psicov"]

    outsinfo = {}
    for n in range(len(sources)):
        outsinfo[sources[n]] = [confiledir[n], confilesuffix[n],
                                conkittype[n], lowth[n]]

    return outsinfo
