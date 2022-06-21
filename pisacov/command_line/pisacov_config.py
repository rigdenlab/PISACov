"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__
__script__ = 'PISACov Configuration script'

from pisacov import command_line as pcl
from pisacov.iomod import _conf_ops as pco
from pisacov.iomod import conf as pconf
from pisacov.iomod import online as pol
from pisacov.iomod import paths as ppaths

import argparse
import os

logger = None

def create_argument_parser():
    """Customise PISACov's configuration options"""

    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Configuration script of '+__prog__+' v.'+__version__+'\n'+__doc__)

    main_args = parser.add_mutually_exclusive_group(required=True)

    main_args.add_argument("-U", "--update_sifts_database", nargs=1, metavar="SIFTS_update",
                           help="Update the SIFTS database. Include path to save the SIFTS file pdb_chain_uniprot.csv.")
    main_args.add_argument("-C", "--conf_file",  action='store_true', default=False,
                           help="Update all options in the configuration file.")
    main_args.add_argument("-G", "--get_confpath", action='store_true', default=False,
                           help="If you prefer to edit the configuration file " +
                           "manually, this option will return the absolute path. " +
                           "Be aware that wrong editing the configuration file " +
                           "will affect the correct execution of PISACov.")
    main_args.add_argument("-V", "--view_configuration", action='store_true', default=False,
                           help="It reproduces the current contents of the configuration file.")
    main_args.add_argument("-s", "--sifts_path", nargs=1, metavar="New_sifts_path",
                           help="Update path to SIFTS database.")
    main_args.add_argument("-p", "--pisa_path", nargs=1, metavar="PISA_executable_path",
                           help="Update path to PISA executable.")
    main_args.add_argument("-b", "--hhblits_path", nargs=1, metavar="HHBLITS_executable_path",
                           help="Update path to HHBLITS executable.")
    main_args.add_argument("-l", "--hhblits_location", nargs=2, metavar=('HHBLITS_name', 'HHBLITS_dir'),
                           help="Update HHBLITS database's name and location.")
    main_args.add_argument("-d", "--dmp_path", nargs=1, metavar="DMP_executable_path",
                           help="Update path to DeepMetaPSICOV executable.")
    main_args.add_argument("-a", "--hhblits_arguments", nargs=5,
                           metavar=("n_iterations", "evalue_cutoff",
                                    "n_nonreduntant", "mincoverage",
                                    "maxpairwaise_seqidentity"),
                           help=("Introduce HHBLITS arguments." + os.linesep +
                           "DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99]. " +
                           "HHBlits DEEFAULT: [2, 0.001, 1000, 0, 90]"))
    main_args.add_argument("-r", "--reset_hhblits_arguments", action='store_true', default=False,
                           help="Reset to DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99].")
    main_args.add_argument("-i", "--intramolecular", action='store_true', default=None,
                           help="Pairs of residues that appear as intermolecular " +
                           "contacts and intramolecular contacts too are removed " +
                           "with NEIGHBOURS_MINDISTANCE set to 2; " +
                           "Intramolecular contacts determine this threshold.")
    main_args.add_argument("-n", "--neighbours", nargs=1, metavar="NEIGHBOURS_MINDISTANCE",
                           help="Set the minimum sequence distance for contacts " +
                           "to be considered. Intramolecular contacts will " +
                           "now not be considered.")
    main_args.add_argument("-u", "--uniclust_path", nargs='?', const="",
                           default=None, metavar="UNICLUST_fasta_path",
                           help="Update path to UniClust database fasta file.")

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser


def _outconffile(aconf):
    with open(pconf.__file__, 'w') as f:
        f.write("## Full path to pdb_chain_uniprot.csv SIFTS database.\n")
        f.write("SIFTS_PATH = '"+ aconf["SIFTS_PATH"] + "'\n")
        f.write("\n")
        f.write("## Uncomment if Uniprot fasta database is to be read locally.\n")
        f.write("## Full path to Uniprot fasta database.\n")
        if aconf['UNICLUST_FASTA_PATH'] is None:
            f.write("#UNICLUST_FASTA_PATH = ''")
        else:
            f.write("UNICLUST_FASTA_PATH = '" + aconf['UNICLUST_FASTA_PATH']+"'")
        f.write("\n")
        f.write("## Full path to PISA script\n")
        f.write("PISA_PATH = '"+ aconf["PISA_PATH"] + "'\n")
        f.write("\n")
        f.write("## NOTE: HHBLITS paths below should be the same provided in DMP_run.sh\n")
        f.write("## Full path to hhblits script\n")
        f.write("HHBLITS_PATH = '"+ aconf["HHBLITS_PATH"] + "'\n")
        f.write("## Name of sequence database to be used by HHBLITS (as requested by DeepMetaPSICOV)\n")
        f.write("HHBLITS_DATABASE_NAME = '"+ aconf["HHBLITS_DATABASE_NAME"] +
                "'  # e.g. uniclust30_2018_08\n")
        f.write("## Full path to directory containing sequence database to be used by HHBLITS\n")
        f.write("HHBLITS_DATABASE_DIR = '"+ aconf["HHBLITS_DATABASE_DIR"] + "'\n")
        f.write("\n")
        f.write("## Path to DeepMetaPsicov script\n")
        f.write("DMP_PATH = '" + aconf["DMP_PATH"] + "'\n")
        f.write("\n")
        f.write("## Input parameters for HHBLITS search.\n")
        f.write("## (#iterations, E-value cutoff, Non-redundant seqs to keep, MinimumCoverageWithMasterSeq(%%),MaxPairwiseSequenceIdentity)\n")
        f.write("## DeepMetaPSICOV & PISACOV DEFAULT: [3, 0.001, 'inf', 50, 99]\n")
        f.write("## HHBLITS DEFAULT: [2, 0.001, 1000, 0, 90]\n")
        f.write("## Uncomment next line to use non-DeepMetaPSICOV-default values\n")
        if aconf['HHBLITS_PARAMETERS'] is None:
            f.write("#HHBLITS_PARAMETERS = [3, 0.001, 'inf', 50, 99]\n")
        else:
            f.write("HHBLITS_PARAMETERS = " + str(aconf['HHBLITS_PARAMETERS']))
        f.write("\n")
        f.write("## CONTACT PARAMETERS")
        f.write("\n")
        f.write("## By default, pairs of residues that appear as intermolecular contacts\n")
        f.write("## and intramolecular contacts too are removed with NEIGHBOURS_MINDISTANCE set to 2;\n")
        f.write("## Intramolecular contacts determine this threshold\n")
        f.write("\n")
        f.write("## Minimum distance within the sequence for neighbours to be accounted for (int)\n")
        f.write("## |pos_res1-pos_res2|>=MinDistance; Minimum value: 2; Default: 2\n")
        f.write("## Uncomment next line if you want to bypass the default behaviour and fix a custom minimum distance while ignoring intramolecular contacts.\n")
        if aconf['REMOVE_INTRA_CONTACTS'] is True:
            f.write("#NEIGHBOURS_MINDISTANCE = 2\n")
        else:
            f.write("NEIGHBOURS_MINDISTANCE = " + str(aconf['NEIGHBOURS_MINDISTANCE']) +"\n")
        f.write("\n")

    return


def _lineformat(string):

    if string.startswith('#'):
        string = '\u001b[34m' + string + '\u001b[0m'
    else:
        stringlist = string.split('#')
        if len(stringlist) == 1:
            string = '\u001b[36m' + string + '\u001b[0m'
        else:
            string = '\u001b[36m' + stringlist[0]
            for n in range(1, len(stringlist)):
                string += '\u001b[34m' + '#' + stringlist[n]
            string += '\u001b[0m'

    return string


def main():

    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    welcomemsg, starttime = pcl.welcome(command=__script__)
    logger.info(welcomemsg)

    kwlist = pco._default_keys()

    configin = {}

    if args.update_sifts_database is not None:
        outpath = ppaths.check_path(args.update_sifts_database[0], 'either')
        if os.path.isdir(outpath) is True:
            outpath = os.path.join(outpath, ('pdb_chain_uniprot' +
                                             os.extsep + 'csv'))
        pol.getsifts(outpath)
        configin['SIFTS_PATH'] = outpath
    else:
        pass

    if args.view_configuration is True:
        print('** ' + __prog__+" v."+__version__+' configuration file **' + os.linesep)
        with open(pconf.__file__) as f:
            for line in f:
                line = _lineformat(line)
                print(line, end='')
        print('** End of configuration file **' + os.linesep)
        return
    else:
        pass

    if args.get_confpath is True:
        print(pconf.__file__)
        return
    else:
        pass

    if args.conf_file is False:
        newconf = pco._parse_conf()
    else:
        newconf = {}

    if args.reset_hhblits_arguments is True:
        newconf['HHBLITS_PARAMETERS'] = pco._default_values('HHBLITS_PARAMETERS')
    else:
        pass

    if args.sifts_path is not None:
        configin['SIFTS_PATH'] = args.sifts_path[0]
    if args.pisa_path is not None:
        configin['PISA_PATH'] = args.pisa_path[0]
    if args.hhblits_path is not None:
        configin['HHBLITS_PATH'] = args.hhblits_path[0]
    if args.dmp_path is not None:
        configin['DMP_PATH'] = args.dmp_path[0]
    if args.uniclust_path is not None:
        configin['UNICLUST_FASTA_PATH'] = args.uniclust_path[0]
    if args.hhblits_arguments is not None:
        configin['HHBLITS_PARAMETERS'] = args.hhblits_arguments[0]
    if args.neighbours is not None:
        configin['NEIGHBOURS_MINDISTANCE'] = args.neighbours[0]
        configin['REMOVE_INTRA_CONTACTS'] = False
    if isinstance(args.hhblits_location, list) is True:
        configin['HHBLITS_DATABASE_NAME'] = args.hhblits_location[0]
        configin['HHBLITS_DATABASE_DIR'] = args.hhblits_location[1]


    compmsg = {'SIFTS_PATH': "Please, enter the path to the SIFTS csv file:\n",
               'PISA_PATH': "Please, enter the path to the PISA executable:\n",
               'HHBLITS_PATH': "Please, enter the path to the HHBLITS executable:\n",
               'HHBLITS_DATABASE_NAME': ("Please, enter the name of the HHBLITS database name\n" +
                                         "as requested by DeepMetaPSICOV (e.g. uniclust30_2018_08):\n"),
               'HHBLITS_DATABASE_DIR': "Please, enter the directory of the HHBLITS database:\n",
               'DMP_PATH': "Please, enter the path to the DeepMetaPSICOV executable:\n",
               'HHBLITS_PARAMETERS': ("Please, enter the 5 HHBLITS parameters.\n" +
                                      "#iterations, E-value cutoff, Non-redundant seqs to keep, " +
                                      " MinimumCoverageWithMasterSeq(%), and MaxPairwiseSequenceIdentity.\n" +
                                      "Leave empty for DMP default (3, 0.001, 'inf', 50, 99).\n"),
               'UNICLUST_FASTA_PATH': ("Please, enter the path to the UniClust fasta file.\n"
                                       "An empty input will deactivate this option.\n"),
               'NEIGHBOURS_MINDISTANCE': ("Press ENTER for intramolecular contacts " +
                                          "to be removed from intermolecular contact lists.\n" +
                                          "Otherwise, enter the minimum distance to be cosidered "+
                                          "between neighbours. This will override the removal of" +
                                          "intramolecular contacts that will be now ignored.\n")}

    for keystr in kwlist[:-1]:
        if keystr in configin or args.conf_file is True:
            if args.conf_file is True:
                while True:
                    newval = input(compmsg[keystr])
                    try:
                        newconf[keystr] = pco._check_input(newval, keystr)
                        if keystr == 'NEIGHBOURS_MINDISTANCE':
                            if newval is None or newval == "":
                                newconf['REMOVE_INTRA_CONTACTS'] = True
                            else:
                                newconf['REMOVE_INTRA_CONTACTS'] = False
                    except Exception:
                        print("Not a valid input. Please, try again:")
                    else:
                        break
            else:
                newconf[keystr] = pco._check_input(configin[keystr], keystr)
                if keystr == 'NEIGHBOURS_MINDISTANCE':
                    newconf['REMOVE_INTRA_CONTACTS'] = pco._check_input(configin['REMOVE_INTRA_CONTACTS'],
                                                            'REMOVE_INTRA_CONTACTS')
                else:
                    pass

    if 'REMOVE_INTRA_CONTACTS' in configin:
        newconf['REMOVE_INTRA_CONTACTS'] = pco._check_input(configin['REMOVE_INTRA_CONTACTS'],
                                                            'REMOVE_INTRA_CONTACTS')
        if newconf['REMOVE_INTRA_CONTACTS'] is True:
            newconf['NEIGHBOURS_MINDISTANCE'] == 2
        else:
            pass
    else:
        pass

    _outconffile(newconf)

    endmsg = pcl.ok(starttime, command=__script__)
    logger.info(endmsg)
    return

if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        # logger.info(pcl.ok())
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
