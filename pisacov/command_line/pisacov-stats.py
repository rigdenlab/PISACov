"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import argparse
import time

logger = None

def create_argument_parser():
    """Create a parser for the command line arguments used in crops-renumber"""

    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)
    # MUTUALLY EXCLUSIVE: PDBID or FASTA path + STR path
    main_args = parser.add_mutually_exclusive_group(required=True)
    main_args.add_argument("-i", "--initialise", nargs=2, metavar="Initial_files",
                           help="Input sequence and structure filepaths.")
    main_args.add_argument("-s", "--skip_conpred", nargs=2, metavar="Initial_files",
                           help="If HHBLITS and DMP files have already been generated in pdbid/deepmetapsicov, they will be read and those processeses bypassed. Input sequence and structure filepaths.")

    parser.add_argument("-a","--add_noncropped", action='store_true', default=False,
                        help="Include results for original sequence even if cropped version exists..")

    parser.add_argument("-o","--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence. If -s option is used, be aware that this directory must already contain the pdbid/deepmetapsicov directory and its files.")
    parser.add_argument("-c","--collection_file", nargs=1, metavar="Collection_File",
                        help="Path to CSV file where pisacov signals will be appended. Default: oudir/pisacov_data.csv.")

    parser.add_argument("-u","--uniprot_threshold", nargs=1, metavar=("Uniprot_ratio_threshold"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. [MIN,MAX)=[0,100).')

    parser.add_argument("-h","--hhparams", nargs=5, metavar=("hhblits_new_parameters"),
                          help='Override default HHBLITS parameters in config file by introducing new ones: #iterations, E-value cutoff, Non-redundant seqs to keep, MinimumCoverageWithMasterSeq(%), MaxPairwiseSequenceIdentity.')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser

def main():
    starttime=time.time()
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    logger.info(pcl.welcome())



    return


if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        logger.info(pcl.ok())
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
