"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov.about import __prog__, __description__, __version__
from pisacov.about import  __author__, __date__, __copyright__

import argparse

from pisacov import command_line as pcl
#from pisacov import io as pio

def create_argument_parser():
    """Customise PISACov's configuration options"""

    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)

    parser.add_argument("-s", "--sifts", nargs=1, metavar="SIFTS_database_path",
                        help="Path to the SIFTS pdb_chain_uniprot.csv")

    # MUTUALLY EXCLUSIVE: PDBID or FASTA path + STR path
    parser.add_argument("input_seqpath",nargs=1, metavar="Sequence_filepath",
                        help="Input sequence filepath.")
    parser.add_argument("input_strpath", nargs=1, metavar="Structure_filepath",
                        help="Input structure filepath or dir. If a directory is inserted, it will act on all structure files in such directory.")
    parser.add_argument("input_database", nargs=1, metavar="Intervals_database",
                        help="Input intervals database filepath.")

    parser.add_argument("-o","--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")

    parser.add_argument("-s","--sort", nargs=1, metavar="Sort_type",
                        help="Sort output sequences in descending order by criteria provided - 'ncrops' or 'percent'. Add 'T' ('ncropsIN', 'percentIN') to ignore numbers from terminals. Only for multiple ID fasta inputs.")

    parser.add_argument("-u","--uniprot_threshold", nargs=2, metavar=("Uniprot_ratio_threshold","Sequence_database"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. [MIN,MAX)=[0,100) uniclust##_yyyy_mm_consensus.fasta-path')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser

def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = pcl.pisacov_logger(level="info")
    logger.info(pcl.welcome())


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