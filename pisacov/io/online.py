"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""
from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__
from pisacov.io import _surl

from urllib import request as request
from contextlib import closing
import gzip
import shutil

# DOWNLOAD AND EXTRACT SIFTS DATABASE
def getsifts(outfile):

    with closing(request.urlopen(_surl)) as handle:
        with gzip.open(handle, 'rb') as f_in:
            with open(outfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    return

# READ FROM UNIPROT ONLINE

# DOWNLOAD PDB STRUCTURE

# DOWNLOAD FASTA FILE
