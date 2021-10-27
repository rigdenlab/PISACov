"""
This is PISACov, a PISA extension to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov.about import __prog__, __description__, __version__
from pisacov.about import  __author__, __date__, __copyright__

import os
import crops
import sys
import logging

def runcrops(seqin, strin, dbin, thin, upin, outdirin):#, loggingfile):
    """
    Run CROPS to produce clean sequence and structure files.

    :param seqin: Sequence filepath.
    :type seqin: str
    :param strin: Structure filepath.
    :type strin: str
    :param dbin: SIFTS database filepath.
    :type dbin: str
    :param outdirin: Output directory's path.
    :type outdirin: str

    """

    cropsdir = os.path.dirname(crops.__file__)
    pythonexec = '"'+sys.executable+'"'

    #RENUMBER PDB (ORIGINAL POSITIONS)
    logging.info('    Calling crops-renumber...')
    cropspy = os.path.join(cropsdir, 'crops', 'command_line', 'crops-renumber.py')
    try:
        os.system(pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + strin + ' ' +
                  '-o ' + outdirin)# + ' > ' + loggingfile)
    except:
        logging.critical('        An error was produced while executing Crops-Renumber')

    logging.info('    Done.')

    #CROP SEQUENCE
    logging.info('    Calling crops-cropseq...')
    cropspy = os.path.join(cropsdir, 'crops', 'command_line', 'crops-cropseq.py')
    try:
        os.system(pythonexec + ' ' + cropspy + ' ' + seqin + ' ' + dbin + ' ' +
                  '-o ' + outdirin)# + ' > ' + loggingfile)
    except:
        logging.critical('        An error was produced while executing Crops-Renumber')



def runhhblits(isitbio, outdir=output_dir(), param=None, spath=SEQUENCE_PATH):
    """
    Run HHSUITE to produce a Multiple Sequence Alignment from a single fasta sequence.

    :param isitbio: Specify whether using the ".bio" label in files or not.
    :type isitbio: bool
    :param outdir: Directory where results of HHSUITE will be printed out, defaults to output_dir()
    :type outdir: str, optional
    :param param: HHBLITS SEQUENCE ALINGMENT PARAMETERS, defaults to None (i.e. [2 , 0.001])
    :type param: list (int, float), dim=2, optional
    :param spath: Sequence path, defaults to SEQUENCE_PATH
    :type spath: str, optional
    :return: ConKit parsed MSA
    :rtype: MSA #CHANGE to format

    """

    # Obtain Multiple Sequence Alignment with HHblits
    msaa3mfile= pdbid()+biofile(isitbio)+".msa.a3m"
    msaa3mpath = os.path.join(outdir, msaa3mfile)

    hhsuite_exec='"'+HHSUITE_PATH+'"'
    if param is None or param == [ 2 , 0.001, 1000, 0, 90 ]:
        #os.system(hhsuite_exec + ' -i '+ spath + ' -d ' + HHBLITS_DATABASE_DIR + HHBLITS_DATABASE_NAME +' -oa3m ' + msaa3mpath + sysfileout())
        try:
            stdoutput=subprocess.check_output(hhsuite_exec + ' -i '+ spath + ' -d ' + HHBLITS_DATABASE_DIR + HHBLITS_DATABASE_NAME +' -oa3m ' + msaa3mpath, stderr=subprocess.STDOUT,shell=True)
        except:
            printout(stdoutput)
            interrupt("ERROR: An error occurred during the execution of HHSUITE.")
    else:
        #os.system(hhsuite_exec + ' -i '+ spath + ' -d ' + HHBLITS_DATABASE_DIR + HHBLITS_DATABASE_NAME + ' -n ' + str(param[0]) + ' -e ' + str(param[1]) + ' -oa3m ' + msaa3mpath + sysfileout())
        try:
            stdoutput=subprocess.check_output(hhsuite_exec + ' -i '+ spath + ' -d ' + HHBLITS_DATABASE_DIR + HHBLITS_DATABASE_NAME + ' -n ' + str(param[0]) + ' -e ' + str(param[1]) + ' -oa3m ' + msaa3mpath + ' -diff ' + str(param[2]) + ' -cov ' + str(param[3]) + ' -id ' + str(param[4]), stderr=subprocess.STDOUT, shell=True)
        except:
            printout(stdoutput)
            interrupt("ERROR: An error occurred during the execution of HHSUITE.")
    printout(stdoutput)
#   with open(stdoutpath(), 'a') as out:
#       out.write(stdoutput)

    # Convert A3M MSA file to Jones format (DMP standard input format)
    parsedmsa, msajonespath = msafilesgen(msaa3mpath)

    return parsedmsa, msajonespath

def msafilesgen(inpath_a3m):
    """
    Convert a3m format alignment to "jones" format (.aln) and print msa coverage graph

    :param inpath_a3m: Multiple Sequence Alignment (.a3m) path
    :type inpath_a3m: str
    :return: ConKit parsed MSA
    :rtype: MSA (Conkit) #CHANGE

    """
    outpath_aln = os.path.join(os.path.splitext(inpath_a3m)[0]+".aln")

    biotag = os.path.splitext(os.path.splitext(inpath_a3m)[0])[1]

    msacoveragefile = pdbid()+biotag+".msa.coverage.png" if biotag == ".bio" else pdbid()+".msa.coverage.png"
    msacoveragepath = os.path.join(output_dir(), msacoveragefile)

    parsemsa=ckio.read(inpath_a3m,'a3m')
    ckio.write(outpath_aln,'jones',parsemsa)

    ckplot.SequenceCoverageFigure(parsemsa, file_name=msacoveragepath)

    return parsemsa, outpath_aln