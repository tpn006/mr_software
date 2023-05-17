"""
Utilities for mrsoftware
"""
import sys
import numpy as np

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg):
	"""
	Print an error message and die

	Parameters
	----------
	msg : str
	   Error message to print
	"""
	sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
	sys.exit(1)


nucs = {"A": 0, "C": 1, "G": 2, "T": 3} # this might be helpful

def GetPFM(sequences):
    """ Compute the PFM for a set of sequences
    
    Parameters
    ----------
    sequences : list of str
        List of sequences (e.g. binding_sites)
    
    Returns
    -------
        pfm : 2d np.array
        
    Assumes all sequences have the same length
    """
    pfm = np.zeros((4, len(sequences[0])))
    for seq in sequences:
        i = 0
        while i < len(seq):
            if seq[i] == 'A':
                pfm[0,i] = pfm[0,i] + 1
            elif seq[i] == 'C':
                pfm[1,i] = pfm[1,i] + 1
            elif seq[i] == 'G':
                pfm[2,i] = pfm[2,i] + 1
            elif seq[i] == 'T':
                pfm[3,i] = pfm[3,i] + 1
            i += 1
    return pfm