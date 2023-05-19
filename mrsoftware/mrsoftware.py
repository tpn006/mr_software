#!/usr/bin/ python3

import argparse
#from . import myutils as myutils
#from mrsoftware import __version__
import os
from pyfaidx import Fasta
import sys
import numpy as np
import math
import pandas as pd
import scipy
import argparse


genome = Fasta("../tests/GRCm38.fa")
bed = None
meme = None
nucs = {"A": 0, "C": 1, "G": 2, "T": 3}

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

def main():
    genome_path = "../tests/GRCm38.fa"
    bed_path = "../tests/Oct4.peaks.bed"
    meme_path ="../tests/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme" 
    print("Hello wolr")
    parser = argparse.ArgumentParser(
		prog="mrsoftware",
		description="Command-line script to perform  motif recognition from BED file"
	)

	# Input files
        # input Bed File
    parser.add_argument("bed", help="BED file", type=str)
        # input meme file of known motifs  
    parser.add_argument("-m", "--meme", help="meme file corresponding to the transcription factor", type=str, metavar="FILE", required=True)
        # input the reference genome to use
    parser.add_argument("-g","--genome", help="Path to the reference genome that peaks will be pulled from", type=str, metavar="FILE", required=True)
    
	# Output file to write to
    parser.add_argument("-o", "--out", help="Write output to file. Default: stdout", metavar="FILE", type=str, required=False)

	# Other options
    #parser.add_argument("--version", help="Print the version and quit", action="version", version = '{version}'.format(version=__version__))
    
    # Parse args
    args = parser.parse_args()

	# Set up output file
    if args.out is None:
        outf = sys.stdout
    else: outf = open(args.out, "w")

    # Load  genome FASTA
    if args.genome is not None:
        if not os.path.exists(args.genome):
            ERROR("{fasta} does not exist".format(fasta=args.genome))
        genome = Fasta(args.genome)
    else:
        print(":(")
        genome = None


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
    # Fill in pfm below. Note: pfm[i,j] can be used to 
    # access and set the value for position j at nucleotide i
    # your code here
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


def GetPWM(binding_sites, background_freqs=[0.25, 0.25, 0.25, 0.25]):
    """ Compute the PWM for a set of binding sites
    
    Parameters
    ----------
    binding_sites : list of str
        List of sequences 
    background_freqs: list of float
        Background frequency of A, C, G, T
    
    Returns
    -------
        pwm : 2d np.array
        
    Assumes all sequences have the same length
    """
    pwm = np.zeros((4, len(binding_sites[0])))
    pfm = GetPFM(binding_sites)
    pfm = pfm + 0.01 # Add pseudocount. Don't change this!
    # Compute pwm below
    # Note: np.sum(pfm[:,j]) will give the sum of counts for column j
    # Note: pfm[i,j]/np.sum(pfm[:,j]) gives p(i,j) (frequency of nucleotide i at position j)
    # your code here
    j = 0
    while j < len(binding_sites[0]):
        i = 0
        while i < 4:
            pwm[i,j] = math.log2((pfm[i,j]/np.sum(pfm[:,j]))/background_freqs[i])
            i+=1
        j+=1        
    return pwm

def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    # your code here
    i = 0
    while i < len(sequence):
        if sequence[i] == 'A':
            score += pwm[0, i]
        elif sequence[i] == 'C':
            score += pwm[1, i]
        elif sequence[i] == 'G':
            score += pwm[2, i]
        elif sequence[i] == 'T':
            score += pwm[3, i]
        i = i+1
    return score


def GetThreshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       % of null_dist that should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 # set this below to be the score threshold to obtain a p-value <0.01
    # your code here
    num = len(null_dist) * pval
    temp_dist = null_dist.copy()
    i = 1
    while i < num:
        temp_dist.pop(temp_dist.index(max(temp_dist)))
        i+=1
    thresh = max(temp_dist)
    return thresh


def ScanSequence(pwm, sequence):
    """ Scan a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
        PWM matrix
    sequence : str
        Long nucleotide string
        
    Returns
    -------
    scores : list of float
        scores[i] should give the score of the substring sequence[i:i+n]
    """
    n = pwm.shape[1]
    scores = [0]*(len(sequence)-n+1) # list of scores. scores[i] should give the score of the substring sequence[i:i+n]
    # your code here
    i = 0
    while i < (len(sequence)-n+1):
        scores[i] = ScoreSeq(pwm, sequence[i:i+n])
        i+=1
    return scores


def ReverseComplement(sequence):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = ""
    # your code here
    i = len(sequence)-1
    while i >= 0:
        if sequence[i] == 'A':
            revcomp = revcomp + 'T'
        elif sequence[i] == 'C':
            revcomp = revcomp + 'G'
        elif sequence[i] == 'G':
            revcomp = revcomp + 'C'
        elif sequence[i] == 'T':
            revcomp = revcomp + 'A'
        i = i-1
    return revcomp


def ComputeNucFreqs(sequences, background_freqs = [0.25,0.25,0.25,0.25]):
    """ Compute nucleotide frequencies of a list of sequences
    
    Parameters
    ----------
    sequences : list of str
       List of sequences
       
    Returns
    -------
    freqs : list of float
       Frequencies of A, C, G, T in the sequences
    """
    freqs = [0.25, 0.25, 0.25, 0.25] # compute frequency of A, C, G, T
    # your code here
    total = 0
    A_count = 0
    C_count = 0
    G_count = 0
    T_count = 0
    for seq in sequences:
        total += len(seq)
        i = 0
        while i < len(seq):
            if seq[i] == 'A':
                A_count = A_count + 1
            elif seq[i] == 'C':
                C_count = C_count + 1
            elif seq[i] == 'G':
                G_count = G_count + 1
            elif seq[i] == 'T':
                T_count = T_count + 1
            i+=1
    freqs[0] = A_count/total
    freqs[1] = C_count/total
    freqs[2] = G_count/total
    freqs[3] = T_count/total
    return freqs


def getSeqFromGenome(chromosome, start, end):
    """
    Returns sequence of peak from genome

    Parameters:
    ------------
    chromosome : string (what chromosome number)
    start: int (start of peak)
    end: int (end of peak)

    Returns: string of nucleotides
    """
    #genes = Fasta(genome)
    return genome[chromosome][start:end]


def getPeaksFromBed(bed_path):
    """
    Returns a list of locations of peaks in the format {chr, start, stop} from a bed file containing peaks
    """
    peaks = []
    file = open(bed_path)
    for line in file.readlines():
        if '#' not in line[0:0]: 
            peaks.append(line.split('\t')[0:3])
    return peaks


def listOfPeakSeq(peak_info):
    """
    Returns a list of all the peak sequences
    """
    peak_list = []
    ent_list = getPeaksFromBed(bed)
    for ent in ent_list:
        peak_list.append(getSeqFromGenome(ent[0], ent[1], ent[2]))
    return peak_list

#may not be necessary
def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ Compute fisher exact test to test whether motif enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    pval = -1
    # your code here
    table = [[peak_motif, peak_total-peak_motif], [bg_motif, bg_total-bg_motif]]
    odds, pval = scipy.stats.fisher_exact(table)
    return pval


main()

print(getSeqFromGenome(17, 50427861, 50427936))

getPeaksFromBed("../tests/Oct4.peaks.bed")    