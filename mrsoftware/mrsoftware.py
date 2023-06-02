#!/usr/bin/ python3
#TODO EVERYTIME THERE ARE ACTUAL CHANGES, RERUN THE INSTALL COMMAND TO MAKE IT ACTUALLY WORK AS A UNIX COMMAND. ALSO CHANGE THE IMPORT TO "from . "
import argparse
#from . import KnownMotif
import KnownMotif
import os
from pyfaidx import Fasta
import sys
import numpy as np
import math
import argparse
import random


genome = None
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

def sysMessage(message):
    sys.stdout.write(bcolors.OKCYAN + message)

def writeMessage(str, output_file = sys.stdout):
    output_file.write(str)

def main():
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

    #Check file integrity
    sysMessage("Initializing by checking argument validity\n")
    # genome
    if args.genome is not None:
        if not os.path.exists(args.genome):
            ERROR("{genome} does not exist".format(genome=args.genome))
        genome_path = args.genome
        global genome
        genome = Fasta(genome_path)
    else:
        ERROR("Path to Genome not provided :O")
        print(":(")
        return # Just quit and give up
    
    # meme 
    if args.meme is not None:
        if not os.path.exists(args.meme):
            ERROR("{meme} does not exist".format(meme=args.meme) )
        meme_path = args.meme
    else:
        ERROR("Path to meme file not provided >:(")
        return # Just quit and give up
    
    # bed
    if args.bed is not None:
        if not os.path.exists(args.bed):
            ERROR("{bed} does not exist".format(meme=args.bed) )
        bed_path = args.bed
    else:
        ERROR("Path to bed file not provided :(")
        return # Just quit and give up
    
    out_file = None
    if args.out is not None:
        out_file = open(args.out, 'w')
    else:
        out_file = sys.stdout

    #~~~~~~~~~~~~~~~~~~~~~~~ BEGIN PIPELINE ~~~~~~~~~~~~~~~~~~~~~~~~#
    #get the sequences of all the peaks
    print("Getting peak sequences")
    peak_seq_list = getListOfPeakSeqFromBed(bed_path)
    writeMessage("There are " + str(len(peak_seq_list)) + " peaks. \n", output_file=out_file)
    print("Creating random sequences from the genome")
    # hardcoding in chr17 and length of 17
    random_seqs_list = generateRandomSeqsFromChr("chr17", len(peak_seq_list), 75) #generateListOfRandomSeqsFromPeaks(getPeaksFromBed(bed_path))
    print("We are using " + str(len(random_seqs_list)) + " random sequences")

    # Get motifs
    print("Calculating random backgrounds")
    backgroundFreqs = calculateBackgroundFrequencies(random_seqs_list)
    print("Background frequency is: " + str(backgroundFreqs))
    print("Motifs from the meme file")
    known_motifs = makeKnownMotifObjs(meme_path, backgroundFreqs)
    print("We are using " + str(len(known_motifs)) + " motifs")
    for motif in known_motifs:
        motif.magicDoAllFunction(peak_seq_list, random_seqs_list) #MAKE THE PWM SOON
    '''

    writeMessage("There are " + str(len(known_motifs)) + " known motifs. \n", output_file=out_file)

    # Score motifs based on peaks
    #printMotfs(known_motifs)
    scoreAllRandomSequencesForMotifs(random_seqs_list, known_motifs)
    scoreAllSequencesForMotifs(peak_seq_list, known_motifs)
    setAllThresholds(known_motifs, random_seqs_list, 0.001)
    '''

    #known_motifs.sort(key = KnownMotif.KnownMotif.getThresh)
    for motif in known_motifs:
        writeMessage(str(motif) + "\n", output_file=out_file)
    
    #returnMotifsWithLowestPval(known_motifs, 4)
    sysMessage("DONE >:[ \n\n")

def generateRandomSeqsFromChr(chr, num, length):
    """
    Get <num> random sequences from chromosome <chr>
    """
    backgroundSeqs = []
    for i in range(num):
        backgroundSeqs.append(getRandomBackgroundSeq(chr, length))
    return backgroundSeqs


def returnMotifsWithLowestPval(motif_array, num_of_motifs):
    """
    will return the top motifs with the lowest p values

    Parameters:
    -----------
        motif_array: list of arrays
        num_of_motifs: number of top motifs to return
    
    Returns:
    --------
        list of top motifs (length num_of_motifs)

    """
    sorted_motifs = motif_array[0]
    #motif_array.sort(key = KnownMotif.KnownMotif.getPval)
    i = 0
    while i < len(motif_array):
        j = 0
        while j < len(sorted_motifs):
            if motif_array[i].getPval < sorted_motifs[j].getPval:
                sorted_motifs.insert(j, motif_array[i])
                break
            j+=1
        if j == len(sorted_motifs):
            sorted_motifs.append(motif_array[i])
        i+=1
    return sorted_motifs[:num_of_motifs]


def scoreAllSequencesForMotifs(peak_seq_list, motif_array):
    """
    Given a list of sequences and motifs, score all motifs with the sequences

    Parameters:
    -----------
        peak_seq_list : List of sequences from peaks
        motif_array : array of all the motifs to be scored

    Returns:
    --------
        Nothing
    """
    for peak in peak_seq_list:
        scoreAllMotifsByPeak(peak, motif_array)
    
def scoreAllRandomSequencesForMotifs(peak_seq_list, motif_array):
    """
    Given a list of sequences and motifs, score all motifs with the sequences

    Parameters:
    -----------
        peak_seq_list : List of sequences from peaks
        motif_array : array of all the motifs to be scored

    Returns:
    --------
        Nothing
    """
    for peak in peak_seq_list:
        scoreAllRandomMotifsByPeak(peak, motif_array)

def scoreAllRandomMotifsByPeak(peak_seq, motif_array):
    """
    For a given sequence, score all motifs using that sequence. Updates all the motifs with new scores

    Parameters:
    -----------
        peak_seq : Sequence of the peak to score motifs with
        motif_array : array of all the known motif objects 

    Returns:
    --------
        No Return
    """
    for motif in motif_array:
        motif.scoreRandSequence(str(peak_seq))

def scoreAllMotifsByPeak(peak_seq, motif_array):
    """
    For a given sequence, score all motifs using that sequence. Updates all the motifs with new scores

    Parameters:
    -----------
        peak_seq : Sequence of the peak to score motifs with
        motif_array : array of all the known motif objects 

    Returns:
    --------
        No Return
    """
    for motif in motif_array:
        motif.scoreSequence(peak_seq)

def getSeqFromPeakList(peaks):
    """
    peaks - 
    returns the sequences of the peaks
    """
    peaks_seq_list = []
    for peak in peaks:
        peaks_seq_list.append(getSeqFromPeak(peak))
    return peaks_seq_list

def getSeqFromPeak(peak):
    """
    returns
    """
    return getSeqFromGenome(peak[0],peak[1],peak[2])

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
    # pfm[i,j]: value for position j at nucleotide i
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

def GetBackgroundFreqs(meme_path):
    """ Get background freqs from meme file
    
    Parameters
    ----------
        meme_path = meme file that was inputed

    Returns
    -------
        list of background freqs in order of A, C, G, T
    """

    background_freqs = []
    file = open(meme_path)
    prev_line = ''
    back_seq_line = ''
    for line in file.readlines():
        if 'Background letter frequencies' in line: 
            prev_line = line
        elif prev_line != '':
            back_seq_line = line
            break
    split_text = back_seq_line.split(" ")
    background_freqs.append(float(split_text[1].strip()))
    background_freqs.append(float(split_text[3].strip()))
    background_freqs.append(float(split_text[5].strip()))
    background_freqs.append(float(split_text[7].strip()))
    return background_freqs

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
       Null distribution of scores (random scores)
    pval : float
       % of null_dist that should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 # set this below to be the score threshold to obtain a p-value <0.01
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
    """ 
    Compute nucleotide frequencies of a list of sequences
    
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

    Returns: 
    string of nucleotides
    """
    #genes = Fasta(genome)
    if chromosome in genome: # if chromosome is in the same format as the indexing
        return str(genome[chromosome][start:end])
    elif ("chr"+str(chromosome)) in genome: # if you need to add "chr" to the chromosome number to get it to match
        return str(genome["chr"+str(chromosome)][start:end])
    elif str(chromosome)[3:] in genome: # if you need to remove "chr" from the chromosome name to get it to match the index
        return str(genome[str(chromosome)[3:]][start:end])
    else:
        print("CHROMOSOME NOT FOUND!!!!! >:[ ") # TODO CHANGE THIS to a REAL ERROR Rather than this print


def getPeaksFromBed(bed_path):
    """
    Helper function
    Returns a list of locations of peaks in the format {chr, start, stop} from a bed file containing peaks
    
    Description

    Parameters:
    -----------
        meme_path : path to the meme file

    Returns:
    --------
        known_motifs_list : list of known motifs
    
    """
    peaks = []
    file = open(bed_path)
    for line in file.readlines():
        if '#' not in line: 
            peaks.append(line.split('\t')[0:3])

    file.close()
    return peaks

def getListOfPeakSeqFromBed(bed_path):
    """
    Returns a list of all the peak sequences
    Calls get peaks from Bed

    Description

    Parameters:
    -----------
        meme_path : path to the meme file

    Returns:
    --------
        known_motifs_list : list of known motifs
    """
    peaks = []
    peak_list = []
    ent_list = getPeaksFromBed(bed_path)
    for ent in ent_list:
        temp_seq = getSeqFromGenome(int(ent[0]), int(ent[1]), int(ent[2]))
        peak_list.append(str(temp_seq)) # add the string of the sequence
    return peak_list

def makeKnownMotifObjs(meme_path, background_freqs):
    """
    makes list of KnownMotif objects from the meme file

    Parameters:
    -----------
        meme_path : path to the meme file

    Returns:
    --------
        known_motifs_list : list of known motifs
    """
    known_motifs_list = []
    file = open(meme_path)
    lines = file.readlines()
    i = 0
    while i < len(lines):
        if 'MOTIF' in lines[i]:
            pwm_index_start = 0
            if 'letter-probability matrix' not in lines[i+1]: # empty line in case motif doesnt start yet
                pwm_index_start = i+3
            else:
                pwm_index_start = i+2
            name = lines[i][6:].strip() # name of the motif
            alength = 4
            w_str = 'w= '
            w_index =  lines[pwm_index_start-1].index(w_str)
            nsites_index = lines[pwm_index_start-1].index(' nsites')
            w = int(lines[pwm_index_start-1][w_index + len(w_str) : nsites_index])
            j = pwm_index_start
            k = 0
            temp_pwm = [[0 for i in range(w)] for j in range(alength)]
            while j < len(lines):
                if 'URL' in lines[j] or len(lines[j]) == 0:
                    break
                temp_line = lines[j].split()
                if len(temp_line) == 0:
                    break
                pseudocount = 0.00001

                A_val = math.log2(float(temp_line[0].strip())/background_freqs[0]+pseudocount)
                C_val = math.log2(float(temp_line[1].strip())/background_freqs[1]+pseudocount)
                G_val = math.log2(float(temp_line[2].strip())/background_freqs[2]+pseudocount)
                T_val = math.log2(float(temp_line[3].strip())/background_freqs[3]+pseudocount)
                
                temp_pwm[0][k] = float(A_val)
                temp_pwm[1][k] = float(C_val)
                temp_pwm[2][k] = float(G_val)
                temp_pwm[3][k] = float(T_val)
                k+=1
                j+=1 
            i = j
            
            known_motifs_list.append(KnownMotif.KnownMotif(name, alength, w, temp_pwm))
        i+=1
    return known_motifs_list

def calculateBackgroundFrequencies(background_seqs):
    freqs = [0,0,0,0] # compute frequency of A, C, G, T
    counts = [0,0,0,0]
    for seq in background_seqs:
        temp_seq = str(seq)
        counts[0] += temp_seq.count("A")
        counts[1] += temp_seq.count("C")
        counts[2] += temp_seq.count("G")
        counts[3] += temp_seq.count("T")
    total =  sum(counts)
    for i in range(len(freqs)):
        freqs[i] = counts[i]/total
    return freqs

def getListOfReversePeakSeq(forward_seqs):
    reverse_seqs = []
    for seq in forward_seqs:
        reverse_seqs.append(ReverseComplement(seq))
    return reverse_seqs

def getRandomBackgroundSeq(chr, length):
    """
    Given a chromosome, get a random sequence of that length

    Parameters:
    -----------
        chr: chromosome number of random seq
        length: length of random seq that will be generated

    Returns:
    --------
        the sequence
    """
    length_of_chromosome = -1
    if chr in genome: # if chromosome is in the same format as the indexing
        length_of_chromosome = len(genome[chr])
    elif ("chr"+str(chr)) in genome: # if you need to add "chr" to the chromosome number to get it to match
        length_of_chromosome = len(str(genome["chr"+str(chr)]))
    elif str(chr)[3:] in genome: # if you need to remove "chr" from the chromosome name to get it to match the index
        length_of_chromosome = len(str(genome[str(chr)[3:]]))
    else:
        print("CHROMOSOME NOT FOUND!!!!! >:[ ") # TODO CHANGE THIS to a REAL ERROR Rather than this print
    start = random.randint(0, length_of_chromosome-length)
    return getSeqFromGenome(chr, start, start + length)

def generateListOfRandomSeqsFromPeaks(peaks_info):
    #print("running generateListOfRandomSeqsFromPeaks")
    #peaks_info from getPeaksFromBed 
    random_seqs = []
    for ent in peaks_info:
        length = int(ent[2])-int(ent[1])
        chr = ent[0]
        random_seqs.append(getRandomBackgroundSeq(chr,length))
    return random_seqs

def setAllThresholds(known_motifs, background_seqs, pval):
    for motif in known_motifs:
        motif.calcAndSetThreshold(background_seqs, pval)

def printMotfs(knownMotifs):
    for motif in knownMotifs:
        motif.printName()
        print(motif.printPWM())

# Run Main
if __name__ == "__main__":
    main()
