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
    sys.stdout.write(bcolors.OKCYAN + message + bcolors.ENDC)

def writeMessage(str, output_file = sys.stdout):
    output_file.write(str)

def main():
    parser = argparse.ArgumentParser(
		prog="mrsoftware",
		description="Command-line script to perform  motif recognition from BED file"
	)

	# Input files
        # input Bed File
    parser.add_argument("bed", help="path to .BED file with the peaks", type=str)
        # input meme file of known motifs  
    parser.add_argument("-m", "--meme", help="path to .meme file containing nucleotide probability matrices for motifs to test for enrichment", type=str, metavar="FILE", required=True)
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
    else: 
        outf = open(args.out, "w")

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
    writeMessage("Optimizing for provided .fa")
    needs_chr = determineIfNeedCHR(genome_path)
    writeMessage("Needs to add \"chr\" is " + str(needs_chr) + "\n")

    writeMessage("Getting peak sequences\n")
    writeMessage("Collecting random sequences from the genome\n")
    peak_seq_list, random_background_seqs_list = getListOfPeakSeqFromBed(bed_path, need_chr = needs_chr)
    writeMessage("There are " + str(len(peak_seq_list)) + " peaks. \n")
    writeMessage("We are using " + str(len(random_background_seqs_list)) + " random sequences\n")

    # Get motifs
    writeMessage("Calculating random backgrounds\n")
    backgroundFreqs = calculateBackgroundFrequencies(random_background_seqs_list)
    writeMessage("Background frequency is: " + str(backgroundFreqs) + "\n")
    writeMessage("Motifs from the meme file\n")
    known_motifs = makeKnownMotifObjs(meme_path, backgroundFreqs)
    writeMessage("We are using " + str(len(known_motifs)) + " motifs\n")
    writeMessage("Creating new randomly generated sequences\n")
    amount = len(peak_seq_list)
    length_of_random = len(peak_seq_list[0])
    random_generated_seqs_list = createArrOfRandomlyGeneratedSeqs(amount, backgroundFreqs, length_of_random)
    writeMessage("Done creating " + str(amount) + " random sequenceis of length " + str(length_of_random) +"\n")
    for motif in known_motifs:
        if motif != None:
            motif.scoreAllGetThresholdScoreCalcPValue(peak_seq_list, random_generated_seqs_list) #MAKE THE PWM SOON

    known_motifs.sort(key = KnownMotif.KnownMotif.getPval)
    for motif in known_motifs:
        writeMessage(str(motif) + "\n", output_file=out_file)
    #returnMotifsWithLowestPval(known_motifs, 4)
    sysMessage("DONE >:[ \n\n")
    out_file.close()

def determineIfNeedCHR(genome_path):
    """
    returns True if you need to add "chr" to the chromosome name 
    """
    file = open(genome_path + ".fai")
    if str.isdigit(file.readline().split("\t")[0]):
        return False
    return True

def generateRandomSeqsFromChr(chr, num, length):
    """
    Get <num> random sequences from chromosome <chr>
    
    Parameters:
    -----------
        chr: Chromosome to get the sequence from
        num: How many to get
        length:
    
    Returns:
    --------
        backgroundSeqs:
    """
    backgroundSeqs = []
    for i in range(num):
        toAdd = getRandomBackgroundSeq(chr, length)
        if toAdd != None:
            backgroundSeqs.append(toAdd)
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

def randomSeqFromFrequencies(frequencies, length):
    """
    Given the frequencies and desiered length, generate a random sequence to fit the request.
    """
    sequence = ""
    nuceleotides = ["A", "C", "G", "T"]
    for i in range(length):
        sequence += random.choices(nuceleotides,weights=frequencies)[0]
    return sequence

def createArrOfRandomlyGeneratedSeqs(num, frequencies,length):
    """
    Given a number of sequences to generate, the frequencies at which those appear, and the length of the sequence, return an array of sequences
    """
    return [randomSeqFromFrequencies(frequencies, length) for i in range(num)]

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
    #chromosome = str(chromosome).removeprefix("chr")
    if chromosome in genome: # if chromosome is in the same format as the indexing
        return str(genome[chromosome][start:end])
    elif ("chr"+str(chromosome)) in genome: # if you need to add "chr" to the chromosome number to get it to match
        return str(genome["chr"+str(chromosome)][start:end])
    else:
        print("CHROMOSOME NOT FOUND!!!!! in getSeqFromGenome >:[ Params were: " + str(chromosome) +" at "+ str(start) +":"+ str(end)) # TODO CHANGE THIS to a REAL ERROR Rather than this print

def getPeaksFromBed(bed_path):
    """
    Helper function
    Returns a list of locations of peaks in the format {chr, start, stop} from a bed file containing peaks
    
    Description

    Parameters:
    -----------
        bed_path : path to the bed file

    Returns:
    --------
        peaks: list of locations of peaks in the format {chr, start, stop}
    
    """
    peaks = []
    file = open(bed_path)
    for line in file.readlines():
        if '#' not in line: 
            peaks.append(line.split('\t')[0:3])

    file.close()
    return peaks

def getListOfPeakSeqFromBed(bed_path, need_chr = True):
    """
    Returns a list of all the peak sequences
    Calls getPeaksFromBed

    Parameters:
    -----------
        bed_path : path to the bed file

    Returns:
    --------
        peak_list : list of peaks (nucleotide sequence)
        random_seq_list : list of random ones from that genome
    """
    peak_list = []
    random_list = []
    ent_list = getPeaksFromBed(bed_path)
    has_chr = "chr" in ent_list[0][0]
    print("There are " + str(len(ent_list)) + " peaks")
    for ent in ent_list:
        if need_chr and has_chr: # can go ahead and access easily
            temp_seq = getSeqFromGenome(ent[0], int(ent[1]), int(ent[2]))
            temp_rand = getRandomBackgroundSeq(ent[0], int(ent[2])-int(ent[1]))
        elif need_chr and not has_chr: # needs chr but doesnt have chr already, must add chr
            temp_seq = getSeqFromGenome("chr"+ent[0], int(ent[1]), int(ent[2]))
            temp_rand = getRandomBackgroundSeq("chr"+ent[0], int(ent[2])-int(ent[1]))
        elif not need_chr and has_chr: # doesnt need chr but has chr, must remove chr
            temp_seq = getSeqFromGenome(ent[0][3:], int(ent[1]), int(ent[2]))
            temp_rand = getRandomBackgroundSeq(ent[0][3:], int(ent[2])-int(ent[1]))
        else: # doesnt need and doesnt have
            temp_seq = getSeqFromGenome(ent[0], int(ent[1]), int(ent[2]))
            temp_rand = getRandomBackgroundSeq(ent[0], int(ent[2])-int(ent[1]))
        if temp_seq != None: # will be None if the chromosome isnt available
            peak_list.append(str(temp_seq)) # add the string of the sequence
            random_list.append(str(temp_rand))
    return peak_list, random_list

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

def getRandomBackgroundSeq(chr, length):
    """
    Given a chromosome, get a random sequence of that length

    Parameters:
    -----------
        chr: chromosome number of random seq
        length: length of random seq that will be generated

    Returns:
    --------
        str: the sequence of nucleotides
    """
    length_of_chromosome = -1
    if chr in genome: # if chromosome is in the same format as the indexing
        length_of_chromosome = len(genome[chr])
    elif ("chr"+str(chr)) in genome: # if you need to add "chr" to the chromosome number to get it to match
        length_of_chromosome = len(str(genome["chr"+str(chr)]))
    elif str(chr).removeprefix("chr") in genome: # if you need to remove "chr" from the chromosome name to get it to match the index
        length_of_chromosome = len(str(genome[str(chr).removeprefix("chr")]))
    else:
        print("CHROMOSOME NOT FOUND!!!!! >:[ ") # TODO CHANGE THIS to a REAL ERROR Rather than this print
    start = random.randint(0, length_of_chromosome-length)
    return getSeqFromGenome(chr, start, start + length)

def generateListOfRandomSeqsFromPeaks(peaks_info):
    """
    Creates a list of random sequences from function getPeaksFromBed

    Parameters:
    -----------
    peaks_info: 
        list of known motifs

    Returns:
    --------

    """
    
    #peaks_info from getPeaksFromBed 
    random_seqs = []
    for ent in peaks_info:
        length = int(ent[2])-int(ent[1])
        chr = ent[0]
        random_seqs.append(getRandomBackgroundSeq(chr,length))
    return random_seqs

# Run Main
if __name__ == "__main__":
    main()