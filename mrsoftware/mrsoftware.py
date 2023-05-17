#!/usr/bin/env/ python

"""
bleh bleh
"""

import argparse
from . import myutils as myutils
from mrsoftware import __version__
import os
import pyfaidx
import sys

def main():
    parser = argparse.ArgumentParser(
		prog="mrsoftware",
		description="Command-line script to perform  motif recognition from BED file"
	)

	# Input files
        # input Bed File
    parser.add_argument("bed", help="BED file", type=str, required=True)
        # input meme file of known motifs  
    parser.add_argument("-m", "--meme", help="meme file corresponding to the transcription factor", type=str, metavar="FILE", required=True)
        # input the reference genome to use
    parser.add_argument("-g","--genome", help="Path to the reference genome that peaks will be pulled from", type=str, metavar="FILE", required=True)
    
	# Output file to write to
    parser.add_argument("-o", "--out", help="Write output to file. Default: stdout", metavar="FILE", type=str, required=False)

	# Other options
    parser.add_argument("--version", help="Print the version and quit", action="version", version = '{version}'.format(version=__version__))
    
    # Parse args
    args = parser.parse_args()

	# Set up output file
    if args.out is None:
        outf = sys.stdout
    else: outf = open(args.out, "w")

    # Load  genome FASTA
    if args.fasta_ref is not None:
        if not os.path.exists(args.fasta_ref):
            myutils.ERROR("{fasta} does not exist".format(fasta=args.fasta_ref))
        reffasta = pyfaidx.Fasta(args.fasta_ref)
    else:
        reffasta = None