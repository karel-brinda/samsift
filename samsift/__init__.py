#! /usr/bin/env python3

"""SAMsift - sift alignments in SAM/BAM

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""


import argparse
import os
import sys
import pysam

# todo: fix pipes
# todo: specify scope
# todo: put sam spec fields into variables
# todo: put tags into variables

def sam_sift(in_sam_fn, out_sam_fn, sieve):
	in_sam=pysam.AlignmentFile(in_sam_fn, "rb") #check_sq=False)
	out_sam = pysam.AlignmentFile(out_sam_fn, "w", template=in_sam)
	for read in in_sam.fetch(until_eof=True):
		if eval(sieve):
			out_sam.write(read)

def main():
	parser = argparse.ArgumentParser(description="Program: samsift (sift your SAM files)\nVersion: 0.0.1\nAuthor: Karel Brinda (kbrinda@hsph.harvard.edu)")

	parser.add_argument('sieve',
			type=str,
			metavar='expr',
			help='sieve (a Python expression)',
		)

	parser.add_argument('-i',
			type=str,
			metavar='str',
			dest='in_sam_fn',
			default='-',
			required=False,
			help='input SAM/BAM file [-]',
		)

	parser.add_argument('-o',
			type=str,
			metavar='str',
			dest='out_sam_fn',
			default='-',
			required=False,
			help='output SAM/BAM file [-]',
		)

	args = parser.parse_args()

	sam_sift(args.in_sam_fn, args.out_sam_fn, args.sieve)

if __name__ == "__main__":
	main()
