#! /usr/bin/env python3

"""SAMsift - sift alignments in SAM/BAM

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import os
import sys
import pysam

sys.path.append(os.path.dirname(__file__))
from version import *

def sam_sift(in_sam_fn, out_sam_fn, sieve, dexprs=[]):
	in_sam=pysam.AlignmentFile(in_sam_fn, "rb") #check_sq=False)
	out_sam = pysam.AlignmentFile(out_sam_fn, "w", template=in_sam)
	for a in in_sam.fetch(until_eof=True):
		vardict1={
			'a': a,
			'QNAME': a.query_name,
			'FLAG': a.flag,
			'RNAME': a.reference_id,
			'POS': a.reference_start+1,
			'MAPQ': a.mapping_quality,
			'CIGAR': a.cigarstring,
			'RNEXT': a.next_reference_id,
			'PNEXT': a.next_reference_start+1,
			'TLEN': a.template_length,
			'SEQ': a.query_sequence,
			'QUAL': a.query_qualities,
		}
		vardict2=dict(a.get_tags())
		vardict={**vardict1, **vardict2}
		p=eval(sieve, vardict)
		if len(dexprs)>0:
			res=[eval(dexpr, vardict) for dexpr in dexprs]
			print(a.query_name, bool(p), *res, file=sys.stderr, sep="\t")
		if p:
			out_sam.write(a)

def main():
	parser = argparse.ArgumentParser(description="Program: samsift (sift your SAM files)\nVersion: {}\nAuthor: Karel Brinda <kbrinda@hsph.harvard.edu>".format(VERSION), formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('sieve',
			type=str,
			metavar='expr',
			help='sieve (a Python expression)',
		)

	parser.add_argument('dexprs',
			type=str,
			metavar='debug_expr',
			help='debugging expression (a Python expression)',
			nargs="+",
		)

	parser.add_argument('-i',
			type=str,
			metavar='in.sam',
			dest='in_sam_fn',
			default='-',
			required=False,
			help="input SAM/BAM file ['-', i.e., stdin]",
		)

	parser.add_argument('-o',
			type=str,
			metavar='out.sam',
			dest='out_sam_fn',
			default='-',
			required=False,
			help="output SAM/BAM file ['-', i.e., stdout]",
		)

	args = parser.parse_args()

	sam_sift(args.in_sam_fn, args.out_sam_fn, args.sieve, args.dexprs)

if __name__ == "__main__":
	main()
