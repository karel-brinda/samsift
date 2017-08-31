#! /usr/bin/env python3

"""SAMsift - sift and enrich SAM/BAM alignments using Python expressions

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import datetime
import os
import pysam
import sys

sys.path.append(os.path.dirname(__file__))
import version

PROGRAM='samsift'
VERSION=version.VERSION
DESC='advanced filtering and tagging of SAM/BAM alignments using Python expressions'


def info(msg):
	dt = datetime.datetime.now()
	fdt = dt.strftime("%Y-%m-%d %H:%M:%S")
	print("[samsift]", fdt, msg, file=sys.stderr)


def sam_sift(in_sam_fn, out_sam_fn, sieve, code, dexpr, dtrig):
	info("Starting.")
	if in_sam_fn=='-':
		info("Reading from standard input. Press Ctrl+D to finish reading or run '{} -h' for help.".format(PROGRAM))

	in_sam=pysam.AlignmentFile(in_sam_fn, "rb") #check_sq=False)
	#print("@PG", "ID:{}".format(PROGRAM), "PN:{}".format(PROGRAM), "VN:{}".format(VERSION), "CL:{}".format(" ".join(sys.argv)), sep="\t")
	header=in_sam.header

	pg={
			"ID":PROGRAM,
			"PN":PROGRAM,
			"VN":VERSION,
			"CL":" ".join(map(lambda x:"'{}'".format(x),sys.argv)),
		}

	try:
		header['PG'].insert(0,pg)
	except KeyError:
		header['PG']=[pg]

	out_sam = pysam.AlignmentFile(out_sam_fn, "w", header=header)

	nt,np,nf=0,0,0
	for a in in_sam.fetch(until_eof=True):
		nt+=1
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
		if dexpr is not None:
			trig=eval(str(dtrig), vardict)
			if trig:
				res=eval(dexpr, vardict)
				print(a.query_name, bool(p), res, file=sys.stderr, sep="\t")
		if p:
			if code is not None:
				vardict_new=vardict.copy()
				exec(code, vardict)
				for k,v in vardict.items():
					if len(k)==2:
						if isinstance(v, float):
							# workaround (see "https://github.com/pysam-developers/pysam/issues/531")
							a.set_tag(k, v, value_type="f")
						else:
							a.set_tag(k, v)
			out_sam.write(a)
			np+=1
		else:
			nf+=1
	info("Finishing. {} alignments processed. {} alignments passed. {} alignments filtered out.".format(nt, np, nf))

def main():

	parser = argparse.ArgumentParser(description=
			"Program: {} ({})\n".format(PROGRAM, DESC)+
			"Version: {}\n".format(VERSION) +
			"Author:  Karel Brinda <kbrinda@hsph.harvard.edu>",
			formatter_class=argparse.RawDescriptionHelpFormatter
			)

	parser.add_argument('-v', '--version',
			action='version',
			version='{} {}'.format(PROGRAM, VERSION),
			)

	parser.add_argument('-i',
			type=str,
			metavar='file',
			dest='in_sam_fn',
			default='-',
			required=False,
			help="input SAM/BAM file [-]",
			)

	parser.add_argument('-o',
			type=str,
			metavar='file',
			dest='out_sam_fn',
			default='-',
			required=False,
			help="output SAM/BAM file [-]",
			)

	parser.add_argument('-f',
			type=str,
			metavar='py_expr',
			help='filter [True]',
			dest='sieve',
			required=False,
			default='True'
			)

	parser.add_argument('-c',
			type=str,
			metavar='py_code',
			dest='code',
			default=None,
			required=False,
			help="code to be executed (e.g., assigning new tags) [None]",
			)

	parser.add_argument('-d',
			type=str,
			metavar='py_expr',
			dest='dexpr',
			help='debugging expression to print [None]',
			default=None,
		)

	parser.add_argument('-t',
			type=str,
			metavar='py_expr',
			dest='dtrig',
			help='debugging trigger [True]',
			default="True",
		)

	args = parser.parse_args()

	sam_sift(
			in_sam_fn=args.in_sam_fn,
			out_sam_fn=args.out_sam_fn,
			sieve=args.sieve,
			code=args.code,
			dexpr=args.dexpr,
			dtrig=args.dtrig,
		)


if __name__ == "__main__":
	main()
