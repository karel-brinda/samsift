#! /usr/bin/env python3

"""SAMsift - advanced filtering and tagging of SAM/BAM alignments using Python expressions

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import datetime
import os
import pysam
import re
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


def sam_sift(in_sam_fn, out_sam_fn, filter, code, dexpr, dtrig, mode):
	info("Starting.")
	if in_sam_fn=='-':
		info("Reading from standard input. Press Ctrl+D to finish reading or run '{} -h' for help.".format(PROGRAM))

	with pysam.AlignmentFile(in_sam_fn, "rb") as in_sam: #check_sq=False)
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

		if out_sam_fn[-4:]==".bam":
			out_mode="wb"
		else:
			out_mode="w"

		with pysam.AlignmentFile(out_sam_fn, out_mode, header=header) as out_sam:

			nt,np,nf,ne=0,0,0,0
			for alignment in in_sam.fetch(until_eof=True):
				err=False
				nt+=1
				#print(alignment.qual)
				vardict={
						'a': alignment,
						'QNAME': alignment.query_name,
						'FLAG': alignment.flag,
						'RNAME': in_sam.get_reference_name(alignment.reference_id),
						'POS': alignment.reference_start+1,
						'MAPQ': alignment.mapping_quality,
						'CIGAR': alignment.cigarstring,
						'PNEXT': alignment.next_reference_start+1,
						'TLEN': alignment.template_length,
						'SEQ': alignment.query_sequence,
						#
						'RNAMEi': alignment.reference_id,
						'RNEXTi': alignment.next_reference_id,
						}
				if isinstance(alignment.qual, str):
					vardict['QUAL']=alignment.qual
					vardict['QUALa']=[ord(x) for x in alignment.qual]
				else:
					vardict['QUAL']=pysam.qualities_to_qualitystring(alignment.qual, offset=0)
					vardict['QUALa']=alignment.qual

				if vardict['RNEXTi']==-1:
					vardict['RNEXT']='*'
				else:
					vardict['RNEXT']=in_sam.get_reference_name(vardict['RNEXTi']),
				vardict.update(alignment.get_tags())
				try:
					passes=eval(filter, vardict)
				except:
					err=True
					if mode=="strict":
						raise
					elif mode=="nonstop-keep":
						passes=True
						ne+=1
					elif mode=="nonstop-remove":
						passes=False
						ne+=1
					else:
						raise NotImplementedError
				if dexpr != "":
					trig=eval(str(dtrig), vardict)
					if trig:
						dbg_res=eval(dexpr, vardict)
						print(alignment.query_name, bool(passes), dbg_res, file=sys.stderr, sep="\t")
				if passes:
					if code is not None:
						try:
							exec(code, vardict)
						except:
							if mode=="strict":
								raise
							else:
								if err==False:
									ne+=1
								err=True
								info("Alignment '{}' - code error ('{}')".format(alignment.query_name, sys.exc_info()[0]))


						for k, v in vardict.items():
							if len(k)==2:
								if isinstance(v, float):
									# workaround (see "https://github.com/pysam-developers/pysam/issues/531")
									alignment.set_tag(k, v, value_type="f")
								else:
									alignment.set_tag(k, v)
					out_sam.write(alignment)
					if not err:
						np+=1
				else:
					if not err:
						nf+=1
	info("Finishing. {} alignments processed. {} alignments passed. {} alignments filtered out. {} alignments caused errors.".format(nt, np, nf, ne))

def main():

	class CustomArgumentParser (argparse.ArgumentParser):
		def print_help(self):
			msg=self.format_help()
			msg=msg.replace("usage:", "Usage:  ")
			for x in 'PY_EXPR', 'PY_CODE':
				msg=msg.replace("[{x} [{x} ...]]\n            ".format(x=x), x)
				msg=msg.replace("[{x} [{x} ...]]".format(x=x), x)
			repl=re.compile(r'\]\s+\[')
			msg=repl.sub("] [",msg)
			print(msg)

		def format_help(self):
			formatter = self._get_formatter()
			formatter.add_text(" \n"+self.description)
			formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
			formatter.add_text(self.epilog)

			# positionals, optionals and user-defined groups
			for action_group in self._action_groups:
				formatter.start_section(action_group.title)
				formatter.add_text(action_group.description)
				formatter.add_arguments(action_group._group_actions)
				formatter.end_section()

			return formatter.format_help()

	parser = CustomArgumentParser (
			formatter_class=argparse.RawTextHelpFormatter,
			description=
			"Program: {} ({})\n".format(PROGRAM, DESC)+
			"Version: {}\n".format(VERSION, pysam.__version__) +
			"Author:  Karel Brinda <kbrinda@hsph.harvard.edu>",
			)
	parser._optionals.title = 'Options'

	parser.add_argument('-v', '--version',
			action='version',
			version='{} {} (using pysam {})'.format(PROGRAM, VERSION, pysam.__version__),
			)

	parser.add_argument('-i',
			type=str,
			metavar='FILE',
			help="input SAM/BAM file [-]",
			dest='in_sam_fn',
			default='-',
			required=False,
			)

	parser.add_argument('-o',
			type=str,
			metavar='FILE',
			help="output SAM/BAM file [-]",
			dest='out_sam_fn',
			default='-',
			required=False,
			)

	parser.add_argument('-f',
			type=str,
			metavar='PY_EXPR',
			help='filter [True]',
			dest='filter_l',
			nargs='*',
			default=['True'],
			)

	parser.add_argument('-c',
			type=str,
			metavar='PY_CODE',
			help="code to be executed (e.g., assigning new tags) [None]",
			dest='code_l',
			nargs='*',
			default=['None'],
			)

	parser.add_argument('-d',
			type=str,
			metavar='PY_EXPR',
			help='debugging expression to print [None]',
			dest='dexpr_l',
			nargs='*',
			default=[],
		)

	parser.add_argument('-t',
			type=str,
			metavar='PY_EXPR',
			help='debugging trigger [True]',
			dest='dtrig_l',
			nargs='*',
			default=["True"],
		)

	parser.add_argument('-m',
			choices=['strict', 'nonstop-keep', 'nonstop-remove'],
			metavar='STR',
			help='mode: strict (stop on first error)\n      nonstop-keep (keep alignments causing errors)\n      nonstop-remove (remove alignments causing errors) [strict]',
			dest='mode',
			default='strict',
		)

	args = parser.parse_args()

	sam_sift(
			in_sam_fn=args.in_sam_fn,
			out_sam_fn=args.out_sam_fn,
			filter=" ".join(args.filter_l),
			code=" ".join(args.code_l),
			dexpr=" ".join(args.dexpr_l),
			dtrig=" ".join(args.dtrig_l),
			mode=args.mode,
		)


if __name__ == "__main__":
	main()
