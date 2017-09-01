#! /usr/bin/env python3

"""SAMsift - sift and enrich SAM/BAM alignments using Python expressions

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


def sam_sift(in_sam_fn, out_sam_fn, filter, code, dexpr, dtrig):
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

		with pysam.AlignmentFile(out_sam_fn, "w", header=header) as out_sam:

			nt,np,nf=0,0,0
			for a in in_sam.fetch(until_eof=True):
				nt+=1
				vardict={
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
				vardict.update(a.get_tags())
				p=eval(filter, vardict)
				if dexpr != "":
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
			formatter_class=argparse.RawDescriptionHelpFormatter,
			description=
			"Program: {} ({})\n".format(PROGRAM, DESC)+
			"Version: {}\n".format(VERSION) +
			"Author:  Karel Brinda <kbrinda@hsph.harvard.edu>",
			)
	parser._optionals.title = 'Options'

	parser.add_argument('-v', '--version',
			action='version',
			version='{} {}'.format(PROGRAM, VERSION),
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

	args = parser.parse_args()

	sam_sift(
			in_sam_fn=args.in_sam_fn,
			out_sam_fn=args.out_sam_fn,
			filter=" ".join(args.filter_l),
			code=" ".join(args.code_l),
			dexpr=" ".join(args.dexpr_l),
			dtrig=" ".join(args.dtrig_l),
		)


if __name__ == "__main__":
	main()
