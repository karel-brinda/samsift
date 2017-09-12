#! /usr/bin/env python3

"""Normalize SAM in order to faciliate SAM files comparison.

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

PROGRAM='samsift-norm-sam'
VERSION=version.VERSION
DESC='normalize SAM in order to faciliate SAM files comparison'


def sam_norm(in_sam_fn, out_sam_fn):
    if in_sam_fn[-4:]==".bam":
        in_mode="rb"
    else:
        in_mode="r"

    if out_sam_fn[-4:]==".bam":
        out_mode="wb"
    else:
        out_mode="w"

    with pysam.AlignmentFile(in_sam_fn, in_mode) as in_sam: #check_sq=False)
        header=in_sam.header


        buffer_alignments=[]
        qname=None

        with pysam.AlignmentFile(out_sam_fn, out_mode, header=header) as out_sam:

            for alignment in in_sam.fetch(until_eof=True):

                #print("alignment", alignment.qname)

                if alignment.qname!=qname and len(buffer_alignments):
                    buffer_alignments.sort(key=lambda x: x.reference_id)
                    for x in buffer_alignments:
                        out_sam.write(alignment)
                    buffer_alignments=[]

                qname=alignment.qname

                tags=alignment.get_tags()
                tags.sort()
                alignment.set_tags(tags)

                buffer_alignments.append(alignment)

            buffer_alignments.sort(key=lambda x: x.reference_id)
            for alignment in buffer_alignments:
                out_sam.write(alignment)
                buffer_alignments=[]



def main():

    class CustomArgumentParser (argparse.ArgumentParser):
        def print_help(self):
            msg=self.format_help()
            msg=msg.replace("usage:", "Usage:  ")
            repl=re.compile(r'\]\s+\[')
            msg=repl.sub("] [",msg)
            msg=msg.replace(" [-h] [-v]","")
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

    args = parser.parse_args()

    sam_norm(
            in_sam_fn=args.in_sam_fn,
            out_sam_fn=args.out_sam_fn,
        )


if __name__ == "__main__":
    main()
