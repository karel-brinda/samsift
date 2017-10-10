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

BASIC_INIT="import random;"


def info(msg):
    dt = datetime.datetime.now()
    fdt = dt.strftime("%Y-%m-%d %H:%M:%S")
    print("[samsift]", fdt, msg, file=sys.stderr)


class SamSift:

    def __init__(self, sam_header, out_sam_fn, filter, code, dexpr, dtrig, mode, initialization):
        self.sam_header=sam_header
        self.output_sam_fn
        self.filter=filter
        self.code=code
        self.dexpr=dexpr
        self.dtrig=dtrig
        self.mode=mode
        self.initialization=initialization

        self.nt=0
        self.np=0
        self.nf=0
        self.ne=0

        info("Starting.")
        self.init_vardict={}
        exec(BASIC_INIT + initialization, self.init_vardict)
        exec(self.initialization, self.init_vardict)

    def _init_alignment(self, alignment):
        self.alignment=alignment

        self.err=False
        self.nt+=1

        #print(alignment.qual)
        self.vardict=self.init_vardict
        self.vardict.update({
                'a': alignment,
                'QNAME': alignment.query_name,
                'FLAG': alignment.flag,
                # RNAME will be set later,
                'POS': alignment.reference_start+1,
                'MAPQ': alignment.mapping_quality,
                'CIGAR': alignment.cigarstring,
                # RNEXT will be set later,
                'PNEXT': alignment.next_reference_start+1,
                'TLEN': alignment.template_length,
                'SEQ': alignment.query_sequence,

                # integer id's
                'RNAMEi': alignment.reference_id,
                'RNEXTi': alignment.next_reference_id,
                })

        # the specific implementation depends on the specific version of PySam, we want the same behaviour
        if isinstance(alignment.qual, str):
            self.vardict['QUAL']=alignment.qual
            self.vardict['QUALa']=[ord(x) for x in alignment.qual]
            self.vardict['QUALs']=alignment.qqual
            self.vardict['QUALsa']=[ord(x) for x in alignment.qqual]
        else:
            self.vardict['QUAL']=pysam.qualities_to_qualitystring(alignment.qual, offset=0)
            self.vardict['QUALa']=alignment.qual
            self.vardict['QUALs']=pysam.qualities_to_qualitystring(alignment.qqual, offset=0)
            self.vardict['QUALsa']=alignment.qqual

        if vardict['RNAMEi']==-1:
            self.vardict['RNAME']='*'
        else:
            self.vardict['RNAME']=self.in_sam.get_reference_name(vardict['RNAMEi'])

        if vardict['RNEXTi']==-1:
            self.vardict['RNEXT']='*'
        else:
            self.vardict['RNEXT']=self.in_sam.get_reference_name(vardict['RNEXTi'])

        # clean cache
        keys_to_delete=[]
        for k in self.vardict:
            if len(k)==2:
                keys_to_delete.append(k)
        for k in keys_to_delete:
            del self.vardict[k]

        self.vardict.update(alignment.get_tags())

    def _filter(self):
        try:
            passes=eval(filter, vardict)
        except:
            self.err=True
            if mode=="strict":
                raise
            elif mode=="nonstop-keep":
                self.passes=True
                self.ne+=1
            elif mode=="nonstop-remove":
                sefl.passes=False
                self.ne+=1
            else:
                raise NotImplementedError

    def _debug(self):
        if self.dexpr != "":
            trig=eval(str(dtrig), vardict)
            if trig:
                try:
                    dbg_res=eval(dexpr, vardict)
                except:
                    # todo: add a better message
                    dbg_res="evaluation_failed"
                print(alignment.query_name, bool(passes), dbg_res, file=sys.stderr, sep="\t")

    def _code(self):
        if self.code is not None:
            try:
                exec(self.code, self.vardict)
            except:
                if self.mode=="strict":
                    raise
                else:
                    if self.err==False:
                        self.ne+=1
                    self.err=True
                    info("Alignment '{}' - code error ('{}')".format(self.alignment.query_name, sys.exc_info()[0]))

            for k, v in vardict.items():
                if len(k)==2:
                    if isinstance(v, float):
                        # workaround (see "https://github.com/pysam-developers/pysam/issues/531")
                        self.alignment.set_tag(k, v, value_type="f")
                    else:
                        self.alignment.set_tag(k, v)


    def process_alignment(self, alignment):
        self._init_alignment(alignment)


        (err, passes)=self._filter()
        self._debug(alignment)
 
        if passes:
            self._code()
            out_sam.write(alignment)
            if not self.err:
                self.np+=1
        else:
            if not self.err:
                self.nf+=1


def sam_sift(in_sam_fn, out_sam_fn, filter, code, dexpr, dtrig, mode, initialization):
    if in_sam_fn=='-':
        info("Reading from standard input. " + ("Press Ctrl+D to finish reading or run '{} -h' for help.".format(PROGRAM) if in_sam_fn=="-" else""))

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
            for alignment in in_sam.fetch(until_eof=True):
                self.process_alignment(alignment)

    info("Finishing. {} alignments processed. {} alignments passed. {} alignments filtered out. {} alignments caused errors.".format(self.nt, self.np, self.nf, self.ne))



def parse_args():

    class CustomArgumentParser (argparse.ArgumentParser):
        def print_help(self):
            msg=self.format_help()
            msg=msg.replace("usage:", "Usage:  ")
            for x in 'PY_EXPR', 'PY_CODE':
                msg=msg.replace("[{x} [{x} ...]]\n            ".format(x=x), x)
                msg=msg.replace("[{x} [{x} ...]]".format(x=x), x)
            repl=re.compile(r'\]\s+\[')
            msg=repl.sub("] [",msg)
            msg=msg.replace("\n  -0","\n\nAdvanced options:\n  -0")
            msg=msg.replace(" [-h] [-v]","")
            msg=msg.replace("[-0","\n                    [-0")
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
    parser._optionals.title = 'Basic options'

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

    parser.add_argument('-m',
            choices=['strict', 'nonstop-keep', 'nonstop-remove'],
            metavar='STR',
            help='mode: strict (stop on first error)\n      nonstop-keep (keep alignments causing errors)\n      nonstop-remove (remove alignments causing errors) [strict]',
            dest='mode',
            default='strict',
        )

    parser.add_argument('-0',
            type=str,
            metavar='PY_CODE',
            help='initialization [None]',
            dest='init_l',
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


def main ():
    args=parse_args()
    sam_sift(
            in_sam_fn=args.in_sam_fn,
            out_sam_fn=args.out_sam_fn,
            filter=" ".join(args.filter_l),
            code=" ".join(args.code_l),
            dexpr=" ".join(args.dexpr_l),
            dtrig=" ".join(args.dtrig_l),
            initialization=" ".join(args.init_l),
            mode=args.mode,
        )


if __name__ == "__main__":
    main()
