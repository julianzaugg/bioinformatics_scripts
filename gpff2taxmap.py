#!/usr/bin/env python3

# Take gpff file and generate a tab-delimited tax map file for building BLAST database

import sys
import argparse
import contextlib
import re
from Bio import SeqIO # The Python Bio package needs to be installed

@contextlib.contextmanager
def _smart_open(filename, mode='Ur'):
    if filename == '-':
        if mode is None or mode == '' or 'r' in mode:
            fh = sys.stdin
        else:
            fh = sys.stdout
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        if filename != '-':
            fh.close()

def _process_args(input_parser, input_args):
    with _smart_open(input_args.input, 'r') as fh:
        records = SeqIO.parse(input_args.input, "genbank")

    for seq_record in records:
        for feature in seq_record.features:
            if 'source' in feature.type:
               taxid=''.join(feature.qualifiers['db_xref'])
               taxid=re.sub(r'.*taxon:','',taxid)
               with _smart_open(args.output, 'a') as fh:
                   print(seq_record.id+'\t'+taxid, file = fh)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate tax map file from gbff file.')
    parser.add_argument('-i', '--input', #type=argparse.FileType('r'),
                        help='GBFF filename. Default STDIN',
                        default="-")
    parser.add_argument('-o', '--output', #type=argparse.FileType('w'),
                        help='Output filename. Default STDOUT',
                        default="-")
    if len(sys.argv[1:])==0:
        parser.print_help()
        # parser.print_usage() # for just the usage line
        parser.exit()
    args = parser.parse_args()
    _process_args(parser, args)
