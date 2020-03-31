#!/usr/bin/env python3

__author__ = 'julianzaugg'


import sys
from pathlib import Path
import contextlib
import argparse

from Bio import SeqIO

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
    if input_args.output_dir == None:
        output_dir = Path(input_args.input).parent
    else:
        output_dir = input_args.output_dir
    with _smart_open(input_args.input, 'r') as fh:
        sequences = SeqIO.parse(fh, "fasta")
        for s in sequences:
            outname="{}.fasta".format(s.name)
            with open(Path(output_dir) / outname, 'w') as fh2:
                SeqIO.write(s, fh2, format="fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split a FASTA file into individual files, one per entry.'
                                                 'Assumes unique names.')
    parser.add_argument('-i', '--input',
                        help='Input fasta file. Default STDIN',
                        default="-")
    parser.add_argument('-o', '--output_dir',
                        help='Output directory. Default is the directory of input.')

    args = parser.parse_args()
    _process_args(parser, args)