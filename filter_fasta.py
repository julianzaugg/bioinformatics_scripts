#!/usr/bin/env python3

__author__ = 'julianzaugg'

import sys
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
    with _smart_open(input_args.input, 'r') as fh:
        sequences = list(SeqIO.parse(fh, "fasta"))
    N_seqs = len(sequences)
    filtered_sequences = sequences
    if input_args.min_length != None:
        print("Filtering by minimum length of: " + str(input_args.min_length))
        filtered_sequences = [s for s in filtered_sequences if len(s) >= input_args.min_length]
    if input_args.max_length != None:
        print("Filtering by maximum length of: " + str(input_args.max_length))
        filtered_sequences = [s for s in filtered_sequences if len(s) <= input_args.max_length]

    print("Filtering complete. Removed " + str(N_seqs - len(filtered_sequences)) + " sequences")
    with _smart_open(input_args.output, 'w') as fh:
        SeqIO.write(filtered_sequences,fh, format = "fasta-2line")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter FASTA file')
    parser.add_argument('-i', '--input', #type=argparse.FileType('r'),
                        help='Input FASTA file. Default STDIN',
                        default="-")
    parser.add_argument('-o', '--output', #type=argparse.FileType('w'),
                        help='Output filename. Default STDOUT',
                        default="-")

    parser.add_argument("--min_length",
                        help = "Sequences with length less than this value are removed",
                        type = int,
                        required=False)
    parser.add_argument("--max_length",
                        help = "Sequences with length more than this value are removed",
                        type=int,
                        required=False)

    args = parser.parse_args()
    _process_args(parser, args)
