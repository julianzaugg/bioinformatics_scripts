#!/usr/bin/env python3

__author__ = 'julianzaugg'

import sys
import argparse

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from collections import Counter

def _process_args(input_parser, input_args):

    align = AlignIO.read(args.input, "fasta")

    n_seqs = len(align)
    alignment_length = align.get_alignment_length()
    columns_to_keep = []

    for position in range(alignment_length):
        if Counter(align[:, position])["-"]/n_seqs >= args.col_gap_fraction:
            continue
        columns_to_keep.append(position)

    new_seqs = []
    for seq in align:
        new_seq = SeqRecord(Seq(''.join([seq[col] for col in columns_to_keep])), id=seq.name,name="",description="")
        fraction_gap = Counter(new_seq)["-"]/len(new_seq)
        if fraction_gap >= args.seq_gap_fraction:
            continue
        new_seqs.append(new_seq)

    new_alignment = MultipleSeqAlignment(records =new_seqs)

    with open(args.output, 'w') as fh:
        SeqIO.write(new_alignment,fh, format = "fasta")

def _restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Remove gappy columns and sequences from an alignment.')
    parser.add_argument('-i', '--input', help='Input fasta file containing aligned sequences.', required=True)
    parser.add_argument('-o', '--output', help='Output filename. Default STDOUT', required=True, default="-")

    parser.add_argument("-c", "--col_gap_fraction",
                        help = "Columns with fraction of gaps higher or equal to this are removed. Default = 1.0",
                        type=_restricted_float, default=1.0,required=False)
    parser.add_argument("-s", "--seq_gap_fraction",
                        help = "Sequences with fraction of gaps higher or equal to this are removed. Default = 1.0",
                        type=_restricted_float, default=1.0,required=False)

    args = parser.parse_args()
    _process_args(parser, args)
