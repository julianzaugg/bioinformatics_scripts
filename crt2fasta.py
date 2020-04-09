#!/usr/bin/env python3

"""
Convert crisper recognition tool (CRT) output to FASTA. Specifically spacer sequences.
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-209

Example input:
ORGANISM:  NC_009656.1 Pseudomonas aeruginosa PA7, complete genome
Bases: 6588339


CRISPR 1   Range: 1204501 - 1204829
POSITION        REPEAT                          SPACER
--------        -----------------------------   -------------------------------
1204501         GTTCACTGCCGTATAGGCAGCTAAGAAAA   CGACGACGGCCAGCAACAGCGTCACCCTCGG [ 29, 31 ]
1204561         GTTCACTGCCGTATAGGCAGCTAAGAAAA   GCGGCAACATCACCACCAGCGCGGACATCAG [ 29, 31 ]
1204621         GTTCACTGCCGTATAGGCAGCTAAGAAAA   TGCGTTCAGAAGTCTTCATGCCGCAGTCCTC [ 29, 31 ]
1204681         GTTCACTGCCGTATAGGCAGCTAAGAAAA   CGGTGTCGCCGCAGCTTTACGCCCTGGCGTT [ 29, 31 ]
1204741         GTTCACTGCCGTATAGGCAGCTAAGAAAA   TCGTGATCAATGCCCGCCATGATTTCGCCGG [ 29, 31 ]
1204801         GTTCACTGCCGTATAGGCAGCTAAGAAAT
--------        -----------------------------   -------------------------------
Repeats: 6      Average Length: 29              Average Length: 31



CRISPR 2   Range: 2886935 - 2887202
POSITION        REPEAT                          SPACER
--------        ----------------------------    --------------------------------
2886935         GTTCACTGCCGTATAGGCAGCTAAGAAA    TAACTGTTTCGACAGTTCGCAATGCCAATGAA        [ 28, 32 ]
2886995         GTTCACTGCCGTATAGGCAGCTAAGAAA    GTTGAAGCGCACCAGGCCCGACTCCTCACCGA        [ 28, 32 ]
2887055         GTTCACTGCCGTATAGGCAGCTAAGAAA    TCTCCGCCGTTACCAACGCCATGAACCACGCG        [ 28, 32 ]
2887115         GTTCACTGCCGTATAGGCAGCTAAGAAA    ATGCTCCTGCGCAACGAGTTGATCGACGTGCT        [ 28, 32 ]
2887175         GTTCACTGCCGTATAGGCAGCTAAGAAA
--------        ----------------------------    --------------------------------
Repeats: 5      Average Length: 28              Average Length: 32
"""

__author__ = 'julianzaugg'

import sys
import contextlib
import argparse
import re

from Bio import SeqIO, Seq


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
        data=fh.read()

    spacer_sequence_records = []
    seq_set = set()
    for crispr_match_block in re.findall("CRISPR.*?Repeats",data,re.DOTALL):
        crispr_id = re.match("CRISPR\s\d*",crispr_match_block)[0].replace(" ", "_")
        crispr_block_values = [line.replace("\t\t", "\t").split("\t") for line in crispr_match_block.split("\n")[4:-2]]
        for entry_idx, entry in enumerate(crispr_block_values):
            position = entry[0]
            repeat_sequence = Seq.Seq(entry[1])
            repeat_name = "{}_REPEAT_{}".format(crispr_id, entry_idx + 1)
            spacer_sequence = Seq.Seq(entry[2])
            spacer_name = "{}_SPACER_{}".format(crispr_id,entry_idx+1)
            spacer_description = "position: {}".format(position)

            spacer_seq_record = SeqIO.SeqRecord(name=spacer_name,
                                                id=spacer_name,
                                                seq=spacer_sequence,
                                                description=spacer_description)
            if spacer_sequence != "":
                if input_args.unique and str(spacer_seq_record.seq) not in seq_set:
                    spacer_sequence_records.append(spacer_seq_record)
                    seq_set.add(str(spacer_seq_record.seq))
                elif not input_args.unique:
                    spacer_sequence_records.append(spacer_seq_record)

    with _smart_open(input_args.output, 'w') as fh:
        SeqIO.write(spacer_sequence_records,fh, format = "fasta-2line")


if  __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract spacer sequences from crisper recognition tool output.')
    parser.add_argument('-i', '--input',
                        help='Input should be output (file or STDOUT) from crisper recognition tool. Default STDIN',
                        default="-")
    parser.add_argument('-o', '--output',
                        help='Output filename. Default STDOUT',
                        default="-")
    parser.add_argument('-u', '--unique',
                        help='Ensure unique sequences',
                        action='store_true')

    args = parser.parse_args()
    _process_args(parser, args)