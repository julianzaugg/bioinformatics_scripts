#!/usr/bin/env python3

""" Download fasta formatted sequence for a provided uniprot accession
"""

__author__ = 'julianzaugg'

import argparse
import contextlib
import sys

#import pandas as pd
import requests as r
from Bio import SeqIO
from io import StringIO


@contextlib.contextmanager
def _smart_open(filename, mode='Ur'):
    """Opens filename and returns a file handle. 
    If filename is None or '-', uses stdout.
    """
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

def _download_uniprot(parser, args):
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl + args.accession + ".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)

    sequence = SeqIO.parse(StringIO(cData), "fasta")
    with _smart_open(args.output, 'w') as fh:
        SeqIO.write(sequence, fh, format = "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-a', '--accession', help='Uniprot accession. Default STDIN', default="-")
    parser.add_argument('-o', '--output', help='Output filename. Default STDOUT', default="-")
    args = parser.parse_args()
    _download_uniprot(parser, args)