#!/usb/bin/env python

from __future__ import print_function

import argparse
import getopt
import re
import sys

from port import open_or_die

usage = '''
Usage:

{} _seq.txt

This script parses _seq.txt to fasta format. Options:
-a    format is qseq.txt

Example of use:
{} mouse_seq.txt > mouse.fa'''.format(sys.argv[0], sys.argv[0])

running = 0


def parse_file_solexa(options, file_solexa):
    global running

    FILE = open_or_die(file_solexa, 'rb',
                       'can not open file: {}'.format(file_solexa))

    while True:
        l = FILE.readline()
        if not l:
            break

        if re.search(r'^\s*$', l):
            continue

        fields = re.split(r'\s+', l)
        seq = ''

        if options.get('-a') == '':
            seq = fields[8]
        else:
            seq = fields[4]

        seq = seq.strip()
        seq = re.sub(r'\.', 'N', seq)
        print('>{}_{}'.format(seq, running))
        print(seq)

        running += 1

    FILE.close()


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(-1)

    parser = argparse.ArgumentParser(usage)
    parser.add_argument('file_solexa', hep='file_solexa')

    args = parser.parse_args(sys.argv[1:2])
    file_solexa = args.file_solexa

    opts, argss = getopt.getopt(sys.argv[2:], 'a')
    options = dict(opts)

    parse_file_solexa(file_solexa)

    sys.exit(0)
