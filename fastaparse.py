#!/usr/bin/env python

from __future__ import print_function

import argparse
import getopt
import re
import sys

from port import open_or_die, print_stderr

usage = '''{} fastafile

This script parses a fastafile.

-a int    only output entries where the sequence is minimum int nts long
-b        remove all entries that have a sequence that contains letters
          other than a,c,g,t,u,n,A,C,G,T,U,N.
-s        output progress
'''.format(sys.argv[0])

running = 0


def resolve(options, _id, seq):
    global running

    running += 1

    if options.get('-s') == '':
        print_stderr('{}\r'.format(running))

    lng = len(seq)

    if options.get('-a') and lng < int(options.get('-a')):
        print_stderr('>{}\n{}\n'.format(_id, seq))
        return

    if options.get('-b') == '' and not re.match(r'^(a|c|g|t|u|n)+$', seq, re.IGNORECASE):
        print_stderr('>{}\n{}\n'.format(_id, seq))
        return

    print('>{}'.format(_id))
    print(seq)


def parse_fasta(options, file_fasta):
    _id = ''
    seq = ''

    FASTA = open_or_die(
        file_fasta, 'rb', 'can not open file {}\n'.format(file_fasta))

    while True:
        l = FASTA.readline()
        if not l:
            break

        l = l.strip()

        m = re.match(r'^>(\S+)', l)
        if m:
            _id = m.groups()[0]
            seq = ''

            while True:
                ll = FASTA.readline()
                if not ll:
                    break

                ll = ll.strip()

                mm = re.match(r'^>(\S+)', ll)
                if mm:
                    resolve(options, _id, seq)
                    _id = mm.groups()[0]
                    seq = ''
                    continue
                seq += ll

    resolve(options, _id, seq)
    FASTA.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(-1)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_fasta', help='fastafile')
    args = parser.parse_args(sys.argv[1:2])
    file_fasta = args.file_fasta

    opt, args = getopt.getopt(sys.argv[2:], 'a:bs')
    options = dict(opt)

    parse_fasta(options, file_fasta)
