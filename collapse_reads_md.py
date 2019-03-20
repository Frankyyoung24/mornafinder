#!/usr/bin/env python

from __future__ import print_function

import argparse
import getopt
import re
import sys

from port import (create_hash_key_chain, hash_sort_key, open_or_die2,
                  print_stderr, tr)


usage = '''
{} file_fasta prefix

Collapses reads in the fasta file to make each sequence entry unique. Each collapsed
entry will have an identifier that follows an underscore '_' separated format. Example:

>mmu_1189273_x10

The first field 'mmu' shows which sample the sequence is from. This prefix is given on
the command line, and must consist of exactly three alphabetic letters. The second field
'118273' is a running number that secures that each identifier is unique. The third
field '_x10' shows how many times the given sequence was present in the dataset.

-a    outputs progress

example use:
{} reads.fa mmu'''.format(sys.argv[0], sys.argv[0])


def find_cnt(_id):
    m = re.search(r'_x(\d+)', _id)
    if m:
        cnt = m.group()
        return cnt
    else:
        return 1


def test_prefix(prefix):
    # unless($prefix=~/^\w\w\w$/ and !($prefix=~/_/))
    if not (re.search(r'^\w\w\w$', prefix) and (not re.search(r'_', prefix))):
        print('prefix {} does not contain exactly three alphabet letters'.format(prefix))
        sys.exit(-1)


def parse_file_fasta_seqkey(file_fasta, hsh, options):
    if options.get('-a') == '':
        print_stderr('reading file into hash\n')

    _id = ''
    seq = ''
    running_1 = 0

    FASTA = open_or_die2(file_fasta, 'rb')

    while True:
        l = FASTA.readline().strip()
        if not l:
            break

        m = re.match(r'^>(\S+)', l)
        if m:
            _id = m.group()
            seq = ''

            while True:
                ll = FASTA.readline().strip()
                if not ll:
                    break

                mm = re.match(r'^>(\S+)', ll)
                if mm:
                    cnt = find_cnt(_id)
                    seq = tr(seq, '[acgtun.]', '[ACGTTNN]')
                    # ATTR: Performance issue below:
                    # create_hash_key_chain(hsh, 0, seq)
                    try:
                        hsh[seq] = (hsh[seq]) + cnt
                    except KeyError:
                        hsh[seq] = cnt

                    running_1 += 1

                    if options.get('-a') == '':
                        print_stderr('{}\r'.format(running_1))

                    _id = mm.group()
                    seq = ''
                    continue

                seq += ll

    cnt = find_cnt(_id)
    seq = tr(seq, '[acgtun.]', '[ACGTTNN]')
    create_hash_key_chain(hsh, 0, seq)
    hsh[seq] += cnt
    running_1 += 1

    if options.get('-a') == '':
        print_stderr('{}\r'.format(running_1))

    FASTA.close()


def print_hash_seqkey(hsh):
    if options.get('-a') == '':
        print_stderr('sorting hash\n')

    running_2 = 0
    if options.get('-a') == '':
        print_stderr('printing hash\n')

    keys = hash_sort_key(hsh, lambda x: (x[1] * -1, x[0]))
    for key in keys:
        cnt = hsh[key]
        # print ">$prefix\_$running_2\_x$cnt\n$key\n";
        print('>{}_{}_x{}\n{}'.format(prefix, running_2, cnt, key))
        running_2 += cnt

        if options.get('-a') == '':
            print_stderr('{}\r'.format(running_2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_fasta', help=usage)
    parser.add_argument('prefix', help=usage)

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(-1)

    args = parser.parse_args(sys.argv[1:3])
    file_fasta = args.file_fasta
    prefix = args.prefix

    opts, argss = getopt.getopt(sys.argv[3:], 'a')
    options = dict(opts)

    test_prefix(prefix)
    hsh = {}

    parse_file_fasta_seqkey(file_fasta, hsh, options)

    print_hash_seqkey(hsh)
