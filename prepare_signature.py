#!/usr/bin/env python

from __future__ import print_function

import argparse
import getopt
import os
import re
import shutil
import sys
import time

from port import die, esplit, open_or_die, print_stderr


usage = '''
{} file_reads file_precursors read_align_edit_distance

This script prepares the signature file for moRNA Finder. Options:

-a file  Fasta file with the sequences of known mature miRNAs for the species.
         These sequences will not influence the moRNA Finder scoring, but will
         subsequently make it easy to estimate sensitivity of the run.
-b       Output progress to screen

'''.format(sys.argv[0])


def cat_files(file_1, file_2, file_out):
    OUT = open_or_die(file_out, 'w+', 'cannot print to {}\n'.format(file_out))
    IN_1 = open_or_die(file_1, 'rb', 'cannot read from {}\n'.format(file_1))
    while True:
        line = IN_1.readline()
        if not line:
            break
        OUT.write(line)
    IN_1.close()

    IN_2 = open_or_die(file_2, 'rb', 'cannot read from {}\n'.format(file_2))
    while True:
        line = IN_2.readline()
        if not line:
            break
        OUT.write(line)

    IN_2.close()
    OUT.close()


def presort(_file):
    IK = open_or_die(_file, 'rb', 'no arf file given\n')
    IKT = open_or_die('{}.tmp'.format(_file), 'w+',
                      'tmp file could not be created\n')
    index = {}
    count = 0
    l = []

    while True:
        line = IK.readline()
        if not line:
            break
        l = esplit(line)
        if l[5] not in index.keys():
            count += 1
            index[l[5]] = count

        IKT.write('{}\t{}'.format(index[l[5]], line))

    IK.close()
    IKT.close()


if __name__ == '__main__':
    if len(sys.argv) < 4:
        die(usage)

    parser = argparse.ArgumentParser(usage)
    parser.add_argument('file_reads', help='file reads')
    parser.add_argument('file_precursors', help='file precursors')
    parser.add_argument('read_align_edit_distance',
                        help='read align edit distance')
    args = parser.parse_args(sys.argv[1:4])

    file_reads = args.file_reads
    file_precursors = args.file_precursors
    read_align_edit_distance = args.read_align_edit_distance

    opts, argss = getopt.getopt(sys.argv[4:], "a:bo:")
    options = dict(opts)

    ltime = long(time.time())
    _dir = 'dir_prepare_signature{}'.format(ltime)

    if '-o' not in options.keys() or options.get('-o') == '':
        die('no outfile specified with option -o\n')

    outfile = options.get('-o')

    if options.get('-b') == '':
        print_stderr('preparing signature file\n')

    os.mkdir(_dir)
    shutil.copy(file_precursors, _dir)

    if options.get('-b') == '':
        print_stderr('constructing index of precursors\n')
    os.system(
        'bowtie-build {} {}/precursors.ebwt > /dev/null'.format(file_precursors, _dir))

    if options.get('-b') == '':
        print_stderr('mapping reads to precursors\n')
    cmd = 'bowtie -f -v {} -a --best --strata --norc {}/precursors.ebwt {} {}/reads_vs_precursors.bwt 2> /dev/null\n'.format(
        read_align_edit_distance,
        _dir,
        file_reads,
        _dir
    )
    print_stderr(cmd)
    os.system(cmd)

    os.system('convert_bowtie_output.py {}/reads_vs_precursors.bwt > {}/reads_vs_precursors.arf'.format(
        _dir,
        _dir
    ))

    if options.get('-a'):
        file_mature = options.get('-a')

        if options.get('-b'):
            print_stderr('mapping reference mature miRNAs to precursors\n')
        os.system('bowtie -f -v 0 -a --best --strata --norc {}/precursors.ebwt {} {}/mature_vs_precursors.bwt 2> /dev/null'.format(
            _dir,
            options.get('-a'),
            _dir
        ))

        os.system(
            'convert_bowtie_output.py {}/mature_vs_precursors.bwt > {}/mature_vs_precursors.arf'.format(_dir, _dir))

        if options.get('-b') == '':
            print_stderr('sorting rows\n')

        cat_files("{}/reads_vs_precursors.arf".format(_dir),
                  "{}/mature_vs_precursors.arf".format(_dir),
                  "{}/signature_unsorted.arf".format(_dir))

        presort('{}/signature_unsorted.arf'.format(_dir))
        os.system(
            'sort -nk1 {}/signature_unsorted.arf.tmp > {}/signature_unsorted.arf.tmp2'.format(_dir, _dir))
        os.system(
            'cut -f2-14 {}/signature_unsorted.arf.tmp2 > {}'.format(_dir, outfile))
    else:
        if options.get('-b') == '':
            print_stderr('sorting rows\n')

        presort('{}/reads_vs_precursors.arf'.foramt(_dir))
        os.system(
            'sort -nk1 {}/reads_vs_precursors.arf.tmp > {}/reads_vs_precursors.arf.tmp2'.format(_dir, _dir))
        os.system(
            'cut -f2-14 {}/reads_vs_precursors.arf.tmp2 > {}'.format(_dir, outfile))

    if options.get('-b') == '':
        print_stderr('signature file prepared\n\n')
