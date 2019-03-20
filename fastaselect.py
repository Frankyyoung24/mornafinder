#!/usr/bin/env python

import argparse
import getopt
import re
import sys

from port import open_or_die


usage = '''{} file_fasta idfile (-a)

This script only prints out the fasta entries that has an id
that is present in the idfile. Using the option -a only prints
out entries that has an id that is not present in the id file.
'''.format(sys.argv[0])


def parse_file_ids(file_ids, _hash):

    FILE = open_or_die(file_ids, 'rb', 'can not open {}'.format(file_ids))
    while True:
        line = FILE.readline()
        if not line:
            break

        m = re.match(r'^(\S+)', line)
        if m:
            _id = m.groups()[0]
            _hash[_id] = 1

    FILE.close()


def parse_fasta(file_fasta, _hash):
    FASTA = open_or_die(
        file_fasta, 'rb', 'can not open file'.format(file_fasta))
    while True:
        l = FASTA.readline()
        if not l:
            break

        l = l.strip()

        m = re.match(r'^>(\S+)(.*)', l)
        if m:
            m = m.groups()
            _id = m[0]
            desc = m[1]
            sequence = ''

            while True:
                ll = FASTA.readline()
                if not ll:
                    break

                ll = ll.strip()

                mm = re.match(r'^>(\S+)(.*)', ll)
                if mm:
                    mm = mm.groups()
                    defined = False
                    try:
                        _hash[_id]
                        defined = True
                    except KeyError:
                        pass

                    if (defined and options.get('-a') is None) or (not defined
                                                                   and options.get('-a') == ''):
                        print('>{}{}\n{}'.format(_id, desc, sequence))

                    _id = mm[0]
                    desc = mm[1]
                    sequence = ''
                    continue

                sequence += ll

    defined = False
    try:
        _hash[_id]
        defined = True
    except KeyError:
        pass
    if (defined and options.get('-a') is None) or (not defined and options.get('-a') == ''):
        print('>{}{}\n{}'.format(_id, desc, sequence))

    FASTA.close()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(usage)
        sys.exit(-1)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_fasta', help='fasta file')
    parser.add_argument('idfile', help='idfile')
    args = parser.parse_args(sys.argv[1:3])

    file_fasta = args.file_fasta
    file_ids = args.idfile

    opts, argss = getopt.getopt(sys.argv[3:], 'a')
    options = dict(opts)

    _hash = {}

    parse_file_ids(file_ids, _hash)

    parse_fasta(file_fasta, _hash)

    sys.exit(0)
