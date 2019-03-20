#!/usr/bin/env python

from __future__ import print_function

import re
import sys

from port import die, open_or_die


usage = '''
usage:

    {} file_fasta

    Removes whitespaces from id line in a fasta file

'''.format(sys.argv[0])


def parse_fasta(file_fasta):
    FASTA = open_or_die(
        file_fasta, 'rb', 'can not open {}\n'.format(file_fasta))
    while True:
        line = FASTA.readline()
        if not line:
            break

        m = re.findall(r'^(>\S+)', line)
        if m:
            print('{}\n'.format(m[0]))
        else:
            print(line.upper())

    FASTA.close()
    return 0


if __name__ == '__main__':
    if len(sys.argv) < 2:
        die(usage)

    file_fasta = sys.argv[1]

    parse_fasta(file_fasta)

    sys.exit(0)
