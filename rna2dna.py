#!/usr/bin/env python

import re
import sys

from port import die, open_or_die, pprint


usage = '''
usage:
    {} file_fasta

    Parses RNA sequences to DNA by substituting 'U's with 'T'.

'''.format(sys.argv[0])


def parse_fasta(file_fasta):
    FASTA = open_or_die(
        file_fasta, 'rb', 'can not open {}\n'.format(file_fasta))
    while True:
        line = FASTA.readline()
        if not line:
            break

        m = re.match(r'^(>\S+)', line)
        if m:
            pprint('{}\n'.format(m.groups()[0]))
        else:
            pprint('{}'.format(re.sub('U', 'T', line).upper()))

    FASTA.close()
    return


if __name__ == '__main__':
    if len(sys.argv) < 2:
        die(usage)

    file_fasta = sys.argv[1]

    parse_fasta(file_fasta)

    sys.exit(0)
