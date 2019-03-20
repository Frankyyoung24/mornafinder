#!/usr/bin/env python

from __future__ import print_function

import re
import sys

from port import ssplit, str_reverse, tr

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('usage: {} reads_mapped.bwt'.format(sys.argv[0]))
        sys.exit(-1)

    line = []
    gseq = []
    changes = []
    pos = None
    ont = None
    mm = 0
    edit = ''
    reverse = 0

    try:
        IN = open(sys.argv[1], 'rb')
    except IOError:
        print('cannot open file {}'.format(sys.argv[1]))
        sys.exit(-1)

    while True:
        l = IN.readline()
        if not l:
            break

        line = re.split(r'\t', l)
        if line[1] == '-':
            line[4] = str_reverse(line[4])
            line[4] = tr(line[4], 'ACGTN', 'TGCAN')

        gseq = ssplit(line[4].lower())
        edit = ssplit('m' * len(line[4]))
        mm = 0

        if line[7]:
            changes = re.split(r',', line[7])
            for change in changes:
                match = re.search(r'(\d+):(\w+)\>\w+', change)
                if match:
                    match = match.groups()
                    mm += 1
                    gseq[int(match[0])] = match[1].lower()
                    edit[int(match[0])] = 'M'

        _id = re.split(r'\s', line[0])
        db = re.split(r'\s', line[2])

        print('{}\t{}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            _id[0],
            len(line[4]),
            len(line[4]),
            line[4].lower(),
            db[0],
            len(line[4]),
            int(line[3]) + 1,
            int(line[3]) + len(line[4]),
            ''.join(gseq),
            line[1],
            mm,
            ''.join(edit)
        ))

    IN.close()
