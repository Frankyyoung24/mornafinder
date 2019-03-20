#!/usr/bin/env python

from __future__ import print_function

import argparse
import re
import sys

from port import esplit, open_or_die

usage = '''
Usage: {} fastq file > outputfile
'''.format(sys.argv[0])

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(-1)

    c = 0
    line = []
    _id = ''
    seq = ''
    processed = 0

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('filename', help='file')

    args = parser.parse_args(sys.argv[1:2])
    filename = args.filename

    IN = open_or_die(filename, 'rb', 'File {} not found.'.format(filename))

    while True:
        l = IN.readline().strip()
        if not l:
            break

        c += 1
        if c == 1:
            processed += 1
            line = esplit(l)
            _id = line[0]
            _id = re.sub(r'\@', '', _id, count=1)
        elif c == 2:
            seq = l
        elif c == 4:
            c = 0
            print('>{}_{}\n{}'.format(seq, processed, seq))
            _id = ''
            line = []
            seq = ''

    IN.close()
    sys.exit(0)
