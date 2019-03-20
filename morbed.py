#!/usr/bin/env python
from __future__ import print_function

import re
import sys

from port import die, esplit, open_or_die

if __name__ == '__main__':
    if len(sys.argv) < 2:
        die('No csv file given for bed conversion\n')

    known, novel, _not, line, thres, score, line, strand, label, end = (
        None, None, None, None, None, None, None, None, None, None)

    IN = open_or_die(sys.argv[1], 'r', 'cannot open {}\n'.format(sys.argv[1]))
    while True:
        line = IN.readline()
        if not line:
            break

        if re.search(r'novel miRNAs predicted by moRNA Finder', line):
            novel = 1
            known = 0
            _not = 0
        elif re.search(r'mature miRBase miRNAs detected', line):
            novel = 0
            known = 1
            _not = 0
        else:
            l = esplit(line)
            if len(l) == 0:
                continue

            m = re.search(r'(\S+):(\d+)\.\.(\d+):(\S)', l[-1])
            if m:
                m = m.groups()
                end = m[2]
                if m[3] == '+':
                    strand = '255,0,0'

                if m[3] == '-':
                    strand = '0,0,255'

                if known:
                    label = 'known'

                if novel:
                    label = 'novel'

                print('{}\t{}\t{}\t{}:{}\t{}\t{}\t{}\t{}\t{}'.format(
                    m[0],
                    m[1],
                    end,
                    label,
                    l[0],
                    l[1],
                    m[3],
                    m[1],
                    end,
                    strand
                ))
