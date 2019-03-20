#!/usr/bin/env python

import re
import sys

from port import Nicenumber, die, open_or_die2

counter = 0


def read_handler(handle):
    global counter
    while True:
        rin = handle.readline().strip()
        if not rin:
            break

        rin = rin.strip()
        counter += 1
        if re.match(r'^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+(\d+)\s*([mDIM]*)$', rin):
            pass
        else:
            die('\nWrong format in line {}: The row\n{}\ndoes not correspond to the format\nread_id_wo_whitespaces\tlength\tstart\tend\tread_sequence         \tgenomicID_wo_whitspaces\tlength\tstart\tend\tgenomic_sequence       \tstrand\t#mismatches\teditstring\ne.g. read_22_x10000 \t22    \t1    \t22 \tagtcgtgactgactgactgacg\tchromosomeIII_x12312312\t22    \t1001 \t1022\tagtcgtgactgactgactgacg\t+-    \t0          \tmmmmmmmmmmmmmmmmmmmmm\nPlease make sure that all lines have the above described format.\n'.format(
                Nicenumber(counter),
                rin
            ))


if __name__ == '__main__':

    if len(sys.argv) == 1:
        read_handler(sys.stdin)
    else:
        for f in sys.argv[1:]:
            IN = open_or_die2(f, 'rb')
            read_handler(IN)
            IN.close()

    sys.exit(0)
