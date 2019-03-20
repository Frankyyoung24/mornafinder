#!/usr/bin/env python
from __future__ import print_function

import argparse
import re
import sys

from port import substr, tr


usage = '''
{} file_fasta seq_adapter

Removes 3' adapters from deep sequencing reads. Works in sequence space.
'''.format(sys.argv[0])


def remove_adapter(_id, seq, prefix):
    seq = tr(seq, '[acgtun.]', '[ACGTTNN]')
    seq_clipped = None

    pattern = r'(\w+)' + prefix
    m = re.search(pattern, seq)
    if m:
        seq_clipped = m.groups()[0]
    elif substr(seq, 0, 6) == prefix:
        seq_clipped = prefix
    else:
        finish = 0

        while not finish and len(prefix) > 0:
            # ATTR: chop $prefix
            prefix = prefix[:-1]
            mm = re.search(r'(\w+){}$'.format(prefix), seq)
            if mm:
                seq_clipped = mm.groups()[0]
                finish = 1

    if not seq_clipped:
        seq_clipped = seq

    # print ">$id\n$seq_clipped\n";
    print('>{}\n{}'.format(_id, seq_clipped))


def remove_adapters(file_fasta, prefix):
    try:
        FASTA = open(file_fasta, 'rb')
    except IOError:
        print('can not open {}'.format(file_fasta))
        sys.exit(-1)

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
                    remove_adapter(_id, seq, prefix)
                    _id = mm.groups()[0]
                    seq = ''
                    continue

                seq += ll

    remove_adapter(_id, seq, prefix)
    FASTA.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_fasta')
    parser.add_argument('seq_adapter')

    if len(sys.argv) != 3:
        print(usage)
        sys.exit(-1)

    args = parser.parse_args(sys.argv[1:3])
    file_fasta = args.file_fasta
    seq_adapter = args.seq_adapter
    seq_test = "TCGTATGCCGTCTTCTGCTTGT"

    prefix = substr(seq_adapter, 0, 6)
    prefix = tr(prefix, '[acgtun.]', '[ACGTTNN]')
    remove_adapters(file_fasta, prefix)
    sys.exit(0)
