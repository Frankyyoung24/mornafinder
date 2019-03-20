#!/usr/bin/env python
from __future__ import print_function

import re
import sys

from port import Nicenumber, create_hash_key_chain, die, open_or_die2

counter = 0
hash_num = {}
_id = None

hint = '''Please check your file for the following issues:\n
I.  Sequences are allowed only to comprise characters [ACGTNacgtn].
II. Identifiers are not allowed to have withespaces.\n'''


def read_handler(handle):
    global counter
    while True:
        rin = handle.readline().strip()
        if not rin:
            break

        counter += 1
        m = re.match(r'^\>(.+)$', rin)
        if m:
            m = m.groups()
            _id = m[0]
            if re.search(r'\s+', _id):
                rin = rin.strip()
                die('''Error in line {}: The identifier\n{}\ncontains white spaces\n{} You could run remove_white_space_in_id.py inputfile > newfile This will remove everything from the id line after the first whitespace'''.format(
                    Nicenumber(counter),
                    rin,
                    hint
                ))
            else:
                create_hash_key_chain(hash_num, 0, _id)
                hash_num[_id] += 1
        elif not (re.search(r'^([A|C|G|T|N|a|c|g|t|n]+)$', rin) or re.search(r'^\s$', rin)):
            rin = rin.strip()
            die('Error in line {},": The sequence\n{}\ncontains characters others than [acgtnACGTN]\n{}'.format(
                Nicenumber(counter),
                rin,
                hint
            ))


if __name__ == '__main__':

    if (len(sys.argv)) == 1:
        read_handler(sys.stdin)
    else:
        for f in sys.argv[1:]:
            IN = open_or_die2(f, 'rb')
            read_handler(IN)
            IN.close()

    counter = 1
    for value in hash_num.values():
        if value > 1:
            counter += 1

    if counter > 1:
        print('\nError: {} read identifiers are not unique'.format(counter))

    sys.exit(0)
