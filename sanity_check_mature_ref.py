#!/usr/bin/env python

import re
import sys

from port import Nicenumber, create_hash_key_chain, die, open_or_die2

counter = 0


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
                die('Error in line {}: The identifier\n {}\n\ncontains white spaces\n\n{}\n\nYou could run remove_white_space_in_id.py inputfile > newfile\nThis will remove everything from the id line after the first whitespace\n'.format(
                    Nicenumber(counter),
                    _id,
                    hint
                ))
            else:
                create_hash_key_chain(hash_num, 0, _id)
                hash_num[_id] += 1
        elif not re.match(r'^([A|C|G|T|U|N|a|c|g|t|u|n]+)$', rin):
            die('Error in line {}: The sequence\n{}\n\ncontains characters others than [acgtunACGTUN]\n\n{}'.format(
                Nicenumber(counter),
                rin,
                hint
            ))


if __name__ == '__main__':
    hash_num = {}
    _id = None

    hint = 'Please check your file for the following issues:\n\nI.  Sequences are allowed only to comprise characters [ACGTNacgtn].\nII. Identifiers are not allowed to have withespaces.\n'

    if len(sys.argv) == 1:
        # from stdin
        read_handler(sys.stdin)
    else:
        # from files
        for f in sys.argv[1:]:
            IN = open_or_die2(f, 'rb')

            read_handler(IN)

            IN.close()

    sys.exit(0)
