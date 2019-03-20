#!/usr/bin/env python

import os
import random
import re
import sys

from port import Nicenumber, die, open_or_die, open_or_die2, tr

counter = 0
hash_seq = {}
hash_num = {}
seq = None

_id = None
tag = None

hint = '''Please check your file for the following issues:

I.  Sequences are allowed only to comprise characters [ACGTNacgtn].
II. Identifiers are not allowed to have withespaces.
'''

usage = '''
Usage:

    {} file

'''.format(sys.argv[0])

lines = 0
aBuffer = None
_id = None
tag = None


def read_handler(handle):
    global counter
    while True:
        rin = handle.readline()
        if not rin:
            break

        rin = rin.strip()

        counter += 1
        m = re.match(r'^\>(\S\S\S)(.+)$', rin)
        if m:
            m = m.groups()
            _id = '{}{}'.format(m[0], m[1])
            tag = m[0]

            if re.search(r'\s+', _id):
                die('Error in line {}: The identifier\n \n{}\n \ncontains white spaces\n\n\nPlease make sure that none of the identifiers contain whitepaces.\nYou could run remove_white_space_in_id.py {} > newfile\nThis will remove everything from the id line after the first whitespace'.format(
                    Nicenumber(counter),
                    rin,
                    in_file
                ))
            elif not re.match(r'^(\S\S\S)_(\d+)_(x\d+)$', _id):
                die('Error in line {}: The identifier\n\n{}\n\nhas to have the format\nname_uniqueNumber_xnumber\n\n\nPlease make sure that all identifiers are unique and have the format described above.'.format(
                    Nicenumber(counter),
                    _id,
                ))
            else:
                mm = re.match(r'^(\S\S\S)_(\d+)_(x\d+)$', _id)
                if mm:
                    mm = mm.groups()

                    try:
                        hash_num[m[1]] += 1
                    except KeyError:
                        hash_num[m[1]] = 1

        else:
            mm = re.match(r'^([A|C|G|T|U|N|a|c|g|t|u|n]{17,})$', rin)
            if mm:
                mm = mm.groups()

                defined = False
                try:
                    hash_seq[tag][seq]
                    defined = True
                except KeyError:
                    pass

                if defined and hash_seq[tag][seq]:
                    die('Error in line {}: The sequence\n\n{}\n\noccures at least twice in sample {} in your reads file.\n\nAt first it occured at line {}\n\nPlease make sure that your reads file only contains unique sequences within each sample.\n'.format(
                        Nicenumber(counter),
                        mm[0],
                        tag,
                        Nicenumber(hash_seq[tag][mm[0]])
                    ))
                else:
                    try:
                        hash_seq[tag]
                    except KeyError:
                        hash_seq[tag] = {}

                    hash_seq[tag][mm[0]] = counter

            else:
                die('Error in line {}: Either the sequence\n\n{}\n\ncontains less than 17 characters or contains characters others than [acgtunACGTUN]\n\n\nPlease make sure that your file only comprises sequences that have at least 17 characters\n\ncontaining letters [acgtunACGTUN]\n'.format(
                    Nicenumber(counter),
                    rin
                ))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        die(usage)

    in_file = sys.argv[1]
    # FILE = open_or_die(in_file, 'rb', "Can't open {}".format(in_file))
    # while True:
    #     aBuffer = FILE.read(4096)
    #     if not aBuffer:
    #         break

    #     lines += int(tr(aBuffer, '\n', ''))
    # FILE.close()
    lines = os.popen("wc -l {}".format(in_file)).read().strip()
    lines = int(re.split(r'\s', lines)[0])

    if lines / 2 > 5000000:
        rhash = {}
        rn = 0
        for i in range(0, 5000000):

            defined = False
            try:
                rhash[rn]
                defined = True
                print i
            except KeyError:
                pass

            while defined and rhash[rn]:
                rn = (2 * int(random.randrange(lines / 2)))
                defined = False
                try:
                    rhash[rn]
                    defined = True
                except KeyError:
                    pass

            rhash[rn] = 1

        IN = open_or_die(in_file, 'rb', 'can not open {}'.format(in_file))
        while True:
            l = IN.readline().strip()
            if not l:
                break

            m = re.match(r'^\>(.+)$', l)
            if m:
                m = m.groups()
                counter += 1
                _id = m[0]
                if re.search(r'\s+', _id):
                    die('Error in line {}: The identifier\n {}\n contains white spaces\n\nPlease make sure that none of the identifiers contain whitepaces.\nYou could run remove_white_space_in_id.py {} > newfile\nThis will remove everything from the id line after the first whitespace'.format(
                        Nicenumber(counter),
                        l,
                        in_file,
                    ))
                elif not re.match(r'^(\S\S\S)_(\d+)_(x\d+)$', _id):
                    die('Error in line {}: The identifier\n\n{}\n\nhas to have the format\nname_uniqueNumber_xnumber\n\n\nPlease make sure that all identifiers are unique and have the format described above.'.format(
                        Nicenumber(counter),
                        _id
                    ))
                else:
                    mm = re.match(r'^(\S\S\S)_(\d+)_(x\d+)$', _id)

                    defined = False
                    try:
                        rhash[counter]
                        defined = True
                        print i
                    except KeyError:
                        pass

                    if mm and defined and rhash[counter]:
                        mm = mm.groups()
                        try:
                            hash_num[mm[1]] += 1
                        except KeyError:
                            hash_num[mm[1]] = 1

                        tag = mm[0]
            else:
                mm = re.match(r'^([A|C|G|T|N|a|c|g|t|n]{17,})$', l)
                if mm:
                    mm = mm.groups()
                    seq = mm[0]

                    defined = False
                    try:
                        rhash[counter]
                        defined = True
                        print i
                    except KeyError:
                        pass

                    if defined and rhash[counter]:

                        defined = False
                        try:
                            hash_seq[tag][seq]
                            defined = True
                            print i
                        except KeyError:
                            pass

                        if defined and hash_seq[tag][seq]:
                            die('Error in line {}: The sequence\n\n{}\n\noccures at least twice in your reads file.\n\nAt first it occured at line {}\n\nPlease make sure that your reads file only contains unique sequences.\n'.format(
                                Nicenumber(counter),
                                mm[0],
                                Nicenumber(hash_seq[seq])
                            ))
                        else:
                            try:
                                hash_seq[tag]
                            except KeyError:
                                hash_seq[tag] = {}
                            hash_seq[tag][seq] = counter
                else:
                    die('Error in line {}: Either the sequence\n\n\n{}\n\ncontains less than 17 characters or contains characters others than [acgtunACGTUN]\n\n\nPlease make sure that your file only comprises sequences that have at least 17 characters\n\ncontaining letters [acgtunACGTUN]\n'.format(
                        Nicenumber(counter),
                        l
                    ))
        IN.close()
        counter = 0

        IN = open_or_die(in_file, 'rb', 'can not open {}\n'.format(in_file))
        while True:
            l = IN.readline()
            if not l:
                break

            m = re.match(r'^\>(\S\S\S)_.+$', l)
            if m:
                m = m.groups()
                tag = m[0]
                counter += 1
            else:

                defined = False
                try:
                    rhash[counter]
                    defined = True
                except KeyError:
                    pass

                if defined and rhash[counter]:
                    pass
                else:
                    mm = re.match(r'^([A|C|G|T|N|a|c|g|t|n]{17,})$', l)
                    if mm:
                        mm = mm.groups()
                        seq = mm[0]
                        defined = False
                        try:
                            hash_seq[tag][seq]
                            defined = True
                        except KeyError:
                            pass

                        if defined and hash_seq[tag][seq]:
                            die('Error in line {}: The sequence\n\n{}\n\noccures at least twice in sample {} in your reads file.\n\nAt first it occured at line {}\n\nPlease make sure that your reads file only contains unique sequences within each sample.\n'.format(
                                Nicenumber(counter),
                                mm[0],
                                tag,
                                Nicenumber(hash_seq[tag][seq])
                            ))
        IN.close()
    else:
        if len(sys.argv) == 1:
            read_handler(sys.stdin)
        else:
            for f in sys.argv[1:]:
                IN = open_or_die2(f, 'rb')
                read_handler(IN)
                IN.close()
