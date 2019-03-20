#!/usr/bin/env python
from __future__ import print_function

import argparse
import getopt
import os
import re
import shutil
import sys
import time

from port import die, open_or_die, pprint, print_stderr


usage = '''
{} file_command_line file_structure rounds_controls

-a   Output progress to screen

'''.format(sys.argv[0])


def perform_controls(_dir, rounds, command_line, options):
    os.mkdir(_dir)
    _round = 1

    while _round <= rounds:
        if options.get('-a') == '':
            print_stderr('{}\r'.format(_round))

        cmd = 'permute_structure.py {} > {}/precursors_permuted.str 2> /dev/null'.format(
            file_structure, _dir)
        # print(cmd)
        os.system(cmd)

        pprint('permutation {}\n\n'.format(_round))
        cmd = '{} 2> /dev/null'.format(command_line)
        os.system(cmd)
        # ret = os.popen(cmd).read().strip()
        # pprint(ret)

        _round += 1

    shutil.rmtree(_dir)

    if options.get('-a') == '':
        print_stderr('controls performed\n\n')


def parse_file_command_line(file_command_line, file_structure, _dir):
    FILE = open_or_die(file_command_line, 'rb',
                       'can not open {}'.format(file_command_line))
    while True:
        line = FILE.readline()
        if not line:
            break

        if re.search(r'(\S+)', line):
            line = line.strip()
            line = re.sub(
                file_structure, '{}/precursors_permuted.str'.format(_dir), line, count=1)
            line = re.sub(r'>.+', '', line, count=1)
            return line

    die('{} is empty\n'.format(file_command_line))


if __name__ == '__main__':
    if len(sys.argv) < 4:
        die(usage)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_command_line', help='command list file')
    parser.add_argument('file_structure', help='structure file')
    parser.add_argument('rounds', help='rounds')
    args = parser.parse_args(sys.argv[1:4])

    file_command_line = args.file_command_line
    file_structure = args.file_structure
    rounds = int(args.rounds)

    opt, argss = getopt.getopt(sys.argv[4:], "a")
    options = dict(opt)

    ltime = long(time.time())
    _dir = 'dir_perform_controls{}'.format(ltime)

    command_line = parse_file_command_line(
        file_command_line, file_structure, _dir)
    if options.get('-a') == '':
        print_stderr('total number of rounds controls={}\n'.format(rounds))

    perform_controls(_dir, rounds, command_line, options)
