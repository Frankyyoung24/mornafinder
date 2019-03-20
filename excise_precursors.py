#!/usr/bin/env python

from __future__ import print_function
import sys
import re
import argparse
import getopt
from port import print_stderr, str_reverse, tr, substr
import os

hash_pos = {}
count_lines = 0
count_excisions = 0
freq_min = 2

usage = '''
{} file_fasta file_arf precursor.coords

This script excises potential miRNA precursor sequences from a genome.
The fasta file should be the relevant genome, and the arf file should
contain the read mappings.
The precursors.coords designates where to write the precursor genome coordinates to.

-a integer   Only excise if the potential mature microRNA is represented
             by a number of reads equal to or higher than the integer
             (default 2).
-b           Output progress to screen'''.format(sys.argv[0])


def insertfeature(db, strand, db_beg, db_end, freq):
    global hash_pos
    hash_pos[db][strand][db_beg][db_end] += freq


def find_freq(freq):
    m = re.search(r'_x(\d+)', freq)
    if m:
        m = m.groups()
        return m[0]
    else:
        print_stderr('Problem with read format\n')
        return 1


def parse_file_arf(file_arf):
    global count_lines

    lines = int(os.popen('cat {} | wc -l'.format(file_arf)).read().strip())

    if options.get('-b') == '':
        print_stderr('reading the mapping file into memory, total lines={}\n'.format(lines))

    try:
        FILENAME = open(file_arf, 'rb')
    except IOError:
        print('Could not open file {}'.format(file_arf))
        sys.exit(-1)

    while True:
        line = FILENAME.read()
        if not line:
            break

        m = re.match(r'^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if m:
            m = m.groups()
            query = m[0]
            query_map_lng = int(m[1])
            query_beg = int(m[2])
            query_end = int(m[3])
            query_seq = m[4]
            db = m[5]
            db_map_lng = int(m[6])
            db_beg = int(m[7])
            db_end = int(m[8])
            db_seq = m[9]
            strand = m[10]
            edits = m[11]
            edit_string = m[12]

            freq = find_req(query)

            insertfeature(db, strand, db_beg, db_end, freq)
            count_lines += 1

            if options.get('b') == '':
                pass

    FILENAME.close()


def find_freq_max_downstream(db, strand, db_beg, db_end):
    global hash_pos

    freq_max = 0
    for pos_beg in range(db_beg + 1, db_end + 70 + 1):
        try:
            hash_pos[db][strand][pos_beg]
        except KeyError:
            pass
        else:
            # Defined
            pos_ends = sorted(hash_pos[db][strand][pos_beg].keys())
            for pos_end in pos_ends:
                freq = hash_pos[db][strand][pos_beg][pos_end]
                if freq > freq_max:
                    freq_max = freq

    return freq_max


def max2(a, b):
    return a if a > b else b


def min2(a, b):
    return a if a < b else b


def excise_position(db_seq, db_lng, excise_beg, excise_end):
    excise_beg_limit = max2(1, excise_beg)
    excise_end_limit = min2(db_lng, excise_end)

    excise_lng = excise_end_limit - excise_beg_limit + 1
    # excise_seq = substr($$db_seq,$excise_beg_limit-1,$excise_lng);
    excise_seq = substr(db_seq, excise_beg_limit - 1, excise_lng)

    return excise_seq


def com(sequence):
    return tr(sequence, 'acgtuACGTU', 'TGCAATGCAA')


def rev(sequence):
    return str_reverse(sequence)


def revcom(sequence):
    return rev(com(sequence))


def print_positions(PF, db, strand, db_seq, db_lng, excise_beg, excise_end):
    global count_excisions

    excise_seq = excise_position(db_seq, db_lng, excise_beg, excise_end)
    if strand == '-':
        excise_seq = revcom(excise_seq)

    # print ">$$db\_$count_excisions\n$excise_seq\n";
    # print PF ">$$db\_$count_excisions\t$$strand\t$$excise_beg\t$$excise_end\n";
    print('>{}_{}\n{}'.format(db, count_excisions, excise_seq))
    PF.write('>{}_{}\t{}\t{}\t{}\n'.format(
        db,
        count_excisions,
        strand,
        excise_beg,
        excise_end
    ))

    count_excisions += 1


def excise(PF, db, db_seq):
    global freq_min
    global hash_pos

    strands = sorted(hash_pos[db])
    for strand in strands:
        db_lng = len(db_seq)
        db_limit = 0

        db_begs = sorted(hash_pos[db][strand].keys())
        for db_beg in db_begs:
            db_ends = sorted(hash_pos[db][strand][db_beg])

            for db_end in db_ends:
                freq = hash_pos[db][strand][db_beg][db_end]
                freq_max_ds = find_freq_max_downstream(db, strand, db_beg, db_end)

                if freq < freq_min or freq < freq_max_ds or db_beg > db_limit:
                    continue

                excise_beg = db_beg - 70
                excise_end = db_end + 20

                # print out in fasta format
                print_positions(PF, db, strand, db_seq, db_lng, excise_beg, excise_end)

                excise_beg = db_beg - 20
                excise_end = db_end + 70

                print_positions(PF, db, strand, db_seq, db_lng, excise_beg, excise_end)

                # the most 3' position that has yet been excised
                db_limit = excise_end


def parse_genome_and_excise(PF, file_fasta):
    try:
        FASTA = open(file_fasta, 'rb')
    except IOError:
        FASTA.close()

    while True:
        line = FASTA.readline().strip()
        if not line:
            break

        m = re.match(r'^>(\S+)(.*)', line)
        if m:
            m = m.groups()
            _id = m[0]
            desc = m[1]
            sequence = ''
            while True:
                ll = FASTA.readline().strip()
                if not ll:
                    break

                mm = re.match(r'^>(\S+)(.*)', ll)
                if mm:
                    mm = mm.groups()
                    excise(PF, _id, sequence)
                    _id = mm[0]
                    desc = mm[1]
                    sequence = ''
                    continue

                sequence += ll

    excise(PF, _id, sequence)
    FASTA.close()


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(usage)
        sys.exit(-1)

    parser = argparse.ArgumentParser()
    parser.add_argument('file_fasta', help=usage)
    parser.add_argument('file_arf', help=usage)
    parser.add_argument('coord_file', help=usage)

    args = parser.parse_args(sys.argv[1:4])
    file_fasta = args.file_fasta
    file_arf = args.file_arf
    coord_file = args.coord_file

    opts, argss = getopt.getopt(sys.argv[4:], 'a:b')
    options = dict(opts)

    try:
        PF = open(coord_file, 'w+')
    except:
        print('cannot create file {}'.format(coord_file))
        sys.exit(-1)

    if options.get('-a'):
        freq_min = int(options.get('-a'))

    if options.get('-b') == '':
        print_stderr('finding lengths of genome contigs\n')

    parse_file_arf(file_arf)

    if options.get('-b') == '':
        print_stderr('reading the genome into memory and excising potential precursors\n')

    parse_genome_and_excise(PF, file_fasta)

    if options.get('-b') == '':
        print_stderr('potential precursors excised\n')

    close(PF)
