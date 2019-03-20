#!/usr/bin/env python
from __future__ import print_function

import argparse
import getopt
import os
import re
import sys

from port import (create_hash_key_chain, esplit, max2, min2, open_or_die,
                  print_stderr, str_reverse, substr, tr)

usage = '''{} file_fasta file_arf precursor.coords file_output pres_max

This script excises potential miRNA precursor sequences from a genome.
The fasta file should be the relevant genome, and the arf file should
contain the read mappings.
The precursors.coords designates where to write the precursor genome coordinates to.

-b           Output progress to screen'''.format(sys.argv[0])

thres_counts = {}
upper_bound = 50
dblimit = {}
hash_db_lng = {}
hash_pos = {}
count_lines = 0
count_excisions = 0
freq_min = 1


def insertfeature(db, strand, db_beg, db_end, freq):
    global hash_pos
    create_hash_key_chain(hash_pos, 0, db, strand, db_beg, db_end)
    hash_pos[db][strand][db_beg][db_end] += freq


def find_freq(query):
    m = re.search(r'_x(\d+)', query)
    if m:
        return int(m.groups()[0])
    else:
        print_stderr('Problem with read format\n')
        return 1


def parse_file_arf(file_arf):
    global count_lines, hash_pos

    lines = int(os.popen('cat {} | wc -l'.format(file_arf)).read().strip())

    if options.get('-b') == '':
        print_stderr(
            'reading the mapping file into memory, total lines=$lines\n'.format(lines))

    FILENAME = open_or_die(
        file_arf, 'rb', 'Could not open file {}'.format(file_arf))

    while True:
        line = FILENAME.readline()
        if not line:
            break

        m = re.match(
            r'^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
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
            edits = int(m[11])
            edit_string = m[12]

            freq = find_freq(query)
            # read into position hash
            insertfeature(db, strand, db_beg, db_end, freq)

            count_lines += 1

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


def rev(sequence):
    return str_reverse(sequence)


def com(sequence):
    return tr(sequence, 'acgtuACGTU', 'TGCAATGCAA')


def revcom(sequence):
    return rev(com(sequence))


def excise_position(db_seq, db_lng, excise_beg, excise_end):
    excise_beg_limit = max2(1, excise_beg)
    excise_end_limit = min2(db_lng, excise_end)
    excise_lng = excise_end_limit - excise_beg_limit + 1
    excise_seq = substr(db_seq, excise_beg_limit - 1, excise_lng)

    return excise_seq


def print_positions(TMP1, TMP2, db, strand, db_seq, db_lng, excise_beg, excise_end, freq, thres):
    global count_excisions

    excise_seq = excise_position(db_seq, db_lng, excise_beg, excise_end)
    if strand == '-':
        excise_seq = revcom(excise_seq)

    # print TMP1 ">$$db\_xyz123_${count_excisions}_xyz123_$${freq}_xyz123_$thres\n$excise_seq\n";
    # print TMP2
    # ">$$db\_xyz123_${count_excisions}_xyz123_$${freq}_xyz123_$thres\t$$strand\t$$excise_beg\t$$excise_end\n";
    TMP1.write('>{}_xyz123_{}_xyz123_{}_xyz123_{}\n{}\n'.format(
        db,
        count_excisions,
        freq,
        thres,
        excise_seq
    ))

    TMP2.write('>{}_xyz123_{}_xyz123_{}_xyz123_{}\t{}\t{}\t{}\n'.format(
        db,
        count_excisions,
        freq,
        thres,
        strand,
        excise_beg,
        excise_end
    ))

    count_excisions += 1


def excise(TMP1, TMP2, db, db_seq):
    global hash_pos
    global dblimit
    global freq_min
    global thres_counts

    try:
        strands = sorted(hash_pos[db].keys())
    except KeyError:
        strands = []

    for strand in strands:
        db_lng = len(db_seq)
        db_limit = 0
        for z in range(1, upper_bound):
            dblimit[z] = 0

        db_begs = sorted(hash_pos[db][strand].keys())
        for db_beg in db_begs:
            # reads with the same 5' position can have different 3' positions
            db_ends = sorted(hash_pos[db][strand][db_beg].keys())
            for db_end in db_ends:
                # height of read stack with those exact 5' and 3' positions
                freq = hash_pos[db][strand][db_beg][db_end]

                # what is the highest read stack downstream (up to 70 nt
                # downstream)?
                freq_max_ds = find_freq_max_downstream(
                    db, strand, db_beg, db_end)

                # if read stack to low, if higher read stack downstream or if this locus has already
                # been excised, then continue to next stack

                for z in range(1, upper_bound):
                    freq_min = z
                    if freq < freq_min or freq < freq_max_ds or db_beg < dblimit[z]:
                        continue

                    # else excise to sequences, corresponding to the read stack being the mature sequence
                    # in the 5' arm or the 3' arm of the precursor hairpin
                    excise_beg = db_beg - 70
                    excise_end = db_end + 20

                    # correction if beg pos is negative
                    if excise_beg < 0:
                        excise_beg = 0

                    # print out in fasta format
                    print_positions(TMP1, TMP2, db, strand, db_seq,
                                    db_lng, excise_beg, excise_end, freq, z)

                    excise_beg = db_beg - 20
                    excise_end = db_end + 70

                    # correction if beg pos is negative
                    if excise_beg < 0:
                        excise_beg = 0

                    print_positions(TMP1, TMP2, db, strand, db_seq,
                                    db_lng, excise_beg, excise_end, freq, z)

                    dblimit[z] = excise_end
                    thres_counts[z] += 2


def parse_genome_and_excise(TMP1, TMP2, file_fasta):
    FASTA = open_or_die(file_fasta, 'rb', 'can not open {}'.format(file_fasta))

    while True:
        line = FASTA.readline()
        if not line:
            break

        line = line.strip()

        m = re.match(r'^>(\S+)(.*)', line)
        if m:
            _id = m.groups()[0]
            desc = m.groups()[1]
            sequence = ''
            while True:
                ll = FASTA.readline()
                if not ll:
                    break
                ll = ll.strip()

                mm = re.match(r'^>(\S+)(.*)', ll)
                if mm:
                    excise(TMP1, TMP2, _id, sequence)
                    _id = mm.groups()[0]
                    desc = mm.groups()[1]
                    sequence = ''
                    continue
                sequence += ll

    excise(TMP1, TMP2, _id, sequence)
    FASTA.close()


def make_files(file_output, coord_file, thres):
    TMP1 = open_or_die(
        '{}_all'.format(file_output),
        'rb',
        'can not open {}'.format('file_output')
    )

    PRES = open_or_die(file_output, 'w+',
                       'can not create '.format(file_output))
    counter = 0

    while True:
        l = TMP1.readline()
        if not l:
            break

        l = l.strip()
        if re.match('^>', l):
            line = l.split('_xyz123_')
            if int(line[3]) == thres:  # last value
                counter += 1
                seq = TMP1.readline()
                PRES.write('{}_{}\n{}'.format(line[0], counter, seq))

    TMP1.close()
    PRES.close()

    COORD = open_or_die(
        coord_file, 'w+', 'can not create {}'.format(coord_file))
    TMP2 = open_or_die('{}_all'.format(coord_file), 'rb',
                       'can not open {}'.format(coord_file))
    counter = 0

    while True:
        l = TMP2.readline()
        if not l:
            break

        l = l.strip()

        if re.match('^>', l):
            # split() => split by space, \t, \n
            line = esplit(l)
            line2 = line[0].split('_xyz123_')
            if int(line2[3]) == thres:
                counter += 1
                COORD.write('{}_{}\t{}\t{}\t{}\n'.format(
                    line2[0],
                    counter,
                    line[1],
                    line[2],
                    line[3]
                ))

    TMP2.close()
    COORD.close()


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print(usage)
        sys.exit(-1)

    parser = argparse.ArgumentParser()
    parser.add_argument('file_fasta', help=usage)
    parser.add_argument('file_arf', help=usage)
    parser.add_argument('file_output', help=usage)
    parser.add_argument('coord_file', help=usage)
    parser.add_argument('pres_max', help=usage)

    args = parser.parse_args(sys.argv[1:6])
    file_fasta = args.file_fasta
    file_arf = args.file_arf
    file_output = args.file_output
    coord_file = args.coord_file
    pres_max = args.pres_max

    opts, argss = getopt.getopt(sys.argv[6:], 'b')
    options = dict(opts)

    if not re.search(r'^[-]*\d+', pres_max):
        print_stderr('{} is not an integer number\n'.format(pres_max))
        sys.exit(-1)

    for z in range(1, upper_bound):
        dblimit[z] = 0
        thres_counts[z] = 0

    TMP1 = open_or_die(
        '{}_all'.format(file_output),
        'w+',
        'cannot create file {}'.format(file_output))

    TMP2 = open_or_die(
        '{}_all'.format(coord_file),
        'w+',
        'cannot create file {}'.format(coord_file))

    if options.get('-b') == '':
        print_stderr('finding lengths of genome contigs\n')

    parse_file_arf(file_arf)

    if options.get('-b') == '':
        print_stderr(
            'reading the genome into memory and excising potential precursors\n')

    parse_genome_and_excise(TMP1, TMP2, file_fasta)

    if options.get('-b') == '':
        print_stderr('')

    TMP1.close()
    TMP2.close()

    for z in range(1, upper_bound):
        print_stderr('{}\t{}\n'.format(z, thres_counts[z]))
        if thres_counts[z] < pres_max or pres_max < 0:
            if options.get('-b') == '':
                print_stderr('creating output file now\n')

            make_files(file_output, coord_file, z)
            OSS = open_or_die('{}_stack'.format(file_output),
                              'w+', 'File could not be created\n')
            OSS.write('{}\n'.format(z))
            OSS.close()
            sys.exit(0)

    sys.exit(0)
