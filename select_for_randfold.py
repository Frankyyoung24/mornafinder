#!/usr/bin/env python

import argparse
import os
import re
import sys

from port import (create_hash_key_chain, die, hash_sort_key, open_or_die,
                  pprint, print_stderr, revcom, substr, tr)

usage = '''
Usage:

    {} file_signature file_structure

    This script identifies potential precursors whose structure is basically consistent
    with Dicer recognition. The relevant ids are outputted one id per line.

'''.format(sys.argv[0])

hash_desc = {}
hash_seq = {}
hash_struct = {}
hash_mfe = {}
hash_query = {}
hash_comp = {}
hash_bp = {}
db_old = None


def excise_struct(struct, beg, end, strand):
    global db_old
    lng = len(struct)

    # begin can be equal to end if only one nucleotide is excised
    if not(beg <= end):
        print_stderr(
            'begin can not be greater than end for {}\n'.format(db_old))
        sys.exit(0)

    # rarely, permuted combinations of signature and structure cause out of bound excision errors.
    # this happens once appr. every two thousand combinations
    if not(beg <= len(struct)):
        return 0

    # the blast parsed format is 1-indexed, substr is 0-indexed
    sub_struct = substr(struct, beg - 1, end - beg + 1)

    return sub_struct


def fill_structure():
    global hash_struct, db_old, hash_bp
    struct = hash_struct[db_old]
    lng = len(struct)

    # local stack for keeping track of basepairings
    bps = []
    for pos in range(1, lng + 1):
        struct_pos = excise_struct(struct, pos, pos, "+")
        if struct_pos == '(':
            bps.append(pos)

        if struct_pos == ')':
            pos_prev = bps.pop()
            hash_bp[pos_prev] = pos
            hash_bp[pos] = pos_prev


def fill_pri():
    '''
    fills basic specifics on the precursor into the 'comp' hash
    '''
    global hash_seq, hash_struct, hash_mfe, db_old, hash_comp
    seq = hash_seq[db_old]
    struct = hash_struct[db_old]
    mfe = hash_mfe[db_old]
    length = len(seq)

    hash_comp["pri_id"] = db_old
    hash_comp["pri_seq"] = seq
    hash_comp["pri_struct"] = struct
    hash_comp["pri_mfe"] = mfe
    hash_comp["pri_beg"] = 1
    hash_comp["pri_end"] = length


def find_mature_query():
    '''
    finds the query with the highest frequency of reads and returns it
    is used to determine the positions of the potential mature sequence
    '''
    global hash_query
    # queries = sort {$hash_query{$b}{"freq"} <=> $hash_query{$a}{"freq"}}
    # keys %hash_query;
    queries = hash_sort_key(hash_query, lambda x: (
        hash_query[x[0]]['freq'] * -1, x[0]))
    mature_query = queries[0]
    # print_stderr(mature_query)
    return mature_query


def find_positions_query(query):
    global hash_query
    beg = hash_query[query]["db_beg"]
    end = hash_query[query]["db_end"]
    return (beg, end)


def find_strand_query(query):
    global hash_query
    strand = hash_query[query]['strand']
    return strand


def excise_seq(seq, beg, end, strand):
    '''
    excise sub sequence from the potential precursor
    '''
    global db_old

    # begin can be equal to end if only one nucleotide is excised
    if not(beg <= end):
        print_stderr('begin can not greater than end for {}\n'.format(db_old))
        sys.exit(0)

    # rarely, permuted combinations of signature and structure cause out of bound excision errors.
    # this happens once appr. every two thousand combinations
    if not(beg <= len(seq)):
        return 0

    # the blast parsed format is 1-indexed, substr is 0-indexed
    sub_seq = substr(seq, beg - 1, end - beg + 1)

    # if on the minus strand, the reverse complement should be returned
    if strand == "-":
        sub_seq = revcom(sub_seq)

    return sub_seq


def arm_mature(beg, end, strand):
    '''
    tests whether the mature sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor
    '''
    global hash_comp

    # mature and star sequences should alway be on plus strand
    if strand == '-':
        return 0

    # there should be no bifurcations and minimum one base pairing
    struct = excise_seq(hash_comp['pri_struct'], beg, end, strand)
    if struct is not None and re.search(r'^(\(|\.)+$', struct) and re.search(r'\(', struct):
        return 'first'
    elif struct is not None and re.search(r'^(\)|\.)+$', struct) and re.search(r'\)', struct):
        return 'second'

    return 0


def fill_mature():
    '''
    fills specifics on the mature sequence into the 'comp' hash
    '''
    global hash_comp
    mature_query = find_mature_query()
    (mature_beg, mature_end) = find_positions_query(mature_query)
    mature_strand = find_strand_query(mature_query)
    mature_seq = excise_seq(
        hash_comp["pri_seq"], mature_beg, mature_end, mature_strand)
    mature_struct = excise_struct(
        hash_comp["pri_struct"], mature_beg, mature_end, mature_strand)
    mature_arm = arm_mature(mature_beg, mature_end, mature_strand)

    hash_comp["mature_query"] = mature_query
    hash_comp["mature_beg"] = mature_beg
    hash_comp["mature_end"] = mature_end
    hash_comp["mature_strand"] = mature_strand
    hash_comp["mature_struct"] = mature_struct
    hash_comp["mature_seq"] = mature_seq
    hash_comp["mature_arm"] = mature_arm

    # print_stderr(mature_struct)


def find_star():
    '''
    uses the 'bp' hash to find the expected star begin and end positions from the mature positions
    '''
    global hash_comp, hash_bp

    # the -2 is for the overhang
    mature_beg = hash_comp["mature_beg"]
    mature_end = hash_comp["mature_end"] - 2
    mature_lng = mature_end - mature_beg + 1

    # in some cases, the last nucleotide of the mature sequence does not form a base pair,
    # and therefore does not basepair with the first nucleotide of the star sequence.
    # In this case, the algorithm searches for the last nucleotide of the mature sequence
    # to form a base pair. The offset is the number of nucleotides searched
    # through.
    offset_star_beg = 0
    offset_beg = 0

    # the offset should not be longer than the length of the mature sequence, then it
    # means that the mature sequence does not form any base pairs
    while not offset_star_beg and offset_beg < mature_lng:
        key = mature_end - offset_beg
        try:
            if hash_bp[key]:
                offset_star_beg = hash_bp[key]
            else:
                offset_beg += 1
        except KeyError:
            offset_beg += 1

    # when defining the beginning of the star sequence, compensate for the
    # offset
    star_beg = offset_star_beg - offset_beg

    # same as above
    offset_star_end = 0
    offset_end = 0
    while not offset_star_end and offset_end < mature_lng:
        key = mature_beg + offset_end
        try:
            if hash_bp[key]:
                offset_star_end = hash_bp[key]
            else:
                offset_end += 1
        except KeyError:
            offset_end += 1

    # the +2 is for the overhang
    star_end = offset_star_end + offset_end + 2

    return(star_beg, star_end)


def arm_star(beg, end):
    '''
    tests whether the star sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor
    '''
    global hash_comp
    # unless the begin and end positions are plausible, test negative
    if not(beg > 0 and beg <= hash_comp['pri_end'] and end > 0 and end <= hash_comp['pri_end'] and beg <= end):
        return 0

    if hash_comp['mature_arm'] == 'first':
        if not(hash_comp['mature_end'] < beg):
            return 0
    elif hash_comp['mature_arm'] == 'second':
        if not(end < hash_comp['mature_beg']):
            return 0

    struct = excise_seq(hash_comp['pri_struct'], beg, end, '+')
    if re.match(r'^(\(|\.)+$', struct) and re.search(r'\(', struct):
        return 'first'
    elif re.match(r'^(\)|\.)+$', struct) and re.search(r'\)', struct):
        return 'second'

    return 0


def fill_star():
    '''
    fills specifics on the expected star strand into 'comp' hash ('component' hash)
    '''
    global hash_comp

    # if the mature sequence is not plausible, don't look for the star arm
    mature_arm = hash_comp['mature_arm']
    if not mature_arm:
        hash_comp['star_arm'] = 0
        return

    # if the star sequence is not plausible, don't fill into the hash
    (star_beg, star_end) = find_star()
    star_arm = arm_star(star_beg, star_end)

    if not star_arm:
        return

    # excise expected star sequence and structure
    star_seq = excise_seq(hash_comp["pri_seq"], star_beg, star_end, "+")
    star_struct = excise_seq(hash_comp["pri_struct"], star_beg, star_end, "+")

    # fill into hash
    hash_comp["star_beg"] = star_beg
    hash_comp["star_end"] = star_end
    hash_comp["star_seq"] = star_seq
    hash_comp["star_struct"] = star_struct
    hash_comp["star_arm"] = star_arm

    # print_stderr(star_struct)


def test_loop(beg, end):
    '''
    tests the loop positions
    '''
    global hash_comp

    # unless the begin and end positions are plausible, test negative
    if not(beg > 0 and beg <= hash_comp['pri_end'] and end > 0 and end <= hash_comp['pri_end'] and beg <= end):
        return 0

    return 1


def fill_loop():
    '''
    fills specifics on the loop sequence into the 'comp' hash
    '''
    global hash_comp

    # unless both mature and star sequences are plausible, do not look for the
    # loop
    try:
        if not(hash_comp['mature_arm'] and hash_comp['star_arm']):
            return
    except KeyError:
        return

    loop_beg = 0
    loop_end = 0

    if hash_comp['mature_arm'] == 'first':
        loop_beg = hash_comp['mature_end'] + 1
    else:
        loop_end = hash_comp['mature_beg'] - 1

    if hash_comp['star_arm'] == 'first':
        loop_beg = hash_comp['star_end'] + 1
    else:
        loop_end = hash_comp['star_beg'] - 1

    # unless the positions are plausible, do not fill into hash
    if not(test_loop(loop_beg, loop_end)):
        return

    loop_seq = excise_seq(hash_comp['pri_seq'], loop_beg, loop_end, '+')
    loop_struct = excise_struct(
        hash_comp['pri_struct'], loop_beg, loop_end, '+')

    hash_comp["loop_beg"] = loop_beg
    hash_comp["loop_end"] = loop_end
    hash_comp["loop_seq"] = loop_seq
    hash_comp["loop_struct"] = loop_struct

    # print_stderr(loop_seq)


def test_components():
    '''
    tests whether potential mature, star, loop and lower flank parts are identifiable
    '''
    global hash_comp
    try:
        if not(hash_comp["mature_struct"]):
            return 0
    except KeyError:
        return 0

    try:
        if not(hash_comp["star_struct"]):
            return 0
    except KeyError:
        return 0

    try:
        if not(hash_comp["loop_struct"]):
            return 0
    except KeyError:
        return 0

    return 1


def no_bifurcations_precursor():
    '''
    tests whether there are bifurcations in the hairpin

    assembles the potential precursor sequence and structure from the expected Dicer products
    this is the expected biological precursor, in contrast with 'pri_seq' that includes
    some genomic flanks on both sides
    '''
    global hash_comp
    pre_struct = ''
    pre_seq = ''

    if hash_comp['mature_arm'] == 'first':
        pre_struct += str(hash_comp['mature_struct']) + \
            str(hash_comp['loop_struct']) + str(hash_comp['star_struct'])
        pre_seq += str(hash_comp['mature_seq']) + \
            str(hash_comp['loop_seq']) + str(hash_comp['star_seq'])
    else:
        pre_struct += str(hash_comp['star_struct']) + str(
            hash_comp['loop_struct']) + str(hash_comp['mature_struct'])
        pre_seq += str(hash_comp['star_seq']) + \
            str(hash_comp['loop_seq']) + str(hash_comp['mature_seq'])

    # read into hash
    hash_comp['pre_struct'] = pre_struct
    hash_comp['pre_seq'] = pre_seq

    # simple pattern matching checks for bifurcations
    if not(re.match(r'^((\.|\()+..(\.|\))+)$', pre_struct)):
        # print_stderr(pre_struct)
        return 0

    return 1


def bp_duplex():
    '''
    fraction of nts of the mature sequence that are base paired in the duplex
    '''
    global hash_comp
    lng = hash_comp['mature_end'] - hash_comp['mature_beg'] + 1
    duplex_bps = 0
    mature_struct = hash_comp['mature_struct']

    # simple pattern matching
    # ATTR: while($mature_struct=~/(\(|\))/g), Global search
    m = re.findall(r'(\(|\))', mature_struct)
    for mm in m:
        duplex_bps += 1

    # print_stderr(duplex_bps)
    return float(duplex_bps) / lng


def diff_lng():
    '''
    find difference between mature and star lengths
    '''
    global hash_comp
    mature_lng = len(hash_comp["mature_struct"])
    star_lng = len(hash_comp["star_struct"])
    diff_lng = mature_lng - star_lng
    return diff_lng


def pass_filtering_structure():
    '''
    The potential precursor must form a hairpin with miRNA precursor-like characteristics
    '''

    ret = 1

    # potential mature, star, loop and lower flank parts must be identifiable
    if not(test_components()):
        return 0

    # no bifurcations
    if not(no_bifurcations_precursor()):
        ret = 0

    # minimum 14 base pairings in duplex
    if not(bp_duplex() >= 0.6):
        ret = 0

    # not more than 6 nt difference between mature and star length
    if not(-6 < diff_lng() and diff_lng() < 6):
        ret = 0

    return ret


def pass_filtering_initial():
    # test if the structure forms a plausible hairpin
    if not(pass_filtering_structure()):
        return 0

    return 1


def reset_variables():
    global hash_query
    global hash_comp
    global hash_bp

    hash_query = {}
    hash_comp = {}
    hash_bp = {}

    return


def resolve_potential_precursor():
    '''
    dissects the potential precursor in parts by filling hashes, and tests if it passes the
    initial filter and the scoring filter
    binary variable whether the potential precursor is still viable
    '''
    global hash_seq, db_old, hash_struct

    fill_structure()

    fill_pri()

    fill_mature()

    fill_star()

    fill_loop()

    # if db_old == 'CHROMOSOME_I_165':
    # print_stderr(hash_comp, "\n\n")
    # print_stderr(hash_query, "\n\n")
    # print_stderr(hash_bp, '\n')

    if pass_filtering_initial():

        seq = hash_seq[db_old]
        struct = hash_struct[db_old]

        pprint("{}\n".format(db_old))

    reset_variables()


def find_freq(query):
    '''
    #finds the frequency of a given read query from its id.
    '''

    m = re.search(r'x(\d+)', query)
    if m:
        freq = int(m.groups()[0])
        return freq
    else:
        return 0


def parse_file_struct(file_struct):
    global db_old
    FILE_STRUCT = open_or_die(
        file_struct, 'rb', 'can not open file {}\n'.format(file_struct))
    while True:
        line = FILE_STRUCT.readline()
        if not line:
            break

        line = line.strip()

        m = re.match(r'^>(\S+)\s*(.*)', line)
        if m:
            m = m.groups()
            _id = m[0]
            desc = m[1]
            seq = ""
            struct = ""
            mfe = ""

            while True:
                line2 = FILE_STRUCT.readline()
                if not line2:
                    break

                line2 = line2.strip()
                mm = re.match(r'^>(\S+)\s*(.*)', line2)
                if mm:
                    hash_desc[_id] = desc
                    hash_seq[_id] = seq
                    hash_struct[_id] = struct
                    hash_mfe[_id] = mfe
                    _id = mm.groups()[0]
                    desc = mm.groups()[1]
                    seq = ""
                    struct = ""
                    mfe = ""
                    continue

                m3 = re.match(r'^\w', line2)
                if m3:
                    line2 = tr(line2, 'uU', 'tT')
                    seq += line2

                m3 = re.search(r'((\.|\(|\))+)', line2)
                if m3:
                    struct += m3.groups()[0]

                m3 = re.search(r'\((\s*-\d+\.\d+)\)', line2)
                if m3:
                    mfe = m3.groups()[0]

    hash_desc[_id] = desc
    hash_seq[_id] = seq
    hash_struct[_id] = struct
    hash_mfe[_id] = mfe

    # print('\n'.join(sorted(hash_struct.values())))
    # print('\n'.join(sorted(hash_desc.keys())))

    FILE_STRUCT.close()


def parse_file_arf(file_arf):
    '''
    read through the signature blastparsed file, fills up a hash with information on queries
    (deep sequences) mapping to the current db (potential precursor) and resolve each
    potential precursor in turn
    '''
    global db_old, hash_query
    FILENAME = open_or_die(
        file_arf, 'rb', 'could not open file {}\n'.format(file_arf))
    while True:
        line = FILENAME.readline()
        if not line:
            break

        m = re.match(
            r'^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if m:
            m = m.groups()

            query = m[0]
            query_lng = int(m[1])
            query_beg = int(m[2])
            query_end = int(m[3])
            query_seq = m[4]
            db = m[5]
            db_lng = int(m[6])
            db_beg = int(m[7])
            db_end = int(m[8])
            db_seq = m[9]
            strand = m[10]
            edits = int(m[11])
            edit_string = m[12]

            # only reads that map sense to the potential precursor are
            # considered
            if strand == "-":
                continue

            # if the new line concerns a new db (potential precursor) then the
            # old db must be resolved
            if db_old and db_old != db:
                # print(db, db_old)
                resolve_potential_precursor()

            # resolve the number of reads that the deep sequence represents
            freq = find_freq(query)

            # read information of the query (deep sequence) into hash
            create_hash_key_chain(hash_query, db_beg, query, 'db_beg')
            create_hash_key_chain(hash_query, db_end, query, 'db_end')
            create_hash_key_chain(hash_query, strand, query, 'strand')
            create_hash_key_chain(hash_query, freq, query, 'freq')

            hash_query[query]["db_beg"] = db_beg
            hash_query[query]["db_end"] = db_end
            hash_query[query]["strand"] = strand
            hash_query[query]["freq"] = freq

            db_old = db

    resolve_potential_precursor()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        die(usage)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_signature', help='signature file')
    parser.add_argument('file_structure', help='structure file')
    args = parser.parse_args(sys.argv[1:3])

    file_signature = args.file_signature
    file_structure = args.file_structure

    # parse structure file outputted from RNAfold
    parse_file_struct(file_structure)

    # parse signature file in arf format and resolve each potential precursor
    parse_file_arf(file_signature)
