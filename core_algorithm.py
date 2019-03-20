#!/usr/bin/env python
from __future__ import print_function

import argparse
import getopt
import math
import re
import sys

from port import (create_hash_key_chain, die, hash_defined, hash_sort_key,
                  hash_value_true, max2, min2, open_or_die, pprint,
                  print_stderr, pround, revcom, ssplit, substr, tr)

usage = '''
{} file_signature file_structure

This is the core algorithm of moRNA Finder. It takes as input a file in blastparsed format with
information on the positions of reads aligned to potential precursor sequences (signature).
It also takes as input an RNAfold output file, giving information on the sequence, structure
and mimimum free energy of the potential precursor sequences.

Extra arguments can be given. -s specifies a fastafile containing the known mature miRNA
sequences that should be considered for conservation purposes. -t prints out the potential
precursor sequences that do _not_ exceed the cut-off (default prints out the sequences that
exceeds the cut-off). -u gives limited output, that is only the ids of the potential precursors
that exceed the cut-off. -v varies the cut-off. -x is a sensitive option for Sanger sequences
obtained through conventional cloning. -z consider the number of base pairings in the lower
stems (this option is not well tested).

-h print this usage
-s fasta file with known miRNAs
-t print filtered
-u limited output (only ids)
-v cut-off (default 1)
-x sensitive option for Sanger sequences
-y file with randfold p-values
-z consider Drosha processing

'''.format(sys.argv[0])


# Global variables

# parameters
seed_lng = 7
mature_lng_max = 25
score_star = 3.9
score_star_not = -1.3
score_seed = 3
score_seed_not = -0.6
score_randfold = 1.6
score_randfold_not = -2.2
scores_stem = (-3.1, -2.3, -2.2, -1.6, -1.5, 0.1, 0.6, 0.8, 0.9, 0.9, 0)
score_min = 1
e = 2.718281828

hash_desc = {}
hash_seq = {}
hash_struct = {}
hash_mfe = {}
hash_seeds = {}
hash_mirs = {}
hash_randfold = {}
hash_query = {}
hash_comp = {}
hash_bp = {}
hash_stars = {}

db_old = None
lines = ''
out_of_bound = None


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
    global hash_seq, hash_struct, hash_mfe, db_old
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
    global db_old, out_of_bound

    # begin can be equal to end if only one nucleotide is excised
    if not(beg <= end):
        print_stderr('begin can not greater than end for {}\n'.format(db_old))
        sys.exit(0)

    # rarely, permuted combinations of signature and structure cause out of bound excision errors.
    # this happens once appr. every two thousand combinations
    if not(beg <= len(seq)):
        out_of_bound += 1
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
        if key in hash_bp.keys() and hash_bp[key]:
            offset_star_beg = hash_bp[key]
        else:
            offset_beg += 1

    # when defining the beginning of the star sequence, compensate for the
    # offset
    star_beg = offset_star_beg - offset_beg

    # same as above
    offset_star_end = 0
    offset_end = 0
    while not offset_star_end and offset_end < mature_lng:
        key = mature_beg + offset_end
        if key in hash_bp.keys() and hash_bp[key]:
            offset_star_end = hash_bp[key]
        else:
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
    if re.search(r'^(\(|\.)+$', struct) and re.search(r'\(', struct):
        return 'first'
    elif re.search(r'^(\)|\.)+$', struct) and re.search(r'\)', struct):
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
    (star_beg_exp, star_end_exp) = find_star()
    star_arm = arm_star(star_beg_exp, star_end_exp)

    if not star_arm:
        return

    # excise expected star sequence and structure
    star_seq = excise_seq(hash_comp["pri_seq"],
                          star_beg_exp, star_end_exp, "+")
    star_struct = excise_seq(
        hash_comp["pri_struct"], star_beg_exp, star_end_exp, "+")

    # fill into hash
    hash_comp["star_beg_exp"] = star_beg_exp
    hash_comp["star_end_exp"] = star_end_exp
    hash_comp["star_seq"] = star_seq
    hash_comp["star_struct"] = star_struct
    hash_comp["star_arm"] = star_arm


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

    loop_beg = None
    loop_end = None

    if hash_comp['mature_arm'] == 'first':
        loop_beg = hash_comp['mature_end'] + 1
    else:
        loop_end = hash_comp['mature_beg'] - 1

    if hash_comp['star_arm'] == 'first':
        loop_beg = hash_comp['star_end_exp'] + 1
    else:
        loop_end = hash_comp['star_beg_exp'] - 1

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


def test_components():
    '''
    tests whether potential mature, star, loop and lower flank parts are identifiable
    '''
    global hash_comp

    try:
        if not(hash_comp["mature_struct"]):
            hash_comp[
                "problem_struct_mature"] = "no mature sequence could be identified\n"
            return 0
    except KeyError:
        return 0

    try:
        if not(hash_comp["star_struct"]):
            hash_comp[
                "problem_struct_star"] = "the candidate mature sequence does not form part of a likely miRNA duplex\n"
            return 0
    except KeyError:
        return 0

    try:
        if not(hash_comp["loop_struct"]):
            hash_comp[
                "problem_struct_loop"] = "no loop sequence could be identified\n"
            return 0
    except KeyError:
        return 0

    try:
        if not(hash_comp["flank_first_struct"]):
            hash_comp[
                "problem_struct_us_flanks"] = "no upstream flanking sequence could be identified\n"
            return 0
    except KeyError:
        return 0

    try:
        if not(hash_comp["flank_second_struct"]):
            hash_comp[
                "problem_struct_ds_flanks"] = "no downstream flanking sequence could be identified\n"
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
    if not(re.search(r'^((\.|\()+..(\.|\))+)$', pre_struct)):
        hash_comp[
            "problem_struct_bifurcation"] = "there are bifurcations in the candidate precursor structure\n"
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
    global hash_comp

    # return value
    ret = 1

    # potential mature, star, loop and lower flank parts must be identifiable
    if not(test_components()):
        return 0

    # no bifurcations
    if not(no_bifurcations_precursor()):
        ret = 0

    # minimum 60% base pairings in duplex
    if not(bp_duplex() >= 0.6):
        ret = 0
        hash_comp[
            "problem_duplex_bp"] = "less than 60% of the nts in the candidate mature seq are base paired in the duplex\n"

    # not more than 6 nt difference between mature and star length
    if not(-6 < diff_lng() and diff_lng() < 6):
        ret = 0
        hash_comp[
            "problem_duplex_lng"] = "the difference between the candidate mature and star sequence is > 6 nts\n"

    return ret


def testbeginend(begin1, end1, begin2, end2):
    '''
    #Are the beginposition numerically smaller than the endposition for each pair?
    '''
    global db_old
    if not(begin1 <= end1 and begin2 <= end2):
        print_stderr("beg can not be larger than end for {}\n".format(db_old))
        sys.exit(0)


def contained(beg1, end1, beg2, end2):
    '''
    # Is the stretch defined by the first positions contained in the stretch defined by the second?
    '''
    testbeginend(beg1, end1, beg2, end2)

    if beg2 <= beg1 and end1 <= end2:
        return 1
    else:
        return 0


def test_query(query):
    '''
    test if deep sequence maps in consistence with Dicer processing
    '''
    global hash_query, hash_comp

    # begin, end, strand and read count
    beg = hash_query[query]["db_beg"]
    end = hash_query[query]["db_end"]
    strand = hash_query[query]["strand"]
    freq = hash_query[query]["freq"]

    # should not be on the minus strand (although this has in fact anecdotally
    # been observed for known miRNAs)
    if strand == '-':
        return 0

    # the deep sequence is allowed to stretch 2 nt beyond the expected 5' end
    fuzz_beg = 2
    # the deep sequence is allowed to stretch 5 nt beyond the expected 3' end
    fuzz_end = 5

    # if in accordance with Dicer processing, return the type of Dicer product
    if contained(beg, end, hash_comp["mature_beg"] - fuzz_beg, hash_comp["mature_end"] + fuzz_end):
        return "mature"
    if contained(beg, end, hash_comp["star_beg_exp"] - fuzz_beg, hash_comp["star_end_exp"] + fuzz_end):
        return "star"
    if contained(beg, end, hash_comp["loop_beg"] - fuzz_beg, hash_comp["loop_end"] + fuzz_end):
        return "loop"

    # if not in accordance, return 0
    return 0


def test_star(query):
    '''
    test if a deep sequence maps in good consistence with 3' overhang
    '''
    global hash_query, hash_comp, hash_stars
    # 5' begin and 3' end positions
    beg = hash_query[query]["db_beg"]
    end = hash_query[query]["db_end"]
    freq = hash_query[query]["freq"]

    create_hash_key_chain(hash_stars, 0, beg, end)
    hash_stars[beg][end] += freq

    # the difference between observed and expected begin positions must be 0
    # or 1
    offset = beg - hash_comp["star_beg_exp"]
    if offset == 0 or offset == 1 or offset == -1:
        return 1

    return 0


def observed_star():
    global hash_stars, hash_comp

    beg_max = None
    end_max = None
    freq_max = 0

    begs = sorted(hash_stars.keys())
    for beg in begs:
        ends = sorted(hash_stars[beg], reverse=True)
        for end in ends:
            freq = hash_stars[beg][end]
            if freq_max < freq:
                beg_max = beg
                end_max = end
                freq_max = freq

    hash_comp["star_beg_obs"] = beg_max
    hash_comp["star_end_obs"] = end_max
    hash_comp["star_freq_consensus"] = freq_max

    return


def pass_filtering_signature():
    global hash_comp, hash_query

    # if the putative mature sequence is longer than the designated number of
    # nts, discard
    mature_lng = hash_comp["mature_end"] - hash_comp["mature_beg"] + 1
    if(mature_lng > mature_lng_max):
        hash_comp[
            "problem_mature_lng"] = "the candidate mature seq is > $mature_lng_max nts\n"
        return 0

    # number of reads that map in consistence with Dicer processing
    consistent = 0

    # number of reads that map inconsistent with Dicer processing
    inconsistent = 0

    #   number of potential star reads map in good consistence with Drosha/Dicer processing
    #   (3' overhangs relative to mature product)
    star_perfect = 0

    # number of potential star reads that do not map in good consistence with
    # 3' overhang
    star_fuzzy = 0

    freq_mature = 0
    freq_loop = 0
    freq_star = 0

    # sort queries (deep sequences) by their position on the hairpin
    queries = hash_sort_key(hash_query, lambda x: hash_query[x[0]]['db_beg'])

    for query in queries:

        # number of reads that the deep sequence represents
        if not(hash_defined(hash_query, query, "freq")):
            continue

        query_freq = hash_query[query]["freq"]

        # test which Dicer product (if any) the deep sequence corresponds to
        product = test_query(query)

        # and add the appropriate read counts
        if product == "mature":
            freq_mature += query_freq
        if product == "loop":
            freq_loop += query_freq
        if product == "star":
            freq_star += query_freq

        if product:
            # if the deep sequence corresponds to a Dicer product, add to the
            # 'consistent' variable
            consistent += query_freq
        else:
            # if the deep sequence do not correspond to a Dicer product, add to
            # the 'inconsistent' variable
            inconsistent += query_freq

        # test a potential star sequence has good 3' overhang
        if product == "star":
            if test_star(query):
                star_perfect += query_freq
            else:
                star_fuzzy += query_freq

    #   if the majority of potential star sequences map in good accordance with 3' overhang
    #    score for the presence of star evidence
    if star_perfect > star_fuzzy:
        hash_comp["star_read"] = 1

    # find the most frequent star positions (this is opposed to star_expm which are the
    # positions that would be expected from the positions of the mature sequence and
    # the model of Dicer processing
    if 0 < star_perfect + star_fuzzy:
        observed_star()

    # number of reads that correspond to mature, loop, star and the sum of
    # these
    hash_comp["freq_mature"] = freq_mature
    hash_comp["freq_loop"] = freq_loop
    hash_comp["freq_star"] = freq_star
    hash_comp["freq_total"] = consistent

    if not(consistent + inconsistent) > 0:
        hash_comp["problem_read_cnt"] = "no reads map to the candidate\n"
        return 0

    # unless >90% of the reads map in consistence with Dicer processing, the
    # hairpin is discarded
    inconsistent_fraction = float(inconsistent) / (inconsistent + consistent)
    if not(inconsistent_fraction <= 0.1):
        hash_comp[
            "problem_Dicer_signature"] = "more than 10% of the reads map inconsistent with Dicer processing\n"
        return 0

    # the hairpin is retained
    return 1


def pass_filtering_initial():
    global hash_comp

    # test if the structure forms a plausible hairpin
    if not(pass_filtering_structure()):
        hash_comp["problem_structure"] = "The candidate has been discarded because the structure is inconsistent with Dicer processing:\n"
        return 0

    # test if >90% of reads map to the hairpin in consistence with Dicer
    # processing
    if not(pass_filtering_signature()):
        hash_comp["problem_signature"] = "The candidate has been discarded because the signature is inconsistent with Dicer processing:\n"
        return 0

    return 1


def print_hash_comp(options):
    global hash_comp
    reads = []
    lr = None
    bef = None
    after = None
    col1_width = 40

    if options.get('-l'):
        col1_width = int(options.get('-l'))

    col1_width += 5
    rseq = None

    hash_comp["pri_seq"] = tr(hash_comp["pri_seq"], 'Tt', 'Uu')
    pseq = ssplit(hash_comp['pri_seq'].lower())

    spacer = " " * (col1_width - len("pri_seq"))
    shift = None
    reads_hash = {}

    pprint(">{}\n".format(hash_comp['pri_id']))
    spacer2 = 14
    spacer2s = 0

    # pprint('hi', hash_comp["problem_structure"], hash_comp["problem_signature"], '\n')
    problem = False
    try:
        if (not hash_comp["problem_structure"] and not hash_comp["problem_signature"]):
            problem = True
    except KeyError:
        hash_comp["score"] = int(hash_comp["score"]) if hash_comp[
            "score"] == int(hash_comp["score"]) else hash_comp["score"]
        spacer2s = " " * (spacer2 - len(str(hash_comp["score"])))
        pprint("score total\t\t{}{}\n".format(spacer2s, hash_comp['score']))

        hash_comp["score_star"] = int(hash_comp["score_star"]) if hash_comp[
            "score_star"] == int(hash_comp["score_star"]) else hash_comp["score_star"]
        spacer2s = " " * (spacer2 - len(str(hash_comp["score_star"])))
        pprint("score for star read(s)\t{}{}\n".format(
            spacer2s, hash_comp['score_star']))

        hash_comp["score_freq"] = int(hash_comp["score_freq"]) if hash_comp[
            "score_freq"] == int(hash_comp["score_freq"]) else hash_comp["score_freq"]
        spacer2s = " " * (spacer2 - len(str(hash_comp["score_freq"])))
        pprint("score for read counts\t{}{}\n".format(
            spacer2s, hash_comp['score_freq']))

        hash_comp["score_mfe"] = int(hash_comp["score_mfe"]) if hash_comp[
            "score_mfe"] == int(hash_comp["score_mfe"]) else hash_comp["score_mfe"]
        spacer2s = " " * (spacer2 - len(str(hash_comp["score_mfe"])))
        pprint("score for mfe\t\t{}{}\n".format(
            spacer2s, hash_comp['score_mfe']))

        if options.get('-y'):
            spacer2s = " " * (spacer2 - len(str(hash_comp["score_randfold"])))
            pprint("score for randfold\t{}{}\n".format(
                spacer2s, hash_comp['score_randfold']))

        if options.get('-s'):
            spacer2s = " " * (spacer2 - len(str(hash_comp["score_seed"])))
            pprint("score for cons. seed\t{}{}\n".format(
                spacer2s, hash_comp['score_seed']))

        if options.get('-s') and hash_defined(hash_comp, "seed_family"):
            spacer2s = " " * (spacer2 - len(str(hash_comp["seed_family"])))
            pprint("miRNA with same seed\t{}{}\n".format(
                spacer2s, hash_comp['seed_family']))

        spacer2s = " " * (spacer2 - len(str(hash_comp["freq_total"])))
        pprint("total read count\t{}{}\n".format(
            spacer2s, hash_comp['freq_total']))
        spacer2s = " " * (spacer2 - len(str(hash_comp["freq_mature"])))
        pprint("mature read count\t{}{}\n".format(
            spacer2s, hash_comp['freq_mature']))
        spacer2s = " " * (spacer2 - len(str(hash_comp["freq_loop"])))
        pprint("loop read count\t\t{}{}\n".format(
            spacer2s, hash_comp['freq_loop']))
        spacer2s = " " * (spacer2 - len(str(hash_comp["freq_star"])))
        pprint("star read count\t\t{}{}\n".format(
            spacer2s, hash_comp['freq_star']))

    if problem:
        # structure problems
        if hash_defined(hash_comp, "problem_structure"):
            pprint(hash_comp["problem_structure"])
        if hash_defined(hash_comp, "problem_mature_lng"):
            pprint(hash_comp["problem_mature_lng"])
        if hash_defined(hash_comp, "problem_duplex_bp"):
            pprint(hash_comp["problem_duplex_bp"])
        if hash_defined(hash_comp, "problem_duplex_lng"):
            pprint(hash_comp["problem_duplex_lng"])
        if hash_defined(hash_comp, "problem_struct_mature"):
            pprint(hash_comp["problem_struct_mature"])
        if hash_defined(hash_comp, "problem_struct_star"):
            pprint(hash_comp["problem_struct_star"])
        if hash_defined(hash_comp, "problem_struct_loop"):
            pprint(hash_comp["problem_struct_loop"])
        if hash_defined(hash_comp, "problem_struct_us_flanks"):
            pprint(hash_comp["problem_struct_us_flanks"])
        if hash_defined(hash_comp, "problem_struct_ds_flanks"):
            pprint(hash_comp["problem_struct_ds_flanks"])
        if hash_defined(hash_comp, "problem_struct_bifurcation"):
            pprint(hash_comp["problem_struct_bifurcation"])

        # signature problems
        if hash_defined(hash_comp, "problem_signature"):
            pprint(hash_comp["problem_signature"])
        if hash_defined(hash_comp, "problem_read_cnt"):
            pprint(hash_comp["problem_read_cnt"])
        if hash_defined(hash_comp, "problem_Dicer_signature"):
            pprint(hash_comp["problem_Dicer_signature"])

    # printing alignment;
    pprint("exp")

    pprint(" " * (col1_width - 3))

    problem = False
    try:
        if hash_comp["problem_structure"]:
            problem = True
    except KeyError:
        pass

    if problem:
        pprint("n" * len(hash_comp["pri_seq"]))
    else:
        pprint("f" * hash_comp["flank_first_end"])
        if hash_comp["mature_beg"] < hash_comp["loop_beg"]:
            bef = "M" * (hash_comp["mature_end"] - hash_comp["mature_beg"] + 1)
            pprint(bef)
            bef = "l" * (hash_comp["loop_end"] - hash_comp["loop_beg"] + 1)
            pprint(bef)
            bef = "S" * (hash_comp["star_end_exp"] -
                         hash_comp["star_beg_exp"] + 1)
            pprint(bef)
        else:
            bef = "S" * (hash_comp["star_end_exp"] -
                         hash_comp["star_beg_exp"] + 1)
            pprint(bef)
            bef = "l" * (hash_comp["loop_end"] - hash_comp["loop_beg"] + 1)
            pprint(bef)

            bef = "M" * (hash_comp["mature_end"] - hash_comp["mature_beg"] + 1)
            pprint(bef)

        bef = "f" * (hash_comp["pri_end"] - hash_comp["flank_second_beg"] + 1)
        pprint(bef)

    if hash_defined(hash_comp, 'star_beg_obs'):
        pprint('\nobs')

        pprint(" " * (col1_width - 3))

        problem = False
        try:
            if hash_comp["problem_structure"]:
                problem = True
        except KeyError:
            pass

        if problem:
            pprint("n" * len(hash_comp["pri_seq"]))
        else:
            if hash_comp["mature_beg"] < hash_comp["loop_beg"]:

                bef = "f" * (hash_comp["mature_beg"] - 1)
                pprint(bef)
                bef = "M" * (hash_comp["mature_end"] -
                             hash_comp["mature_beg"] + 1)
                pprint(bef)
                bef = "l" * (hash_comp["star_beg_obs"] -
                             hash_comp["mature_end"] - 1)
                pprint(bef)
                bef = "S" * (hash_comp["star_end_obs"] -
                             hash_comp["star_beg_obs"] + 1)
                pprint(bef)
                bef = "f" * (hash_comp["pri_end"] - hash_comp["star_end_obs"])
                pprint(bef)
            else:
                bef = "f" * (hash_comp["star_beg_obs"] - 1)
                pprint(bef)
                bef = "S" * (hash_comp["star_end_obs"] -
                             hash_comp["star_beg_obs"] + 1)
                pprint(bef)
                bef = "l" * (hash_comp["mature_beg"] -
                             hash_comp["star_end_obs"] - 1)
                pprint(bef)
                bef = "M" * (hash_comp["mature_end"] -
                             hash_comp["mature_beg"] + 1)
                pprint(bef)
                bef = "f" * (hash_comp["pri_end"] - hash_comp["mature_end"])
                pprint(bef)

    pprint("\npri_seq{}{}\n".format(spacer, hash_comp['pri_seq'].lower()))
    spacer = " " * (col1_width - len("pri_struct"))
    pprint("pri_struct{}{}\t#MM\n".format(spacer, hash_comp['pri_struct']))
    reads = re.split(r'\n', lines)
    lr = len(reads)
    rrc = 0
    for r in reads:
        m = re.match(
            r'^(\S+)\s+\d+\s+\d+\s+\d+\s+(\S+)\s+\S+\s+\d+\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\d+).+$', r)
        if m:
            m = m.groups()
            rrc += 1
            create_hash_key_chain(reads_hash, m[0], rrc, 'id')
            create_hash_key_chain(reads_hash, m[1], rrc, 'seq')
            create_hash_key_chain(reads_hash, int(m[2]), rrc, 'beg')
            create_hash_key_chain(reads_hash, int(m[3]), rrc, 'end')
            create_hash_key_chain(reads_hash, m[4], rrc, 'mm')
            reads_hash[rrc]["id"] = m[0]
            reads_hash[rrc]["seq"] = m[1]
            reads_hash[rrc]["beg"] = int(m[2])
            reads_hash[rrc]["end"] = int(m[3])
            reads_hash[rrc]["mm"] = m[4]

    # sorted keys by begin postion
    skeys = hash_sort_key(reads_hash, lambda x: reads_hash[x[0]]['beg'])
    elist = []  # final sorted array

    # all keys that have same begin position should match this value
    first = reads_hash[skeys[0]]["beg"]
    rorder = {}  # temporary hash to store all keys with same begin position

    for j in range(0, len(skeys)):
        if reads_hash[skeys[j]]["beg"] == first:
            # insert key and end position to hash
            rorder[j] = reads_hash[skeys[j]]["end"]
        else:   # if new begin position
            first = reads_hash[skeys[j]]["beg"]
            # sort hash keys by end position
            for k in hash_sort_key(rorder, lambda x: rorder[x[0]]):
                elist.append(skeys[k])
            for k in rorder.keys():
                del rorder[k]
            rorder[j] = reads_hash[skeys[j]]["end"]

    for k in hash_sort_key(rorder, lambda x: rorder[x[0]]):
        elist.append(skeys[k])

    for l in elist:
        rseq = reads_hash[l]['seq'].lower()
        rseq = tr(rseq, 't', 'u')
        bef = "." * (reads_hash[l]['beg'] - 1)
        after = "." * (hash_comp['pri_end'] - reads_hash[l]["end"])
        spacer = " " * (col1_width - len(reads_hash[l]['id']))
        sread = ssplit(rseq)

        bshift = 0
        rseq = ""
        for i in range(0, len(sread)):
            if not pseq[i + reads_hash[l]['beg'] - 1]:    # read is longer than sequence
                pass
            else:
                if pseq[i + reads_hash[l]['beg'] - 1] == sread[i]:
                    rseq += sread[i].lower()
                else:
                    sread[i] = sread[i].upper()
                    rseq += sread[i].upper()

        pprint("{}{}{}{}{}\t{}\n".format(
            reads_hash[l]['id'],
            spacer,
            bef,
            rseq,
            after,
            reads_hash[l]['mm']
        ))

    pprint("\n\n\n")


def print_results(ret, options):
    '''
    print out if the precursor is accepted and accepted precursors should be printed out
    or if the potential precursor is discarded and discarded potential precursors should
    be printed out
    '''
    # print_stderr(ret, options)
    if ((not options.get('-t') == '' and ret) or (options.get('-t') == '' and not ret)):
        # full output
        if not(options.get('-u') == ''):
            print_hash_comp(options)
            # print $lines."\n\n";
            return

        # limited output (only ids)
        _id = hash_comp["pri_id"]
        pprint("{}\n".format(_id))


def reset_variables():
    global hash_query
    global hash_comp
    global hash_bp
    global hash_stars
    global lines

    hash_query = {}
    hash_comp = {}
    hash_bp = {}
    hash_stars = {}
    lines = ''


def test_randfold():
    global hash_comp, hash_randfold
    pri_id = hash_comp["pri_id"]
    p_value = hash_randfold[pri_id]
    if p_value <= 0.05:
        return 1
    return 0


def cdf_gumbel(var, scale, location):
    '''
    calculates the cumulative distribution function of the Gumbel distribution
    '''
    global e

    cdf = e**(-(e**(-(var - location) / scale)))

    return cdf


def prob_gumbel_discretized(var, scale, location):
    '''
    discretized Gumbel distribution, probabilities within windows of 1 kCal/mol
    uses the subroutine that calculates the cdf to find the probabilities
    '''

    bound_lower = var - 0.5
    bound_upper = var + 0.5
    cdf_lower = cdf_gumbel(bound_lower, scale, location)
    cdf_upper = cdf_gumbel(bound_upper, scale, location)
    prob = cdf_upper - cdf_lower
    return prob


def score_mfe(mfe):
    '''
    scores the minimum free energy in kCal/mol of the potential precursor
    Assumes Gumbel distribution as described in methods section of manuscript
    '''
    # numerical value, minimum 1
    mfe_adj = max2(1, -mfe)

    # parameters of known precursors and background hairpins, scale and
    # location
    prob_test = prob_gumbel_discretized(mfe_adj, 5.5, 32)
    prob_background = prob_gumbel_discretized(mfe_adj, 4.8, 23)

    odds = prob_test / prob_background
    log_odds = math.log(odds)

    return log_odds


def score_freq(freq, options):
    '''
    scores the count of reads that map to the potential precursor
    Assumes geometric distribution as described in methods section of manuscript
    '''
    # parameters of known precursors and background hairpins
    parameter_test = 0.999
    parameter_control = 0.6

    # log_odds calculated directly to avoid underflow
    intercept = math.log((1 - parameter_test) / (1 - parameter_control))
    slope = math.log(parameter_test / parameter_control)
    log_odds = slope * freq + intercept

    # if no strong evidence for 3' overhangs, limit the score contribution to 0
    try:
        if not(options.get('-x') == '' or hash_comp["star_read"]):
            log_odds = min2(log_odds, 0)
    except KeyError:
        log_odds = min2(log_odds, 0)

    return log_odds


def test_seed_conservation():
    global hash_comp, seed_lng, hash_seeds, hash_mirs
    seed = substr(hash_comp["mature_seq"], 1, seed_lng)
    seed = tr(seed, '[T]', '[U]')
    if hash_value_true(hash_seeds, seed):
        seed_family = hash_mirs[seed]
        m = re.match(r'^(\S+)', seed_family)
        if m:
            # example of a miRNA with the same seed
            hash_comp["seed_family"] = m.groups()[0]

        return 1

    return 0


def pass_threshold_score(options):
    '''
    this is the scoring
    '''
    global seed_lng, hash_comp, score_seed, score_seed_not, scores_stem, score_min, score_star, score_star_not, score_randfold, score_randfold_not
    # minimum free energy of the potential precursor
    _score_mfe = score_mfe(hash_comp["pri_mfe"])

    # count of reads that map in accordance with Dicer processing
    _score_freq = score_freq(hash_comp["freq_total"], options)

    # basic score
    score = _score_mfe + _score_freq

    # scoring of conserved seed/seed (optional)
    if options.get('-s'):

        # if the seed is conserved
        if(test_seed_conservation()):

            # seed from position 2-8
            seed = substr(hash_comp["mature_seq"], 1, seed_lng)

            # resolve DNA/RNA ambiguities
            seed = tr(seed, '[T]', '[U]')

            # print score contribution
            hash_comp["score_seed"] = score_seed

            # add to score
            score += score_seed

        # if the seed is not conserved
        else:
            # print (negative) score contribution
            hash_comp["score_seed"] = score_seed_not

            # add (negative) score contribution
            score += score_seed_not

    # if the majority of potential star reads fall as expected from Dicer
    # processing
    try:
        if hash_comp['star_read']:
            hash_comp["score_star"] = score_star
            score += score_star
        else:
            hash_comp["score_star"] = score_star_not
            score += score_star_not
    except KeyError:
        hash_comp["score_star"] = score_star_not
        score += score_star_not

    # score lower stems for potential for Drosha recognition (highly optional)
    if options.get('-z') == '':
        stem_bp = hash_comp["stem_bp"]
        score_stem = scores_stem[stem_bp]
        score += score_stem
        hash_comp["score_stem"] = score_stem

    # score for randfold (optional)
    if options.get('-y'):
        # randfold p-value not defined for this potential precursor
        defined = False
        try:
            hash_randfold[hash_comp['pri_id']]
            defined = True
        except KeyError:
            pass

        if not defined:
            hash_comp["score_randfold"] = "0"
        elif test_randfold():
            # randfold value<0.05
            score += score_randfold
            hash_comp["score_randfold"] = score_randfold
        else:
            # randfold value>0.05
            score += score_randfold_not
            hash_comp["score_randfold"] = score_randfold_not

    # round off values to one decimal
    round_mfe = pround(_score_mfe * 10) / float(10)
    round_freq = pround(_score_freq * 10) / float(10)
    _round = pround(score * 10) / float(10)

    # print scores
    hash_comp["score_mfe"] = round_mfe
    hash_comp["score_freq"] = round_freq
    hash_comp["score"] = _round
    # print_stderr(_round, '\n')

    # return 1 if the potential precursor is accepted, return 0 if discarded
    if not(score >= score_min):
        # print(score, score_min, score >= score_min)
        return 0

    # print(1)
    return 1


def resolve_potential_precursor(options):
    '''
    #    dissects the potential precursor in parts by filling hashes, and tests if it passes the
    #    initial filter and the scoring filter
    #    binary variable whether the potential precursor is still viable
    '''

    ret = 1

    fill_structure()

    fill_pri()

    fill_mature()

    fill_star()

    fill_loop()

    fill_lower_flanks(options)

    #    do_test_assemble();
    #    this is the actual classification
    if not(pass_filtering_initial() and pass_threshold_score(options)):
        ret = 0

    # if db_old == 'CHROMOSOME_I_705':
    #    print_stderr(hash_comp, "\n")

    print_results(ret, options)

    reset_variables()


def find_freq(query):
    m = re.search(r'_x(\d+)', query)
    if m:
        m = m.groups()
        freq = int(m[0])
        return freq
    else:
        return 0


def parse_file_arf(file_arf, options):
    '''
    read through the signature blastparsed file, fills up a hash with information on queries
    (deep sequences) mapping to the current db (potential precursor) and resolve each
    potential precursor in turn
    '''
    global db_old, hash_query, lines
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
            edits = m[11]
            edit_string = m[12]

            # only reads that map sense to the potential precursor are
            # considered
            if strand == "-":
                continue

            # if the new line concerns a new db (potential precursor) then the
            # old db must be resolved
            if db_old and db_old != db:
                resolve_potential_precursor(options)

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
            lines += line

    resolve_potential_precursor(options)


def test_flanks(beg, end):
    '''
    tests the positions of the lower flanks
    '''
    global hash_comp

    # unless the begin and end positions are plausible, test negative
    if not(beg > 0 and beg <= hash_comp["pri_end"] and end > 0 and end <= hash_comp["pri_end"] and beg <= end):
        return 0

    return 1


def fill_stems_drosha():
    '''
    scores the number of base pairings formed by the first ten nt of the lower stems
    in general, the more stems, the higher the score contribution
    warning: this options has not been thoroughly tested
    '''
    global hash_comp

    flank_first_struct = hash_comp["flank_first_struct"]
    flank_second_struct = hash_comp["flank_second_struct"]
    stem_first = substr(flank_first_struct, -10)
    stem_second = substr(flank_second_struct, 0, 10)
    stem_bp_first = 0
    stem_bp_second = 0

    # find base pairings by simple pattern matching
    while re.search(r'\(', stem_first):
        stem_bp_first += 1

    while re.search(r'\)', stem_second):
        stem_bp_second += 1

    stem_bp = min2(stem_bp_first, stem_bp_second)

    hash_comp["stem_first"] = stem_first
    hash_comp["stem_second"] = stem_second
    hash_comp["stem_bp_first"] = stem_bp_first
    hash_comp["stem_bp_second"] = stem_bp_second
    hash_comp["stem_bp"] = stem_bp


def fill_lower_flanks(options):
    '''
    fills specifics on the lower flanks and unpaired strands into the 'comp' hash
    unless both mature and star sequences are plausible, do not look for the flanks
    '''
    global hash_comp

    try:
        if not(hash_comp['mature_arm'] and hash_comp['star_arm']):
            return
    except KeyError:
        return

    flank_first_end = None
    flank_second_beg = None

    # defining the begin and end positions of the flanks from the mature and star positions
    # excision depends on whether the mature or star sequence is 5' in the
    # potenitial precursor ('first')
    if hash_comp["mature_arm"] == "first":
        flank_first_end = hash_comp["mature_beg"] - 1
    else:
        flank_second_beg = hash_comp["mature_end"] + 1

    if hash_comp["star_arm"] == "first":
        flank_first_end = hash_comp["star_beg_exp"] - 1
    else:
        flank_second_beg = hash_comp["star_end_exp"] + 1

    # unless the positions are plausible, do not fill into hash
    if not(test_flanks(flank_first_end, flank_second_beg)):
        return

    hash_comp["flank_first_end"] = flank_first_end
    hash_comp["flank_second_beg"] = flank_second_beg
    hash_comp["flank_first_seq"] = excise_seq(hash_comp["pri_seq"], hash_comp[
                                              "pri_beg"], hash_comp["flank_first_end"], "+")
    hash_comp["flank_second_seq"] = excise_seq(hash_comp["pri_seq"], hash_comp[
                                               "flank_second_beg"], hash_comp["pri_end"], "+")
    hash_comp["flank_first_struct"] = excise_struct(
        hash_comp["pri_struct"], hash_comp["pri_beg"], hash_comp["flank_first_end"], "+")
    hash_comp["flank_second_struct"] = excise_struct(
        hash_comp["pri_struct"], hash_comp["flank_second_beg"], hash_comp["pri_end"], "+")

    if options.get('-z') == '':
        fill_stems_drosha()


def create_hash_seeds(infile):
    '''
    parses a fasta file with sequences of known miRNAs considered for conservation purposes
    reads the seeds into a hash
    '''
    global seed_lng, hash_seeds, hash_mirs
    _id, desc, sequence, seed = None, None, None, None
    FASTA = open_or_die(infile, 'rb', 'can not open {}\n'.format(infile))
    while True:
        line = FASTA.readline()
        if not line:
            break

        line = line.strip()
        m = re.match(r'^>(\S+)(.*)', line)
        if m:
            m = m.groups()
            _id = m[0]
            desc = m[1]
            sequence = ""
            seed = ""
            while True:
                line2 = FASTA.readline()
                if not line2:
                    break

                line2 = line2.strip()
                mm = re.match(r'^>(\S+)(.*)', line2)
                if mm:
                    mm = mm.groups()
                    seed = substr(sequence, 1, seed_lng)
                    seed = tr(seed, '[T]', '[U]')
                    create_hash_key_chain(hash_mirs, '', seed)
                    hash_mirs[seed] += "{}\t".format(_id)
                    create_hash_key_chain(hash_seeds, 0, seed)
                    hash_seeds[seed] += 1

                    _id = mm[0]
                    desc = mm[1]
                    sequence = ""
                    seed = ""
                    continue

                sequence += line2

    seed = substr(sequence, 1, seed_lng)
    seed = tr(seed, '[T]', '[U]')
    create_hash_key_chain(hash_mirs, '', seed)
    hash_mirs[seed] += "{}\t".format(_id)
    create_hash_key_chain(hash_seeds, 0, seed)
    hash_seeds[seed] += 1
    FASTA.close()


def parse_file_randfold(infile):
    global hash_randfold

    FILE = open_or_die(infile, 'rb', 'can not open {}\n'.format(infile))

    while True:
        line = FILE.readline()
        if not line:
            break
        m = re.match(r'^(\S+)\s+(\S+)\s+(\S+)', line)
        if m:
            m = m.groups()
            _id = m[0]
            randfold = m[2]
            hash_randfold[_id] = float(randfold)

    FILE.close()


def parse_file_struct(file_structure):
    global hash_desc, hash_seq, hash_struct, hash_mfe

    FILE_STRUCT = open_or_die(file_structure, 'rb',
                              'can not open {}\n'.format(file_structure))
    _id, desc, seq, struct, mfe = None, None, None, None, None
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
                    mm = mm.groups()
                    hash_desc[_id] = desc
                    hash_seq[_id] = seq
                    hash_struct[_id] = struct
                    hash_mfe[_id] = mfe

                    _id = mm[0]
                    desc = mm[1]
                    seq = ""
                    struct = ""
                    mfe = ""

                    continue

                if re.match(r'^\w', line2):
                    line2 = tr(line2, 'uU', 'tT')
                    seq += line2

                mm = re.search(r'((\.|\(|\))+)', line2)
                if mm:
                    mm = mm.groups()
                    struct += mm[0]

                mm = re.search(r'\((\s*-\d+\.\d+)\)', line2)
                if mm:
                    mm = mm.groups()
                    mfe = float(mm[0])

    hash_desc[_id] = desc
    hash_seq[_id] = seq
    hash_struct[_id] = struct
    hash_mfe[_id] = mfe

    # print(hash_mfe)
    FILE_STRUCT.close()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        die(usage)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_signature', help='signature file')
    parser.add_argument('file_structure', help='structure file')
    args = parser.parse_args(sys.argv[1:3])

    file_arf = args.file_signature
    file_struct = args.file_structure

    opts, argss = getopt.getopt(sys.argv[3:], 'hs:tuv:xy:zl:')
    options = dict(opts)

    if options.get('-v'):
        score_min = float(options.get('-v'))

    if options.get('-x') == '':
        score_min = -5

    # print help if that option is used
    if options.get('-h') == '':
        die(usage)

    # parse structure file outputted from RNAfold
    parse_file_struct(file_struct)

    # if conservation is scored, the fasta file of known miRNA sequences is
    # parsed
    if options.get('-s'):
        create_hash_seeds(options.get('-s'))

    # if randfold p-value is scored, a randfold output file of the precursors
    # is parsed
    if options.get('-y'):
        parse_file_randfold(options.get('-y'))

    # parse signature file in arf format and resolve each potential precursor
    parse_file_arf(file_arf, options)
