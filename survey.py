#!/usr/bin/env python
from __future__ import print_function

import argparse
import getopt
import math
import re
import sys

from port import (create_hash_key_chain, die, div, hash_defined, max2, mean,
                  min2, open_or_die, pprint, print_stderr, pround,
                  round_decimals, sd)

usage = '''
Usage:

    {}  output_mor

    Options:

    -a file  file outputted by controls
    -b file  mature miRNA fasta reference file for the species
    -c file  signature file
    -d int   read stack height necessary for triggering excision

'''.format(sys.argv[0])

hash_con = {}
hash_ref = {}
hash_sig = {}
hash_out = {}


def parse_file_ref(filename, hash_ref):
    _id, desc, seq = None, None, None

    FASTA = open_or_die(filename, 'rb', 'can not open {}\n'.format(filename))
    while True:
        l = FASTA.readline()
        if not l:
            break

        l = l.strip()

        m = re.match(r'^>(\S+)(.*)', l)
        if m:
            m = m.groups()
            _id = m[0]
            desc = m[1]
            seq = ''
            while True:
                ll = FASTA.readline()
                if not ll:
                    break

                ll = ll.strip()

                mm = re.match(r'^>(\S+)(.*)', ll)
                if mm:
                    mm = mm.groups()
                    hash_ref[_id] = seq
                    _id = mm[0]
                    desc = mm[1]
                    seq = ''
                    continue
                seq += ll

    hash_ref[_id] = seq
    FASTA.close()


def parse_file_sig(file_sig, hash_sig):
    global hash_ref
    FILE = open_or_die(file_sig, 'rb', 'can not open {}\n'.format(file_sig))
    while True:
        line = FILE.readline()
        if not line:
            break

        m = re.match(r'^(\S+)', line)

        # Test isdefined
        defined = False
        if m:
            try:
                hash_ref[m.groups()[0]]
                defined = True
            except KeyError:
                pass

        if m and defined:
            # dummy value
            hash_sig[m.groups()[0]] = -50


def resolve_entry_file_mrd_permuted(score, permutation, refs, _hash, options):
    if permutation is None:
        print_stderr('The {} file is not properly formatted.\nMaybe it does not contain the lines with \"permutation int\"?\n'.format(
            options.get('-a')
        ))
        sys.exit(0)

    floor = int(math.floor(score))
    create_hash_key_chain(_hash, 0, 'total', permutation, floor)
    _hash['total'][permutation][floor] += 1

    if refs:
        create_hash_key_chain(_hash, 0, 'known', permutation, floor)
        _hash['known'][permutation][floor] += 1
    else:
        create_hash_key_chain(_hash, 0, 'novel', permutation, floor)
        _hash['novel'][permutation][floor] += 1

    for i in range(len(refs)):
        refs.pop()


def parse_file_mrd_permuted(infile, _hash, options):
    global hash_sig
    score, permutation, refs = None, None, []
    FILE = open_or_die(infile, 'rb', 'can not open {}\n'.format(infile))
    while True:
        line = FILE.readline()
        if not line:
            break

        m = re.match(r'^permutation\s+(\d+)', line)
        if m:
            permutation = m.groups()[0]
        else:
            mm = re.match(r'^score\s+total\s+(\S+)', line)
            if mm:
                score = float(mm.groups()[0])
            else:
                m3 = re.match(r'^(\S+)', line)

                defined = False
                if m3:
                    try:
                        hash_sig[m3.groups()[0]]
                        defined = True
                    except KeyError:
                        pass

                if m3 and defined:
                    refs.append(m3.groups()[0])
                elif re.search(r'^>', line) and score is not None:
                    resolve_entry_file_mrd_permuted(
                        score, permutation, refs, _hash, options)
                    # ATTR: ($$i, @$refs) = () in
                    # resolve_entry_file_mrd_permuted
                    score = None

    resolve_entry_file_mrd_permuted(score, permutation, refs, _hash, options)
    score = None


def resolve_entry_file_mrd(score, refs, _hash):
    global hash_sig
    floor = int(math.floor(score))

    create_hash_key_chain(_hash, 0, 'total', floor)
    _hash['total'][floor] += 1

    if refs:
        create_hash_key_chain(_hash, 0, 'known', floor)
        _hash['known'][floor] += 1

        for ref in refs:
            if hash_sig[ref] < floor:
                hash_sig[ref] = floor
    else:
        create_hash_key_chain(_hash, 0, 'novel', floor)
        _hash['novel'][floor] += 1

    for i in range(len(refs)):
        refs.pop()  # clear refs array


def parse_file_mrd(file_out, _hash):
    global hash_sig
    score = None
    refs = []
    FILE = open_or_die(file_out, 'rb', 'can not open {}\n'.format(file_out))
    while True:
        line = FILE.readline()
        if not line:
            break

        m = re.match(r'^score\s+total\s+(\S+)', line)
        if m:
            m = m.groups()
            score = float(m[0])
        else:
            m2 = re.match(r'^(\S+)', line)

            defined = False
            if m2:
                try:
                    hash_sig[m2.groups()[0]]
                    defined = True
                except KeyError:
                    pass

            if m2 and defined:
                refs.append(m2.groups()[0])
            elif re.search(r'^>', line) and score is not None:
                resolve_entry_file_mrd(score, refs, _hash)
                score = None  # ATTR: ($$score,@$refs)=();

    resolve_entry_file_mrd(score, refs, _hash)


def print_header(options):
    pprint("moRNA Finder score")

    if options.get('-b'):
        pprint("\tnovel miRNAs reported by moRNA Finder")

    if options.get('-a'):
        pprint("\tnovel miRNAs, estimated false positives")

        pprint("\tnovel miRNAs, estimated true positives")

    pprint("\tknown miRNAs in species")

    pprint("\tknown miRNAs in data")

    pprint("\tknown miRNAs detected by moRNA Finder")

    if options.get('-a'):
        pprint("\testimated signal-to-noise")

    if options.get('-d'):
        pprint("\texcision gearing")

    pprint("\n")


def hairpins_cnt(htype, score):
    global hash_out
    hairpins_cnt = 0

    scores_cnt = sorted(hash_out[htype].keys())
    for score_cnt in scores_cnt:
        if score <= score_cnt:
            count = hash_out[htype][score_cnt]
            hairpins_cnt += count

    return hairpins_cnt


def mean_sd(htype, score, norm):
    '''
    calculate all summary statistics (estimated number of false and true positives, standard deviations, percentages of total hairpins that are true)
    these statistics are calculated from the permuted controls
    '''
    global hash_con
    counts = []
    counts_above_norm = []
    percentages = []

    permutations = sorted(hash_con[htype].keys())

    # iterate over runs of permuted controls. Each run is separated by
    # 'permutation int' in the output_permuted.mrd file.
    for permutation in permutations:
        cnt_total = 0

        scores_cnt = sorted(hash_con[htype][permutation].keys())
        for score_cnt in scores_cnt:
            if score <= score_cnt:
                count = hash_con[htype][permutation][score_cnt]
                cnt_total += count

        counts.append(cnt_total)

    fp_mean = mean(counts)
    fp_sd = sd(counts)

    # iterate over runs of permuted controls again, using the @counts array
    # just generated
    for count in counts:

        # estimated number of true positives t for this run
        # t = number of hairpins reported - estimated number of false positives
        # for this run
        cnt_total_above_norm = max2(0, norm - count)
        percentage = max2(0, 100 * (div(cnt_total_above_norm, norm)))

        counts_above_norm.append(cnt_total_above_norm)
        percentages.append(percentage)

    # calculate mean, sd, and percentages
    # error
    true_mean = mean(counts_above_norm)
    true_sd = sd(counts_above_norm)
    percentage_mean = mean(percentages)
    percentage_sd = sd(percentages)

    return (fp_mean, fp_sd, true_mean, true_sd, percentage_mean, percentage_sd)


def survey_hairpins(score):
    '''
    summary statistics for miRNA hairpin precursors
    '''

    # partitioning of hairpins into known and novel
    hairpins_known = hairpins_cnt("known", score)
    hairpins_novel = hairpins_cnt("novel", score)

    # print_stderr(hairpins_novel, "\t")

    pprint('\t{}'.format(hairpins_novel))

    # estimation of false positives for the set of novel hairpins
    if options.get('-a'):
        (hairpins_fp_mean, hairpins_fp_sd, estimated_true_mean,
         estimated_true_sd, percent_mean, percent_sd) = (0, 0, 0, 0, 0, 0)

        if hairpins_novel:  # check if novel hairpins detected
            (hairpins_fp_mean, hairpins_fp_sd, estimated_true_mean, estimated_true_sd,
             percent_mean, percent_sd) = mean_sd("novel", score, hairpins_novel)

        hairpins_fp_mean_round = round_decimals(hairpins_fp_mean, 0)
        hairpins_fp_sd_round = round_decimals(hairpins_fp_sd, 0)
        pprint("\t{} +/- {}".format(hairpins_fp_mean_round, hairpins_fp_sd_round))

        estimated_true_mean_round = round_decimals(estimated_true_mean, 0)
        pprint("\t{}".format(estimated_true_mean_round))

        estimated_true_sd_round = round_decimals(estimated_true_sd, 0)
        pprint(" +/- {}".format(estimated_true_sd_round))

        percent_mean_round = round_decimals(percent_mean, 0)
        pprint(" ({}".format(percent_mean_round))

        percent_sd_round = round_decimals(percent_sd, 0)
        pprint(" +/- {}%)".format(percent_sd_round))


def survey_known(score):
    '''
    summary statistics for mature miRNAs
    '''
    global hash_ref, hash_sig

    # mature miRNAs for the species in reference (miRBase) file
    matures_cnt = len(hash_ref)
    pprint('\t{}'.format(matures_cnt))

    # matures present in data
    _matures_present = len(hash_sig.keys())
    pprint("\t{}".format(_matures_present))

    # matures recovered
    matures_recov_cnt = 0
    matures_present = hash_sig.keys()
    for mature_present in matures_present:
        score_mature_present = hash_sig[mature_present]
        if score <= score_mature_present:
            matures_recov_cnt += 1
    pprint("\t{}".format(matures_recov_cnt))

    # matures recovered in percent
    percent = 100 * round_decimals(div(matures_recov_cnt, _matures_present), 2)
    pprint(" ({:.0f}%)".format(percent))


def survey_signal_to_noise(score):
    # total hairpins reported
    hairpins_total = hairpins_cnt("total", score)

    (hairpins_total_fp_mean, hairpins_total_fp_sd, estimated_total_true_mean, estimated_total_true_sd,
     percent_total_mean, percent_total_sd) = mean_sd("total", score, hairpins_total)

    hairpins_total_fp_mean_round = round_decimals(hairpins_total_fp_mean, 0)

    hairpins_total_fp_sd_round = round_decimals(hairpins_total_fp_sd, 0)

    # print "\t$hairpins_total_fp_mean_round +/- $hairpins_total_fp_sd_round";

    signal_to_noise_total = 0
    if hairpins_total_fp_mean == 0:
        signal_to_noise_total = 0
    else:
        signal_to_noise_total = round_decimals(
            div(hairpins_total, hairpins_total_fp_mean), 1)

    pprint("\t{}".format(signal_to_noise_total))


def survey(options):
    # print(hash_con)
    # print_stderr(hash_out, '\n')
    # print(hash_sig)
    # print(hash_out)
    for score in range(10, -11, -1):
        pprint(score)

        if options.get('-b'):
            survey_hairpins(score)

            survey_known(score)

        if options.get('-a'):
            survey_signal_to_noise(score)

        if options.get('-d'):
            read_stack_min = options.get('-d')
            pprint('\t{}'.format(read_stack_min))

        pprint('\n')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        die(usage)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_out', help='output file')
    args = parser.parse_args(sys.argv[1:2])

    file_out = args.file_out
    opts, argss = getopt.getopt(sys.argv[2:], 'a:b:c:d:')
    options = dict(opts)

    if (options.get('-b') and not options.get('-c')) or (not options.get('-b') and options.get('-c')):
        die('options -b and -c must be used in conjunction\n')

    if options.get('-b') and options.get('-c'):
        parse_file_ref(options.get('-b'), hash_ref)

        parse_file_sig(options.get('-c'), hash_sig)

    if options.get('-a'):
        parse_file_mrd_permuted(options.get('-a'), hash_con, options)

    parse_file_mrd(file_out, hash_out)

    print_header(options)

    survey(options)
