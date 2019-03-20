#!/usr/bin/env python
from __future__ import print_function

import argparse
import getopt
import os
import re
import sys

from port import (chop, create_hash_key_chain, die, open_or_die, pprint,
                  print_stderr)

usage = '''
{} file_arf

Parses an arf file, and discards select mappings. Any files inputted with
the options should be single-column format, with a single id on each line.

-a int     Discard mappings of edit distance higher than this
-b int     Discard mappings of read queries shorter than this
-c int     Discard mappings of read queries longer than this
-d file    Discard read queries not in this file
-e file    Discard read queries in this file
-f file    Discard reference dbs not in this file
-g file    Discard reference dbs in this file
-h         Discard remaining suboptimal mappings
-i int     Discard remaining suboptimal mappings and discard any
           reads that have more remaining mappings than this
-j         Remove any unmatched nts in the very 3' end
-k         Output progress to standard output

'''.format(sys.argv[0])

# global hashes
hash_fasta = {}
hash_edits = {}
hash_edit_best = {}
hash_nr_mappings = {}

hash_queries_incl = {}
hash_queries_excl = {}
hash_dbs_incl = {}
hash_dbs_excl = {}

# running index
running = 0

# scan mode - a flag variable
gscan = 0


def parse_file_ids(_file, _hash):
    # read id file into hash
    if options.get('-k') == '':
        print_stderr('reading id file into memory\n')

    FILE = open_or_die(_file, 'rb', 'can not open {}\n'.format(_file))
    while True:
        line = FILE.readline()
        if not line:
            break

        m = re.match(r'^(\S+)', line)
        if m:
            _id = m.groups()[0]
            _hash[_id] = 1


def remove_trailing_nts(query_map_lng, query_end, query_seq, db_map_lng, db_end, db_seq, edits, edit_string):
    while re.search(r'M$', edit_string):
        query_map_lng -= 1
        query_end -= 1
        query_seq = chop(query_seq)
        db_map_lng -= 1
        db_end -= 1
        db_seq = chop(db_seq)
        edits -= 1
        edit_string = chop(edit_string)

    return (query_map_lng, query_end, query_seq, db_map_lng, db_end, db_seq, edits, edit_string)


def evaluate_query(query, edits, options):
    global hash_edit_best, hash_nr_mappings
    if hash_edit_best[query] < edits:
        return 0

    if options.get('-i') and options.get('-i') < hash_nr_mappings[query]:
        return 0

    return 1


def parse(options, file_arf):
    if options.get('-k') == '':
        print_stderr('parsing and printing mappings\n')

    parse_file_arf(file_arf, options)

    if options.get('-k') == '':
        print_stderr('mappings printed\n\n')


def parse_file_arf(file_arf, options):
    global running, gscan, hash_edits
    FILE_ARF = open_or_die(
        file_arf, 'rb', 'can not open {}\n'.format(file_arf))
    while True:
        line = FILE_ARF.readline()
        if not line:
            break

        m = re.match(
            r'^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if m:
            m = m.groups()
            query = m[0]
            query_map_lng = int(m[1])
            query_beg = m[2]
            query_end = int(m[3])
            query_seq = m[4]
            db = m[5]
            db_map_lng = int(m[6])
            db_beg = m[7]
            db_end = int(m[8])
            db_seq = m[9]
            strand = m[10]
            edits = int(m[11])
            edit_string = m[12]

            running += 1
            if options.get('-j') == '':
                (query_map_lng, query_end, query_seq, db_map_lng, db_end, db_seq, edits, edit_string) = remove_trailing_nts(
                    query_map_lng, query_end, query_seq, db_map_lng, db_end, db_seq, edits, edit_string)

            if '-a' in options.keys() and int(options.get('-a')) < edits:
                continue

            if options.get('-b') and query_map_lng < int(options.get('-b')):
                continue

            if options.get('-c') and int(options.get('-c')) < query_map_lng:
                continue

            if options.get('-d') and query not in hash_queries_incl.keys():
                continue

            if options.get('-e') and query in hash_queries_excl.keys():
                continue

            if options.get('-f') and db not in hash_dbs_incl.keys():
                continue

            if options.get('-g') and db in hash_dbs_excl.keys():
                continue

            if not(options.get('-h') == '' or options.get('-i')):
                pprint('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    query,
                    query_map_lng,
                    query_beg,
                    query_end,
                    query_seq,
                    db,
                    db_map_lng,
                    db_beg,
                    db_end,
                    db_seq,
                    strand,
                    edits,
                    edit_string
                ))
                continue

            if gscan:
                create_hash_key_chain(hash_edits, 0, query, edits)
                hash_edits[query][edits] += 1
            else:
                evaluation = evaluate_query(query, edits, options)
                if evaluation:
                    pprint("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        query,
                        query_map_lng,
                        query_beg,
                        query_end,
                        query_seq,
                        db,
                        db_map_lng,
                        db_beg,
                        db_end,
                        db_seq,
                        strand,
                        edits,
                        edit_string
                    ))


def fill_hash():
    global hash_edits, hash_edit_best, hash_nr_mappings
    queries = sorted(hash_edits.keys())
    for query in queries:

        # find the best mapping (lowest edit distance) for the query read
        edits = sorted(hash_edits[query].keys())
        edit_best = edits[0]

        # find the number of mappings of this edit distance
        hash_edit_best[query] = edit_best
        hash_nr_mappings[query] = hash_edits[query][edit_best]


def scan(file_arf, options):
    global gscan, running
    if options.get('-k') == '':
        lines = os.popen('cat {} | wc -l'.format(file_arf)).read().strip()
        print_stderr('scanning mappings, total={}\n'.format(lines))

    gscan = 1
    parse_file_arf(file_arf, options)

    gscan = 0
    if options.get('-k') == '':
        print_stderr('resolving best mappings for each read\n')
    fill_hash()

    running = 0


if __name__ == '__main__':
    if len(sys.argv) < 2:
        die(usage)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('file_arf', help='arf file')
    args = parser.parse_args(sys.argv[1:2])
    file_arf = args.file_arf

    opts, argss = getopt.getopt(sys.argv[2:], "a:b:c:d:e:f:g:hi:jk")
    options = dict(opts)

    if options.get('-d'):
        parse_file_ids(options.get('-d'), hash_queries_incl)

    if options.get('-e'):
        parse_file_ids(options.get('-e'), hash_queries_excl)

    if options.get('-f'):
        parse_file_ids(options.get('-f'), hash_dbs_incl)

    if options.get('-g'):
        parse_file_ids(options.get('-g'), hash_dbs_excl)

    if options.get('-h') == '' or options.get('-i'):
        scan(file_arf, options)

    parse(options, file_arf)
