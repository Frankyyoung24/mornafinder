#!/usr/bin/env python
from __future__ import print_function

import getopt
import os
import re
import sys
import time

from port import (create_hash_key_chain, die, file_s, fileparse, hash_sort_key,
                  hash_value_true, open_or_die, open_or_die2, print_stderr,
                  ssplit, substr, tr)

usage = '''
usage:
  \tpython {} [options] -p precursor.fa -m mature.fa -r reads.fa -s star.fa -t species -y [timestamp] -d [pdfs] -o [sort] -k [stringent] -c config.txt -g [number of mismatches in reads vs precursor mappings]

[options]

[mandatory parameters]
  \t-u\tlist all values allowed for the species parameter that have an entry at UCSC

  \t-p precursor.fa  miRNA precursor sequences from miRBase
  \t-m mature.fa     miRNA sequences from miRBase
  \t-P               specify this option of your mature miRNA file contains 5p and 3p ids only
  \t-r reads.fa      your read sequences

[optional parameters]
  \t-c [file]    config.txt file with different sample ids... or just the one sample id
  \t-s [star.fa] optional star sequences from miRBase
  \t-t [species] e.g. Mouse or mmu
  \t             if not searching in a specific species all species in your files will be analyzed
  \t             else only the species in your dataset is considered
  \t-y [time]    optional otherwise its generating a new one
  \t-d           if parameter given pdfs will not be generated, otherwise pdfs will be generated
  \t-o           if parameter is given reads were not sorted by sample in pdf file, default is sorting
  \t-k           also considers precursor-mature mappings that have different ids, eg let7c
  \t             would be allowed to map to pre-let7a
  \t-n           do not do file conversion again
  \t-x           do not do mapping against precursor again
  \t-g [int]     number of allowed mismatches when mapping reads to precursors, default 1
  \t-e [int]     number of nucleotides upstream of the mature sequence to consider, default 2
  \t-f [int]     number of nucleotides downstream of the mature sequence to consider, default 5
  \t-j           do not create an output.mrd file and pdfs if specified\n
  \t-w           considers the whole precursor as the 'mature sequence'
  \t-W           read counts are weighed by their number of mappings. e.g. A read maps twice so each position gets 0.5 added to its read profile
    \t-U           use only unique read mappings; Caveat: Some miRNAs have multiple precursors. These will be underestimated in their expression since the multimappers are excluded
\n'''.format(sys.argv[0])

__DATA = '''hsa                 Human
ptr                 Chimp
na                  Orangutan
na                  Rhesus
na                  Marmoset
mmu                 Mouse
rno                 Rat
na                  Guinea Pig
lca                 Cat
cfa                 Dog
eca                 Horse
bta                 Cow
na                  Opossum
na                  Platypus
gga                 Chicken
na                  Zebra finch
na                  Lizard
xtr                 X.tropicalis
dre                 Zebrafish
tni                 Tetraodon
fru                 Fugu
na                  Stickleback
na                  Medaka
na                  Lamprey
bfl                 Lancelet
cin                 C.intestinalis
spu                 S.purpuratus
cel                 C.elegans
na                  C.brenneri
cbr                 C.briggsae
na                  C.remanei
sja                 C.japonica
na                  P.pacificus
dme                 D.melanogaster
dsi                 D.simulans
dse                 D.sechellia
dya                 D.yakuba
der                 D.erecta
dan                 D.ananassae
dps                 D.pseudoobscura
dpe                 D.persimilis
dvi                 D.virilis
dmo                 D.mojavensis
dgr                 D.grimshawi
aga                 A.gambiae
ame                 A.mellifera
na                  S.cerevisiae
cel                 worm'''

_hash = {}
hash_star = {}
hash_sample = {}
hash_star_sample = {}
total = {}
total_t = 0
mapcounts = {}
seen = {}
species = "none"

ctime = long(time.time())

organisms = {}
rorganisms = {}
u = None
v = None

name1 = None
name2 = None
name3 = None
name4 = None

upstream = 0
downstream = 0
samples = {}


def ConvertFastaFile(ifile, ofile, des, sp):
    global outdir
    IN = open_or_die(ifile, 'rb', 'File {} not found\n'.format(ifile))

    filename = ''
    if des == '':
        filename = '{}/{}.converted'.format(outdir, ofile)
    else:
        filename = '{}/{}.converted'.format(outdir, des)
    OUT = open_or_die(
        filename, 'w+', 'file {} could not be created'.format(filename))

    line, _id, tempid, seq, first = None, None, None, None, 1
    sp_hits = 0

    while True:
        line = IN.readline().strip()
        if not line:
            break
        m = re.match(r'^(>\S+)\s*(\S*)', line)
        if m:
            tempid = m.groups()[0]
            if not first:
                if sp == 'none':
                    if des == '':
                        if not re.search(r'N', seq):
                            OUT.write('{}\n{}\n'.format(_id, seq))
                            sp_hits += 1
                    else:
                        if not re.search(r'N', seq):
                            OUT.write('{}\n{}\n'.format(_id, seq))
                            sp_hits += 1
                else:
                    if re.search(sp, _id, re.IGNORECASE):
                        if des == '':
                            if not re.search(r'N', seq):
                                OUT.write('{}\n{}\n'.format(_id, seq))
                                sp_hits += 1
                        else:
                            if not re.search(r'N', seq):
                                OUT.write('{}\n{}\n'.format(_id, seq))
                                sp_hits += 1
            else:
                first = 0

            seq = ''
            _id = tempid
        else:
            line = line.upper()
            line = re.sub('U', 'T', line)
            seq += line

    if sp == 'none':
        if des == '':
            if not re.search(r'N', seq):
                OUT.write('{}\n{}\n'.format(_id, seq))
                sp_hits += 1
        else:
            if not re.search(r'N', seq):
                OUT.write('{}\n{}\n'.format(_id, seq))
            seq = ''
            sp_hits += 1
    elif re.search(sp, _id, re.IGNORECASE):
        if des == '':
            if not re.search('N', seq):
                OUT.write('{}\n{}\n'.format(_id, seq))
                sp_hits += 1
        else:
            if not re.search('N', seq):
                OUT.write('{}\n{}\n'.format(_id, seq))
                seq = ''
                sp_hits += 1

    IN.close()
    OUT.close()

    if not sp_hits:
        die('\nError: No entrys for species "{} or {}" found in file {}\nPlease make sure that the given species argument matches the species id in your file {} or say none\n\n\n'.format(
            options.get('-t'),
            species,
            ifile,
            ifile
        ))


def read_stats(f1, f2):
    hsh = {}
    count = 0
    k2 = {}
    total = 0

    IN = open_or_die(
        f1, 'rb', 'No reads file {} in fasta format given\n'.format(f1))

    while True:
        line = IN.readline()
        if not line:
            break

        m = re.match(r'^>*((\S\S\S)\S+_x(\d+))', line)
        if m:
            m = m.groups()
            try:
                if hsh[m[0]]:
                    continue
            except KeyError:
                pass

            hsh[m[0]] = 1
            count += int(m[2])
            if m[1] not in k2.keys():
                k2[m[1]] = 0

            k2[m[1]] += int(m[2])

    IN.close()

    hsh2 = {}
    count2 = 0
    k22 = {}

    print_stderr('Mapping statistics\n')
    IN = open_or_die(f2, 'rb', 'No mapping file given\n')
    while True:
        line = IN.readline()
        if not line:
            break

        m = re.match(r'^>*((\S\S\S)\S+_x(\d+))', line)
        if m:
            m = m.groups()
            try:
                if hsh2[m[0]]:
                    continue
            except KeyError:
                pass

            hsh2[m[0]] = 1
            count2 += int(m[2])

            if m[1] not in k22.keys():
                k22[m[1]] = 0
            k22[m[1]] += int(m[2])
    IN.close()

    print_stderr('\n#desc\ttotal\tmapped\tunmapped\t%mapped\t%unmapped\n')
    print_stderr('total: {}\t{}\t{}\t'.format(
        count,
        count2,
        count - count2
    ))
    print_stderr('{0:.3f}\t{1:.3f}\n'.format(
        float(count2) / count,
        1 - (float(count2) / count)
    ))

    for key in k2.keys():
        print_stderr('{}: {}\t{}\t{}\t'.format(
            key,
            k2[key],
            k22[key],
            k2[key] - k22[key]
        ))
        print_stderr('{0:.3f}\t{1:.3f}\n'.format(
            float(k22[key]) / k2[key],
            1 - (float(k22[key]) / k2[key])
        ))


def ReadinPrecursorFile():
    global _hash, hash_star, samples, hash_star_sample
    _id = None
    IN = open_or_die('precursor.converted', 'rb',
                     'Precursor file precursor.converted not found\n')
    while True:
        line = IN.readline().strip()
        if not line:
            break

        m = re.match(r'^>(\S+)', line)
        if m:
            _id = m.groups()[0]
            create_hash_key_chain(_hash, '', _id, 'seq')
            create_hash_key_chain(hash_star, '', _id, 'seq')

            for sample in samples.keys():
                create_hash_key_chain(hash_sample, '', sample, _id, 'seq')
                create_hash_key_chain(hash_star_sample, '', sample, _id, 'seq')
                create_hash_key_chain(hash_sample, 0, sample, _id, 'c')
                create_hash_key_chain(hash_sample, _hash[_id][
                                      'end'], sample, _id, 'end')
        else:
            _hash[_id]['seq'] = '{}{}'.format(_hash[_id]['seq'], line)
            hash_star[_id]['seq'] = '{}{}'.format(_hash[_id]['seq'], line)

        _hash[_id]['c'] = 0
        _hash[_id]['end'] = len(_hash[_id]['seq'])

        hash_star[_id]['c'] = 0
        hash_star[_id]['end'] = len(hash_star[_id]['seq'])

    IN.close()


def Mapping(options):
    global threads, name1, name2, mismatches, name3

    # bulid bowtie index
    print_stderr('building bowtie index\n')
    os.system('bowtie-build precursor.converted miRNA_precursor')

    # map mature sequences against precursors
    print_stderr('mapping mature sequences against index\n')

    # do not map mature if options are
    if options.get('-w') != '':
        os.system('bowtie -p {} -f -v 0 -a --best --strata --norc miRNA_precursor mature.converted {}_mapped.bwt 2>bowtie_mature.out'.format(
            threads,
            name1
        ))

    # map reads against precursors
    print_stderr('mapping read sequences against index\n')
    os.system('bowtie -p {} -f -v {} -a --best --strata --norc miRNA_precursor {}.converted {}_mapped.bwt 2>bowtie_reads.out'.format(
        threads,
        mismatches,
        name2,
        name2,
    ))

    read_stats("{}.converted".format(name2), "{}_mapped.bwt".format(name2))

    if options.get('-s'):
        print_stderr('mapping star sequences against index\n')
        os.system('bowtie -p {} -f -v 0 -a --best --strata --norc miRNA_precursor star.converted {}_mapped.bwt 2>bowtie_star.out'.format(
            threads,
            name3
        ))


def ReadinMatureMappingFile(options):
    global name1, _hash, hash_sample, hash_star_sample, hash_star, samples
    line = []
    matches = 0
    OUT = open_or_die('mature2hairpin', 'w+',
                      'cannot create file mature2hairpin\n')
    IN = open_or_die('{}_mapped.bwt'.format(
        name1), 'rb', 'Mature mapping file {}_mapped.bwt not found \n'.format(name1))
    cx = 0
    id1 = ''
    id2 = ''

    while True:
        l = IN.readline()
        if not l:
            break

        id1 = ''
        id2 = ''
        line = re.split(r'\t', l)

        id1 = line[0]  # this is the mature ID
        id2 = line[2]  # precursor id

        # remove multiple endings if ambigous just for matching with precursor
        id1 = re.sub(r'-5p', '', id1)
        id1 = re.sub(r'-3p', '', id1)

        # here is assumed that multiple precursor ids have 3 - in their id,
        # seems to be ok so far
        m = re.match(r'^(\w+\-\w+\-\w+)\-\d+$', id2)
        if m:
            id2 = m.groups()[0]

        if options.get('-k') != '' and not re.search(id2, id1, re.IGNORECASE) and not re.search(id1, id2, re.IGNORECASE):
            continue

        cx += 1
        # ATTR: read before assign
        create_hash_key_chain(_hash, 0, line[2], 'c')
        _hash[line[2]]['c'] += 1

        for sample in samples.keys():
            # ATTR: read before assign
            create_hash_key_chain(hash_sample, 0, sample, line[2], 'c')
            hash_sample[sample][line[2]]['c'] += 1

        matches = _hash[line[2]]['c']

        # there is a problem, Hash id is from precursor sequence, mature 7a and
        # 7b map to same precursor
        create_hash_key_chain(_hash, None, line[2], matches, 'beg')
        _hash[line[2]][matches]['beg'] = int(line[3]) - upstream
        if _hash[line[2]][matches]['beg'] < 0:
            _hash[line[2]][matches]['beg'] = 0

        create_hash_key_chain(_hash, None, line[2], matches, 'end')
        _hash[line[2]][matches]['end'] = int(
            line[3]) + len(line[4]) - 1 + downstream

        create_hash_key_chain(_hash, None, line[2], matches, 'score')
        _hash[line[2]][matches]['score'] = 0

        create_hash_key_chain(_hash, None, line[2], matches, 'mature')
        _hash[line[2]][matches]['mature'] = line[0]  # assign unique mature seq

        for sample in samples.keys():
            create_hash_key_chain(hash_sample, _hash[line[2]][matches][
                                  'beg'], sample, line[2], matches, 'beg')
            create_hash_key_chain(hash_sample, _hash[line[2]][matches][
                                  'end'], sample, line[2], matches, 'end')
            create_hash_key_chain(hash_sample, _hash[line[2]][matches][
                                  'score'], sample, line[2], matches, 'score')
            create_hash_key_chain(hash_sample, _hash[line[2]][matches][
                                  'mature'], sample, line[2], matches, 'mature')

        OUT.write('{}\t{}\n'.format(line[2], line[0]))

    print("\n{} mature mappings to precursors\n".format(cx))
    IN.close()
    OUT.close()


def ReadinStarMappingFile(options):
    global name3, _hash, hash_sample, hash_star_sample, hash_star, samples
    line = []
    matches = 0

    IN = open_or_die("{}_mapped.bwt".format(
        name3), 'rb', "Mature mapping file {}_mapped.bwt not found \n".format(name3))
    cx = 0
    ltmp = 'qwertyuiop'
    id1 = ''
    id2 = ''

    while True:
        l = IN.readline()
        if not l:
            break

        id1 = ''
        id2 = ''
        line = re.split('\t', l)

        id1 = line[0]  # mature ID
        id2 = line[2]  # precursor ID

        # remove multiple endings
        id1 = re.sub('*', '', id1)
        id1 = re.sub('-5p', '', id1)
        id1 = re.sub('-3p', '', id1)

        m1 = re.match(r'^(\w+\-\w+\-\w+)\-\d+$', id1)
        if m1:
            id1 = m1.groups()[0]

        m2 = re.match(r'^(\w+\-\w+\-\w+)\-\d+$', id2)
        if m2:
            id2 = m.groups()[0]

        if options.get('-k') == '' and not re.search(id2, id1, re.IGNORECASE) and not re.search(id1, id2, re.IGNORECASE):
            continue

        cx += 1
        create_hash_key_chain(hash_star, 0, line[2], 'c')
        hash_star[line[2]]['c'] += 1

        for sample in samples:
            create_hash_key_chain(hash_star_sample, 0, sample, line[2], 'c')
            hash_star_sample[sample][line[2]]['c'] += 1

        matches = hash_star[line[2]]['c']

        hash_star[line[2]][matches]['beg'] = int(line[3]) - upstream
        if hash_star[line[2]][matches]['beg'] < 0:
            hash_star[line[2]][matches]['beg'] = 0

        create_hash_key_chain(hash_star, line[2], matches, 'end')
        hash_star[line[2]][matches]['end'] = int(
            line[3]) + len(line[4]) - 1 + downstream

        create_hash_key_chain(hash_star, line[2], matches, 'score')
        hash_star[line[2]][matches]['score'] = 0

        create_hash_key_chain(hash_star, line[2], matches, 'mature')
        hash_star[line[2]][matches]['mature'] = line[0]

        for sample in samples:
            create_hash_key_chain(hash_star_sample, sample,
                                  line[2], matches, 'beg')
            hash_star_sample[sample][line[2]][matches][
                'beg'] = hash_star[line[2]][matches]['beg']

            create_hash_key_chain(hash_star_sample, sample,
                                  line[2], matches, 'end')
            hash_star_sample[sample][line[2]][matches][
                'end'] = hash_star[line[2]][matches]['end']

            create_hash_key_chain(hash_star_sample, sample,
                                  line[2], matches, 'score')
            hash_star_sample[sample][line[2]][matches][
                'score'] = hash_star[line[2]][matches]['score']

            create_hash_key_chain(hash_star_sample, sample, line[
                                  2], matches, 'mature')
            hash_star_sample[sample][line[2]][matches][
                'mature'] = hash_star[line[2]][matches]['mature']

    print('\n{} star mappings to precursors\n'.format(cx))


def ReadinReadsMappingFile(options):
    global mapcounts, species, _hash, hash_sample, hash_star, hash_star_sample, total, total_t
    line = []
    rb = 0
    _re = 0
    scores = []
    len_sc = 0
    ids = {}

    IN = open_or_die('{}_mapped.bwt'.format(
        name2), 'rb', 'Reads mapping File {}_mapped.bwt not found \n'.format(name2))
    while True:
        l = IN.readline()
        if not l:
            break

        m = re.match(r'^(\S+)', l)
        if m:
            m = m.groups()
            create_hash_key_chain(mapcounts, 0, m[0])
            mapcounts[m[0]] += 1
    IN.close()

    OUT = open_or_die('read_occ', 'w+',
                      'Could not create file with read_occ\n')
    for k in sorted(mapcounts.keys()):
        OUT.write('{}\t{}\n'.format(k, mapcounts[k]))
    OUT.close()

    IN = open_or_die('{}_mapped.bwt'.format(
        name2), 'rb', 'Reads mapping File {}_mapped.bwt not found \n'.format(name2))
    matched = 0
    sample = None

    while True:
        matched = 0
        l = IN.readline()
        if not l:
            break
        line = re.split(r'\t', l)
        if species != 'none':
            if not re.search(species, line[2]):
                continue

        if options.get('-U') == '' and mapcounts[line[0]] > 1:
            continue

        rb = int(line[3])
        _re = int(line[3]) + len(line[4]) - 1

        if not hash_value_true(_hash, line[2], 'c'):
            continue

        for i in range(1, _hash[line[2]]['c'] + 1):
            if options.get('-w') == '':
                scores = re.split(r'x', line[0])

                mm = re.match(r'^(\S\S\S)_', scores[0])
                if mm:
                    sample = mm.groups()[0]

                len_sc = int(scores[-1])

                if options.get('-W') == '':
                    len_sc /= mapcounts[line[0]]

                create_hash_key_chain(_hash, 0, line[2], i, 'score')
                _hash[line[2]][i]['score'] += len_sc

                create_hash_key_chain(
                    hash_sample, 0, sample, line[2], i, 'score')
                hash_sample[sample][line[2]][i]['score'] += len_sc

                create_hash_key_chain(total, 0, sample)
                total[sample] += len_sc

                total_t += len_sc

                matched = 1

                create_hash_key_chain(_hash, 0, line[2], 'r')
                _hash[line[2]]['r'] += len_sc

                create_hash_key_chain(hash_sample, 0, sample, line[2], 'r')
                hash_sample[sample][line[2]]['r'] += len_sc
            else:

                if rb >= _hash[line[2]][i]['beg'] and _re <= _hash[line[2]][i]['end']:

                    scores = re.split(r'x', line[0])

                    mm = re.match(r'^(\S\S\S)_', scores[0])
                    if mm:
                        sample = mm.groups()[0]

                    len_sc = int(scores[-1])

                    if options.get('-W') == '':
                        len_sc /= mapcounts[line[0]]

                    create_hash_key_chain(_hash, 0, line[2], i, 'score')
                    _hash[line[2]][i]['score'] += len_sc

                    create_hash_key_chain(
                        hash_sample, 0, sample, line[2], i, 'score')
                    hash_sample[sample][line[2]][i]['score'] += len_sc

                    create_hash_key_chain(total, 0, sample)
                    total[sample] += len_sc
                    total_t += len_sc

                    matched = 1

                    create_hash_key_chain(_hash, 0, line[2], 'r')
                    _hash[line[2]]['r'] += len_sc

                    create_hash_key_chain(hash_sample, 0, sample, line[2], 'r')
                    hash_sample[sample][line[2]]['r'] += len_sc

        for i in range(1, hash_star[line[2]]['c'] + 1):

            if options.get('-w') == '':
                scores = re.split('x', line[0])

                m = re.match(r'^(\S\S\S)_', scores[0])
                if m:
                    sample = m.groups()[0]

                len_sc = int(scores[-1])

                if options.get('-W') == '':
                    len_sc /= mapcounts[line[0]]

                create_hash_key_chain(hash_star, 0, line[2], i, 'score')
                hash_star[line[2]][i]['score'] += len_sc

                create_hash_key_chain(
                    hash_star_sample, 0, sample, line[2], i, 'score')
                hash_star_sample[sample][line[2]][i]['score'] += len_sc

                create_hash_key_chain(_hash, 0, line[2], 'r')
                _hash[line[2]]['r'] += len_sc

                create_hash_key_chain(hash_sample, 0, sample, line[2], 'r')
                hash_sample[sample][line[2]]['r'] += len_sc

            else:

                if rb >= hash_star[line[2]][i]['beg'] and _re <= hash_star[line[2]][i]['end']:

                    scores = re.split(r'x', line[0])

                    mm = re.match(r'^(\S\S\S)_', scores[0])
                    if mm:
                        sample = mm.groups()[0]

                    len_sc = int(scores[-1])

                    if options.get('-W') == '':
                        len_sc /= mapcounts[line[0]]

                    create_hash_key_chain(hash_star, 0, line[2], i, 'score')
                    hash_star[line[2]][i]['score'] += len_sc

                    create_hash_key_chain(
                        hash_star_sample, 0, sample, line[2], i, 'score')
                    hash_star_sample[sample][line[2]][i]['score'] += len_sc

                    create_hash_key_chain(total, 0, sample)
                    total[sample] += len_sc
                    total_t += len_sc

                    create_hash_key_chain(_hash, 0, line[2], 'r')
                    _hash[line[2]]['r'] += len_sc

                    create_hash_key_chain(hash_sample, 0, sample, line[2], 'r')
                    hash_sample[sample][line[2]]['r'] += len_sc

                    matched = 1

        if not matched:
            scores = re.split('x', line[0])
            m = re.match(r'^(\S\S\S)_', scores[0])
            if m:
                sample = m.groups()[0]

            len_sc = int(scores[-1])

            if options.get('-W') == '':
                len_sc /= mapcounts[line[0]]

            create_hash_key_chain(_hash, 0, line[2], 'r')
            _hash[line[2]]['r'] += len_sc

            create_hash_key_chain(hash_sample, 0, sample, line[2], 'r')
            hash_sample[sample][line[2]]['r'] += len_sc

    IN.close()
    # print_stderr(_hash, "\n\n")


def PrintExpressionValues():
    global outdir, _hash, species, hash_star, hash_sample, hash_star_sample
    OUT1 = open_or_die('{}/miRNA_expressed.csv'.format(outdir),
                       'w+', 'can ont create {}/miRNA_expressed.csv'.format(outdir))
    OUT2 = open_or_die('{}/miRNA_not_expressed.csv'.format(outdir),
                       'w+', 'can ont create {}/miRNA_not_expressed.csv'.format(outdir))

    OUT1.write('#miRNA\tread_count\tprecursor\n')
    OUT2.write('#miRNA\tread_count\n')

    seen = None
    not_seen = None

    for k in sorted(_hash.keys()):
        if species != 'none':
            if not re.search(species, k):
                continue

        for i in range(1, _hash[k]['c'] + 1):
            if _hash[k][i]['score'] > 0:
                OUT1.write('{}\t{}\t{}\n'.format(
                    _hash[k][i]['mature'],
                    _hash[k][i]['score'],
                    k
                ))
            else:
                OUT2.write('{}\t0\n'.format(_hash[k][i]['mature']))

            if _hash[k][i]['score'] == 0:
                OUT1.write('{}\t{}\t{}\n'.format(
                    _hash[k][i]['mature'],
                    _hash[k][i]['score'],
                    k
                ))

        for k in sorted(hash_star.keys()):
            if species != 'none':
                if not re.search(species, k):
                    continue

            for i in range(1, hash_star[k]['c'] + 1):
                if hash_star[k][i]['score'] > 0:
                    OUT1.write('{}\t{}\t{}\n'.format(
                        hash_star[k][i]['mature'],
                        hash_star[k][i]['score'],
                        k
                    ))
                else:
                    OUT2.write('{}\t0\n'.format(hash_star[k][i]['mature']))

        print_stderr('Expressed miRNAs are written to {}/miRNA_expressed.csv not expressed miRNAs are written to {}/miRNA_not_expressed.csv\n'.format(outdir, outdir))
    OUT1.close()
    OUT2.close()


def PrintExpressionValuesSamples(options):
    global total, total_t, _hash, hash_star, hash_sample, hash_star_sample, ctime

    total_t = 1000000
    OUTG = open_or_die2(
        'miRNAs_expressed_all_samples_{}.csv'.format(ctime), 'w+')
    OUTG.write('#miRNA\tread_count\tprecursor\ttotal')

    for sample in sorted(hash_sample.keys()):
        if re.search('config', sample):
            continue
        OUTG.write('\t{}'.format(sample))

    for sample in sorted(hash_sample.keys()):
        if re.search('config', sample):
            continue
        OUTG.write('\t{}(norm)'.format(sample))

    OUTG.write('\n')

    for k in _hash.keys():
        if species != 'none':
            if not re.search(species, k):
                continue

        if options.get('-w') == '':
            i = 1
            OUTG.write('{}\t{}\t{}\t{}'.format(
                _hash[k][i]['mature'],
                _hash[k][i]['score'],
                k,
                _hash[k][i]['score']
            ))

            for sample in sorted(hash_sample.keys()):
                if re.search(r'config', sample):
                    continue
                if hash_sample[sample][k][i]['score'] > 0:
                    OUTG.write('\t{}'.format(
                        hash_sample[sample][k][i]['score']))
                else:
                    OUTG.write('\t0')

            for sample in sorted(hash_sample.keys()):
                if re.search(r'config', sample):
                    continue
                if hash_sample[sample][k][i]['score'] > 0:
                    OUTG.write('\t{0:.2f}'.format(
                        total_t *
                        hash_sample[sample][k][i][
                            'score'] / float(total[sample])
                    ))
                else:
                    OUTG.write('\t0')

            OUTG.write('\n')

        else:
            for i in range(1, _hash[k]['c'] + 1):

                OUTG.write('{}\t{:.2f}\t{}\t{:.2f}'.format(
                    _hash[k][i]['mature'],
                    _hash[k][i]['score'],
                    k,
                    _hash[k][i]['score']
                ))

                for sample in sorted(hash_sample.keys()):
                    # ATTR: Perl will create a dict automatically when
                    # accessing via key chain
                    create_hash_key_chain(
                        hash_sample, 0, sample, k, i, 'score')
                    create_hash_key_chain(hash_sample, 0, sample, k, i, 'r')

                    if re.search('config', sample):
                        continue

                    if hash_sample[sample][k][i]['score'] > 0:
                        OUTG.write('\t{:.2f}'.format(
                            hash_sample[sample][k][i]['score']
                        ))
                    else:
                        OUTG.write('\t0')

                for sample in sorted(hash_sample.keys()):
                    if re.search('config', sample):
                        continue

                    if hash_sample[sample][k][i]['score'] > 0:
                        OUTG.write('\t{:.2f}'.format(
                            total_t *
                            hash_sample[sample][k][i][
                                'score'] / float(total[sample])
                        ))
                    else:
                        OUTG.write('\t0')

                OUTG.write('\n')

    if options.get('-s'):
        for k in sorted(hash_star.keys()):
            if species != 'none':
                if re.search(species, k):
                    continue

            star_keys = len(hash_star_sample.keys())
            if star_keys > 0:
                if options.get('-w') == '':
                    i = 1
                    OUTG.write('{}\t{}\t{}\t{}'.format(
                        hash_star[k][i]['mature'],
                        hash_star[k][i]['score'],
                        k,
                        hash_star[k][i]['score']
                    ))

                    for sample in sorted(hash_star_sample.keys()):
                        if re.search('config', sample):
                            continue

                        if hash_star_sample[sample][k][i]['score'] > 0:
                            OUTG.write('\t{}'.format(
                                hash_star_sample[sample][k][i]['score']
                            ))
                        else:
                            OUTG.write('\t0')

                    for sample in sorted(hash_star_sample.keys()):
                        if re.search('config', sample):
                            continue

                        if hash_star_sample[sample][k][i]['score'] > 0:
                            OUTG.write('\t{0:.2f}'.format(
                                total_t *
                                hash_star_sample[sample][k][i][
                                    'score'] / float(total[sample])
                            ))
                        else:
                            OUTG.write('\t0')

                    OUTG.write('\n')

                else:
                    for i in range(1, hash_star[k]['c'] + 1):
                        OUTG.write('{}\t{}\t{}\t{}'.format(
                            hash_star[k][i]['mature'],
                            hash_star[k][i]['score'],
                            k,
                            hash_star[k][i]['score']
                        ))

                        for sample in sorted(hash_star_sample.keys()):
                            if re.search('config', sample):
                                continue

                            if hash_star_sample[sample][k][i]['score'] > 0:
                                OUTG.write('\t{}'.format(
                                    hash_star_sample[sample][k][i]['score']
                                ))
                            else:
                                OUTG.write('\t0')

                        for sample in sorted(hash_star_sample.keys()):
                            if re.search('config', sample):
                                continue

                            if hash_star_sample[sample][k][i]['score'] > 0:
                                OUTG.write('\t{0:.2f}'.format(
                                    total_t *
                                    hash_star_sample[sample][k][i][
                                        'score'] / float(total[sample])
                                ))
                            else:
                                OUTG.write('\t0')

                        OUTG.write('\n')

    OUTG.close()


def CreateOutputMRD():
    global mapcounts, outdir, name1, name2, name3, name0, name4, _hash, hash_star_sample, hash_star, hash_sample, species
    exprs = {}
    ex = []

    IN = open_or_die('{}/miRNA_expressed.csv'.format(outdir), 'rb',
                     'file {}/miRNA_expressed.csv not found\n'.format(outdir))
    while True:
        l = IN.readline().strip()
        if not l:
            break

        if re.search('precursorID', l):
            continue

        ex = re.split('\t', l)
        # precursor-entitiy = count
        create_hash_key_chain(exprs, ex[1], ex[2], ex[0])

    IN.close()

    os.chdir(outdir)
    OUT = open_or_die('miRBase.mrd', 'w+',
                      'could not create file {}/miRBase.mrd\n'.format(outdir))
    mature, star, reads = (None, None, None)

    mature = os.system(
        'convert_bowtie_output.py {}_mapped.bwt > {}_mapped.arf'.format(name1, name1))

    tmp = []
    line = []
    ltmp = None
    id1 = ''
    id2 = ''

    IN = open_or_die('{}_mapped.arf'.format(name1), 'rb',
                     'can not open {}_mapped.arf'.format(name1))
    while True:
        l = IN.readline()
        if not l:
            break

        line = re.split(r'\t', l)
        if species != 'none':
            if not re.search(species, line[5]):
                continue

        id1 = line[0]  # mature id
        id2 = line[5]  # precursor id

        id1 = re.sub('-5p', '', id1)
        id2 = re.sub('-3p', '', id2)

        m = re.match(r'^(\w+\-\w+\-\w+)\-\d+$', id2)
        if m:
            id2 = m.groups()[0]

        if options.get('-k') != '' and not re.search(id2, id1, re.IGNORECASE) and not re.search(id1, id2, re.IGNORECASE):
            continue

        if line[5] not in _hash.keys() or 'struct' not in _hash[line[5]].keys() or not _hash[line[5]]['struct']:
            tmp = ssplit('f' * len(_hash[line[5]]['seq']))
        else:
            tmp = ssplit(_hash[line[5]]['struct'])

        if re.search(r'5p', line[0]):
            for i in range(int(line[7]) - 1, int(line[8])):
                tmp[i] = '5'
        elif re.search(r'3p', line[0]):
            for i in range(int(line[7]) - 1, int(line[8])):
                tmp[i] = '3'
        else:
            for i in range(int(line[7]) - 1, int(line[8])):
                tmp[i] = 'M'

        _hash[line[5]]['struct'] = ''.join(tmp)

    IN.close()

    if options.get('-s'):
        star = os.popen('convert_bowtie_output.py {}_mapped.bwt > {}_mapped.arf'.format(
            name3, name3)).read()
        IN = open_or_die('{}_mapped.arf'.format(name3), 'rb',
                         'can not open {}_mapped.arf'.format(name3))
        while True:
            l = IN.readline()
            if not l:
                break

            line = re.split('\t', l)
            if species != 'none':
                if not re.search(species, line[5]):
                    continue

            id1 = ''
            id2 = ''

            id1 = line[0]
            id2 = line[5]

            id1 = re.sub(r'\*', '', id1)
            id1 = re.sub(r'-5p', '', id1)
            id1 = re.sub(r'-3p', '', id1)

            m = re.match(r'^(\w+\-\w+\-\w+)\-\d+$', id1)
            if m:
                id1 = m.groups()[0]

            m = re.match(r'^(\w+\-\w+\-\w+)\-\d+$', id2)
            if m:
                id2 = m.groups()[0]

            if options.get('-k') != '' and not re.search(id2, id1, re.IGNORECASE) and not re.search(id1, id2, re.IGNORECASE):
                continue

            tmp = ssplit(_hash[line[5]]['struct'])

            for i in range(int(line[7]) - 1, int(line[8])):
                tmp[i] = 'S'

            _hash[line[5]]['struct'] = ''.join(tmp)

        IN.close()

    found = 0
    loop = 0

    if options.get('-s'):
        for k in _hash.keys():
            loop = ssplit(_hash[k]['struct'])
            for i in range(0, len(loop)):
                if loop[i] != 'f' and found == 0:
                    found = 1
                elif loop[i] == 'f' and found == 1:
                    loop = 1
                elif loop[i] != 'f' and found == 1 and loop == 1:
                    found = 2

            if found == 2:
                found = 0
                loop = 0
                for i in range(0, len(loop)):
                    if loop[i] != 'f' and found == 0:
                        found = 1
                    elif loop[i] == 'f' and found == 1 and loop < 2:
                        loop[i] = 'l'
                        loop = 1
                    elif loop[i] != 'f' and found == 1 and loop > 0:
                        found = 2

            _hash[k]['struct'] = ''.join(loop)

    reads = os.popen('convert_bowtie_output.py {}_mapped.bwt > {}_mapped.arf'.format(
        name2, name2)).read()
    IN = open_or_die('{}_mapped.arf'.format(name2), 'rb',
                     'can not open {}_mapped.arf'.format(name2))

    counter = 0
    col1_width = 40
    spacer = None
    rseq = []
    qseq = []
    rc = None

    while True:
        l = IN.readline().strip()
        if not l:
            break

        if re.match(r'^\s*$', l):
            continue

        line = re.split('\t', l)
        if options.get('-U') == '' and mapcounts[line[0]] > 1:
            continue

        if species != 'none':
            if not re.search(species, line[5]):
                continue

        rc = re.split('_x', line[0])

        create_hash_key_chain(_hash, 0, line[5], 'reads', 'c')
        # sum up total read count for this precursor
        _hash[line[5]]['reads']['c'] += float(rc[1])
        counter = _hash[line[5]]['reads']['c']

        create_hash_key_chain(_hash, None, line[5], 'reads', counter)
        _hash[line[5]]['reads'][counter] = l

    print('after READS READ IN thing\n\n')

    for k1 in sorted(_hash.keys()):
        if species != 'none':
            if not re.search(species, k1):
                continue

        _hash[k1]['seq'] = tr(_hash[k1]['seq'], 'ACGTU', 'acguu')
        _str = os.popen('echo {} | RNAfold'.format(
            _hash[k1]['seq'])).read().split('\n')
        _str1 = re.split('\s+', _str[1])

        OUT.write('>{}\n'.format(k1))
        spacer = ' ' * (col1_width - len('total read count'))

        # ATTR: Perl will create a dict automatically when accessing via key
        # chain
        create_hash_key_chain(_hash, 0, k1, 'r')

        OUT.write('total read count{}{}\n'.format(spacer, _hash[k1]['r']))

        for k2 in _hash[k1].keys():
            if not re.match(r'^\d+$', str(k2)):
                continue

            mat = _hash[k1][k2]['mature']
            spacer = ' ' * (col1_width - len('{} read count'.format(mat)))
            OUT.write('{} read count{}{}\n'.format(
                mat, spacer, _hash[k1][k2]['score']))

        for k2 in hash_star[k1].keys():
            if not re.match(r'^\d+$', k2):
                continue

            mat = hash_star[k1][k2]['mature']
            spacer = ' ' * (col1_width - len('{} read count'.format(mat)))
            OUT.write('{} read count{}{}\n'.format(
                mat, spacer, hash_star[k1][k2]['score']))

        spacer = ' ' * (col1_width - len('remaining read count'))

        rrc = _hash[k1]['r']
        for i in range(1, _hash[k1]['c'] + 1):
            rrc -= _hash[k1][i]['score']

        for i in range(1, hash_star[k1]['c'] + 1):
            rrc -= hash_star[k1][i]['score']

        OUT.write('remaining read count{}{}\n'.format(spacer, rrc))

        spacer = ' ' * (col1_width - len('exp'))
        OUT.write('exp{}{}\n'.format(spacer, _hash[k1]['struct']))

        spacer = ' ' * (col1_width - len('pri_seq'))
        OUT.write('pri_seq{}{}\n'.format(spacer, _hash[k1]['seq']))

        spacer = ' ' * (col1_width - len('pri_struct'))
        OUT.write('pri_struct{}{}\t#MM\n'.format(spacer, _str1[0]))

        reads_arr = []
        reads_hash = {}

        pseq = ssplit(_hash[k1]['seq'].lower())

        # ATTR: Perl will create a dict automatically when accessing via key
        # chain
        create_hash_key_chain(_hash, {}, k1, 'reads')

        for k2 in _hash[k1]['reads'].keys():
            if k2 == 'c':
                continue
            reads_arr.append(_hash[k1]['reads'][k2])

        lr = len(reads_arr)
        rrc = 0
        for item in reads_arr:
            m = re.match(
                r'^(\S+)\s+\d+\s+\d+\s+\d+\s+(\S+)\s+\S+\s+\d+\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\d+).+$', item)
            if m:
                rrc += 1
                reads_hash[rrc] = {}
                reads_hash[rrc]['id'] = m.groups()[0]
                reads_hash[rrc]['seq'] = m.groups()[1]
                reads_hash[rrc]['beg'] = m.groups()[2]
                reads_hash[rrc]['end'] = m.groups()[3]
                reads_hash[rrc]['mm'] = m.groups()[4]

        sitems = sorted(reads_hash.items(), key=lambda x: x[1]['beg'])
        skeys = [x[0] for x in sitems]
        elist = []

        # ATTR : perl did not check the bound of array, if no elem for the
        # specified index, then undef is return
        if len(skeys) == 0:
            first = None
        else:
            # all keys that have same begin position should match this value
            first = reads_hash[skeys[0]]['beg']
        rorder = {}  # temporary hash to store all keys with same begin position

        for j in range(0, len(skeys)):
            if reads_hash[skeys[j]]['beg'] == first:
                # insert key and end position to hash
                rorder[j] = reads_hash[skeys[j]]['end']
            else:
                first = reads_hash[skeys[j]]['beg']
                # sort hash keys by end position
                for rk in hash_sort_key(rorder, lambda x: x[1]):
                    elist.append(skeys[rk])  # attend keys to elist

                for rk in rorder.keys():  # delete hash
                    del rorder[rk]

                rorder[j] = reads_hash[skeys[j]]["end"]

        ritems = sorted(rorder.items(), key=lambda x: x[1])
        rkeys = [x[0] for x in ritems]
        for rk in rkeys:
            elist.append(skeys[rk])

        for ei in elist:
            rseq = reads_hash[ei]['seq'].lower()
            rseq = tr(rseq, 't', 'u')
            bef = '.' * (int(reads_hash[ei]['beg']) - 1)
            after = '.' * (int(_hash[k1]['end']) - int(reads_hash[ei]['end']))
            spacer = ' ' * (col1_width - len(reads_hash[ei]['id']))
            sread = ssplit(rseq)

            bshift = 0
            rseq = ''

            for i in range(0, len(sread)):
                if not pseq[i + int(reads_hash[ei]['beg']) - 1]:
                    pass
                else:
                    if pseq[i + int(reads_hash[ei]['beg']) - 1] == sread[i]:
                        rseq += sread[i]
                    else:
                        sread[i] = sread[i].upper()
                        rseq += sread[i].upper()

            OUT.write('{}{}{}{}{}\t{}\n'.format(
                reads_hash[ei]['id'],
                spacer,
                bef,
                rseq,
                after,
                reads_hash[ei]['mm']
            ))

        OUT.write('\n\n\n')

    OUT.close()
    os.chdir('../../')


if __name__ == '__main__':

    opts, argss = getopt.getopt(
        sys.argv[1:], 'p:m:r:s:t:y:dokuc:nxg:e:f:vjwT:PWU')
    options = dict(opts)

    for line in __DATA.split('\n'):
        line = line.strip()
        m = re.match('^(\S+)\s+(\S+)$', line)
        if m:
            u = m.groups()[0].lower()
            v = m.groups()[1].lower()
            u = re.sub(' ', '', u)
            v = re.sub(' ', '', v)
            organisms[u] = v
            rorganisms[v] = u

    mismatches = 1
    if options.get('-g'):
        mismatches = options.get('-g')

    threads = 1
    if options.get('-T'):
        threads = options.get('-T')

    upstream = 2
    downstream = 5

    if options.get('-e'):
        upstream = options.get('-e')

    if options.get('-f'):
        downstream = options.get('-f')

    if options.get('-u') == '':
        print_stderr(
            '\n\nAllowed species arguments that have an entry at UCSC\n\n')
        for key in organisms:
            print_stderr('{}\t{}\n'.format(key, organisms[key]))

        die('\n')

    if not options.get('-p') or not options.get('-r'):
        die(usage)

    if options.get('-w') != '' and not options.get('-m'):
        die(usage)

    if options.get('-w') == '':
        options['-m'] = '{}.dummy'.format(options.get('-p'))
        IN = open_or_die(options.get('-p'), 'rb',
                         'cannot open {}'.format(options.get('-p')))
        OUT = open_or_die(options.get('-m'), 'w+',
                          'cannot open {}'.format(options.get('-m')))
        while True:
            line = IN.readline()
            if not line:
                break

            if re.search(r'>', line):
                OUT.write('')
            else:
                OUT.write(substr(line, 0, 18))
                OUT.write('\n')

        IN.close()
        OUT.close()

        options['-e'] = 0
        options['-f'] = 0

    opt_m = ''
    if options.get('-t'):
        species = options.get('-t').lower()
        species = re.sub(' ', '', species)

        if species in rorganisms.keys() and rorganisms[species]:
            species = rorganisms[species]
        elif species in organisms.keys() and organisms[species]:
            pass
        else:
            print_stderr('\n\nThe species {} you specified is not available\nallowed species are\n'.format(
                options.get('-t')))
            os.system('quantifier.py -u')
            sys.exit(1)

        opt_m = '-m {}'.format(species)

    if options.get('-y'):
        ctime = options.get('-y')

    opt_d = ''
    if options.get('-d') == '':
        opt_d = '-d'

    # sort pdf reads sample
    opt_o = ''
    if not options.get('-o') == '':
        opt_o = '-o'

    name0, path0, extension0 = fileparse(options.get('-p'), '\..*')
    name1, path1, extension1 = (None, None, None)
    name1, path1, extension1 = fileparse(options.get('-m'), '\..*')
    name2, path2, extension2 = fileparse(options.get('-r'), '\..*')
    name3, path3, extension3 = (None, None, None)

    name0 += extension0
    name1 += extension1
    name2 += extension2

    if options.get('-s'):
        if file_s(options.get('-s')):
            name3, path3, extension3 = fileparse(options.get('-s'), '\..*')
            name3 += extension3
        else:
            print_stderr('The file {} is empty or not found. It will be ignored for this analysis'.format(
                options.get('-s')))
            options['-s'] = 0

    _dir = 'expression_analyses'
    if not os.path.isdir(_dir):
        os.mkdir(_dir)

    outdir = '{}/{}_{}'.format(_dir, _dir, ctime)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    IN = open_or_die(options.get('-r'), 'rb',
                     'File {} not found\n'.format(options.get('-r')))
    line = IN.readline()
    if not re.match(r'^>\S\S\S_\d+_x\d+', line):
        die('''\n{} ids do not have the correct format

it must have the id line >SSS_INT_xINT\n
SSS is a three letter code indicating the sample origin
INT is just a running number
xINT is the number of read occurrences\n\n

You can use the mapper.py module to create such a file from a fasta file with
mapper.py {} -c -m -s {}.collapsed

See also the mapper.py help for more information on preprocessing input files.

'''.format(options.get('-r'), options.get('-r'), options.get('-r')))

        IN.close()

        samples = {}
        print_stderr('getting samples and corresponding read numbers\n\n')

    # if 0:
    #     if os.path.isfile(options.get('-c')):
    #         IN = open_or_die(options.get('-c'), 'rb', 'can not open {}'.format(options.get('-c')))
    #         while True:
    #             line = IN.readline().strip()
    #             if not line:
    #                 break

    #             m = re.findall(r'^(\S+)\s+(\S+)$', line)
    #             if m:
    #                 FIN = open_or_die(m[0], 'rb', 'file {} not found\n'.format(m[0]))
    #                 while True:
    #                     line2 = FIN.readline()
    #                     if not line2:
    #                         break

    #                     if re.search(r'^>', line2):
    #                         continue

    #                     mm = re.findall(r'^>(\S\S\S)_\S+_x(\d+)$', line2)
    #                     if mm:
    #                         samples[m[0]] = m[1]

    #                 FIN.close()
    #         IN.close()
    #     else:
    #         if options.get('-c') not in samples.keys():
    #             samples[options.get('-c')] = 0

    #         samples[options.get('-c')] += 1

    # for key in samples.keys():
    #     print_stderr('{}\t{} reads\n'.format(key, samples[key]))

    # print_stderr('\n\n')

    # convert input files to bowtie accepting format
    if options.get('-n') != '':
        print_stderr('Converting input files\n')
        ConvertFastaFile(options.get('-p'), name0, 'precursor', species)
        ConvertFastaFile(options.get('-m'), name1, 'mature', species)
        ConvertFastaFile(options.get('-r'), name2, "", "")

        if options.get('-s'):
            ConvertFastaFile(options.get('s'), name3, 'star', species)

    os.chdir(outdir)
    if options.get('-x') != '':
        Mapping(options)

    # now analyze expression file
    print_stderr('analyzing data\n')
    ReadinPrecursorFile()

    ReadinMatureMappingFile(options)

    if options.get('-s'):
        ReadinStarMappingFile(options)

    ReadinReadsMappingFile(options)
    os.chdir('../../')

    PrintExpressionValues()
    PrintExpressionValuesSamples(options)

    print_stderr('\nCreating miRBase.mrd file\n\n')
    if options.get('-j') == '':
        die('exiting here and not creating mor.mrd file\nif you want this created do not specify option -j\n')

    CreateOutputMRD()

    opt_l = '-l'
    if options.get('-k') == '':
        opt_l = ''

    command = ''
    opt_P = ''
    if options.get('-P') == '':
        opt_P = '-P'

    opt_W = ''
    if options.get('-W') == '':
        opt_W = '-W {}/read_occ'.format(outdir)

    startf = ''
    if options.get('-s'):
        startf = '-j {}/{}_mapped.arf'.format(outdir, name3)

    if species in organisms.keys() and organisms[species]:
        command = 'make_html2.py -q {}/miRBase.mrd -k {} -t {} -y {} {} {} -i {}/{}_mapped.arf {} {} {} -M miRNAs_expressed_all_samples_{}.csv {} {}'.format(
            outdir,
            name1,
            organisms[species],
            ctime,
            opt_d,
            opt_o,
            outdir,
            name1,
            startf,
            opt_l,
            opt_m,
            ctime,
            opt_P,
            opt_W
        )

        print_stderr('{}\n'.format(command))

        os.system(command)

    elif species in rorganisms.key() and rorganisms[species]:

        command = "make_html2.py -q {}/miRBase.mrd -k {} -t {} -y {} {} {} -i {}/{}_mapped.arf {} {} {} -M miRNAs_expressed_all_samples_{}.csv {} {}".format(
            outdir,
            name1,
            species,
            ctime,
            opt_d,
            opt_o,
            outdir,
            name1,
            startf,
            opt_l,
            opt_m,
            ctime,
            opt_P,
            opt_W
        )

        print_stderr("{}\n".format(command))

        os.system(command)

    else:

        command = "make_html2.py -q {}/miRBase.mrd -k {} -y {} {} {} -i {}/{}_mapped.arf {} {} {} -M miRNAs_expressed_all_samples_{}.csv {} {}".format(
            outdir,
            name1,
            ctime,
            opt_d,
            opt_o,
            outdir,
            name1,
            startf,
            opt_l,
            opt_m,
            ctime,
            opt_P,
            opt_W
        )

        print_stderr("{}\n".format(command))

        os.system(command)

    sys.exit(0)
