#!/usr/bin/env python
from __future__ import print_function

import argparse
import datetime
import getopt
import os
import random
import re
import sys
from shutil import rmtree

from port import (die, esplit, open_or_die, print_stderr, print_stderr_red,
                  substr)

usage = '''{} input_file_reads

This script takes as input a file with deep sequencing reads (these can be in
different formats, see the options below). The script then processes the reads
and/or maps them to the reference genome, as designated by the options given.
Options:

Read input file:
-a              input file is seq.txt format
-b              input file is qseq.txt format
-c              input file is fasta format
-e              input file is fastq format
-d              input file is a config file (see moRNA Finder documentation).
                options -a, -b, -c or -e must be given with option -d.

Preprocessing/mapping:
-g              three-letter prefix for reads (by default 'seq')
-h              parse to fasta format
-i              convert rna to dna alphabet (to map against genome)
-j              remove all entries that have a sequence that contains letters
                other than a,c,g,t,u,n,A,C,G,T,U,N
-k seq          clip 3' adapter sequence
-l int          discard reads shorter than int nts, default = 18
-m              collapse reads

-p genome       map to genome (must be indexed by bowtie-build). The 'genome'
                string must be the prefix of the bowtie index. For instance, if
                the first indexed file is called 'h_sapiens_37_asm.1.ebwt' then
                the prefix is 'h_sapiens_37_asm'.
-q              map with one mismatch in the seed (mapping takes longer)

-r int          a read is allowed to map up to this number of positions in the genome
                default is 5

Output files:
-s file         print processed reads to this file
-t file         print read mappings to this file

Other:
-u              do not remove directory with temporary files
-v              outputs progress report

-n              overwrite existing files

-o              number of threads to use for bowtie

Example of use:

{} reads_seq.txt -a -h -i -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p
    h_sapiens_37_asm -s reads.fa -t reads_vs_genome.arf -v
'''.format(sys.argv[0], sys.argv[0])

# GLOBAL VARIABLES
threads = 1
orig_file_reads = 0
mismatches_seed = 0
_dir = ''


def check_line(line):
    if re.search(r'-h\s+\d/', line) or re.search(r'-h\s+\w/', line):
        die("option -h should not be given with an integer or string\n")

    if re.search(r'-i\s+\d/', line) or re.search(r'-i\s+\w/', line):
        die("option -i should not be given with an integer or string\n")

    if re.search(r'-j\s+\d/', line) or re.search(r'-j\s+\w/', line):
        die("option -j should not be given with an integer or string\n")

    if re.search(r'-m\s+\d/', line) or re.search(r'-m\s+\w/', line):
        die("option -m should not be given with an integer or string\n")

    if re.search(r'-q\s+\d/', line) or re.search(r'-q\s+\w/', line):
        die("option -q should not be given with an integer or string\n")


def check_file_format_and_option(file_reads, aFormat):
    print_stderr('\n')
    warning = '''\n\n***** Please check if the option you used (options $format) designates the correct format of the supplied reads file $file *****\n\n
[options]
-a              input file is seq.txt format
-b              input file is qseq.txt format
-c              input file is fasta format
-e              input file is fastq format
-d              input file is a config file (see moRNA Finder documentation).
                options -a, -b, -c or -e must be given with option -d.
'''
    line = None
    if aFormat == 'a':
        i = 0
        IN = open_or_die(
            file_reads, 'rb', 'Cannot open file {} supplied by option -a\n'.format(file_reads))
        while True:
            l = IN.readline().strip()
            if not l:
                break
            i += 1
            line = esplit(l)
            # $#line != 4
            if len(line) != 5:
                die('The seq.txt file does not contain 5 columns. Please make sure to follow the _seq.txt file format conventions\n{}'.format(warning))

            if i == 4:
                break
        IN.close()
    elif aFormat == 'b':
        IN = open_or_die(
            file_reads, 'rb', 'Cannot open qseq.txt file {} supplied by option -b\n'.format(file_reads))
        i = 0
        mes = 'Please make sure your file is in accordance with the qses.txt format specifications\n'
        while True:
            l = IN.readline().strip()
            if not l:
                break
            i += 1
            line = esplit(l)

            if len(line) != 11:
                die('The qseq.txt file does not contain 11 columns but {}. Please make sure to follow the qseq.txt file format conventions\n{}'.format(
                    len(line), warning))

            if not re.search(r'^\S+', line[9]):
                die('The sequence field in the qseq.txt file is invalid. Please make sure to follow the qseq.txt file format conventions\n{}'.format(warning))

            if i == 4:
                break
        IN.close()
    elif aFormat == '-c':
        IN = open_or_die(file_reads, 'rb',
                         'Cannot open FASTA file supplied by option -c\n')
        i = 0
        mes = 'Please make sure your file is in accordance with the fasta format specifications and does not contain whitespace in IDs or sequences'
        while True:
            l = IN.readline().strip()
            if not l:
                break
            i += 1
            if i == 1:
                if not re.search(r'^>\S+$', l):
                    die("First line of FASTA reads file is not in accordance with the fasta format specifications\n{}\n{}".format(
                        mes, warning))
            if i == 2:
                if not re.search(r'^\S+$', l):
                    die("Second line of FASTA reads file contains whitespace in sequence\n{}\n".format(
                        mes))
            if i == 3:
                if not re.search(r'^>\S+$', l):
                    die("Second ID line of FASTA reads file is not in accordance with the fasta format specifications\n{}\n{}".format(
                        mes, warning))
            if i == 4:
                if not re.search(r'^\S+$', l):
                    die("Secdond sequence line of FASTA reads file contains whitespace in sequence\n{}\n{}".format(
                        mes, warning))

            if i == 4:
                break
        IN.close()
    elif aFormat == '-e':
        IN = open_or_die(file_reads, 'rb',
                         'Cannot open FASTQ file supplied by option -e\n')
        i = 0
        mes = 'Please make sure your file is in accordance with the FASTQ format specifications'
        while True:
            l = IN.readline().strip()
            if not l:
                break
            i += 1
            if i == 1:
                if not re.search(r'^@\S+', l):
                    die("First line of FASTQ reads file is not in accordance with the fastq format specifications\n{}\n{}".format(
                        mes, warning))
            if i == 2:
                if re.search(r'^\S+$', l):
                    die("Second line of FASTQ reads file contains whitespace in sequence\n{}\n{}".format(
                        mes, warning))
            if i == 3:
                if re.search(r'^\+', l):
                    die("Third line of FASTQ reads file does not start with a '+' character.\n{}\n{}".format(mes, warning))
            if i == 4:
                if re.search(r'^\S+$', l):
                    die("Fourth line of FASTQ reads file contains whitespace\n{}\n{}".format(
                        mes, warning))

            if i == 4:
                break


def test_prefix(prefix):
    if not (re.search(r'^\w\w\w$', prefix) and not re.search(r'_', prefix)):
        die('prefix $prefix does not contain exactly three alphabet letters\n')


def printErr():
    print_stderr_red('\nError: ')


def check_options(options, file_reads):
    formats = 0
    if options.get('-a') == '':
        formats += 1
        if options.get('-d') != '':
            check_file_format_and_option(file_reads, 'a')

    if options.get('-b') == '':
        formats += 1
        if options.get('-d') != '':
            check_file_format_and_option(file_reads, 'b')

    if options.get('-c') == '':
        formats += 1
        if options.get('-d') != '':
            check_file_format_and_option(file_reads, 'c')

    if options.get('-e') == '':
        formats += 1
        if options.get('-d') != '':
            check_file_format_and_option(file_reads, 'e')

    if formats != 1:
        die('exactly one input format (-a, -b , -e or -c) must be designated\n')

    # check if file supplied matches option, otherwise quit
    processing_steps = 0

    if options.get('-h') == '':
        processing_steps += 1

    if options.get('-i') == '':
        processing_steps += 1

    if options.get('-j') == '':
        processing_steps += 1

    if options.get('-k'):
        processing_steps += 1

    if options.get('-l'):
        processing_steps += 1

    if options.get('-m') == '':
        processing_steps += 1

    if options.get('-p'):
        processing_steps += 1

    if processing_steps <= 0:
        die('at least one processing/mapping step (-h, -i, -j, -k, -l, -m or -p) must be designated\n')

    file_output = 0
    if '-o' in options.keys():
        if not(re.search(r'\d+', options.get('-o')) and int(options.get('-o')) > 0):
            die('options -o must be positive integer\n')

    if options.get('-s'):
        file_output += 1

    if options.get('-t'):
        file_output += 1

    if file_output <= 0:
        die('at least one output file (-s or -t) must be designated\n')

    if options.get('-s') and os.path.exists(options.get('-s')) and not options.get('-n') == '':
        die("file {} already exists\n".format(options.get('-s')))

    if options.get('-t') and os.path.exists(options.get('-t')) and not options.get('-n') == '':
        die("file {} already exists\n".format(options.get('-t')))

    if options.get('-a') == '' or options.get('-b') == '' or options.get('-e') == '' and options.get('-h') != '':
        die("raw illumina output must be parsed to fasta format with options -h\n")

    if options.get('-c') == '' and options.get('-h') == '':
        die("input file is already designated as a fasta file, so option -h should not be used\n")

    if options.get('-c') == '' and not(options.get('-i') == '' or options.get('-j') == '' or options.get('-k') or options.get('-l') or options.get('-m') == '' or options.get('-p')):
        die("at least one processing/mapping step (-i, -j, -k, -l, -m or -p) must be designated\n")

    if options.get('-d') == '' and not(options.get('-a') != '' or options.get('-b') == '' or options.get('-c') == '' or options.get('-e')) == '':
        die("option -d must be given with option -a, -b, -c or -e \n")

    if options.get('-d') == '' and options.get('-g'):
        die("option -d and -g are mutually exclusive. If -d is given, the prefixes must be contained in the config file\n")

    if options.get('-g'):
        test_prefix(options.get('-g'))

    if options.get('-i') == '' and not(options.get('-c') == '' or options.get('-h') == ''):
        die("option -i must be used on reads in fasta format or with option -h \n")

    if options.get('-j') == '' and not(options.get('-c') == '' or options.get('-h') == ''):
        die("option -j must be used on reads in fasta format or with option -h \n")

    if options.get('-k') and not(options.get('-c') == '' or options.get('-h') == ''):
        die("option -k must be used on reads in fasta format or with option - h\n")

    if options.get('-l') and not(options.get('-c') == '' or options.get('-h') == ''):
        die("option -l must be used on reads in fasta format or with option -h \n")

    if options.get('-m') == '' and not(options.get('-c') == '' or options.get('-h') == ''):
        die("option -m must be used on reads in fasta format or with option -h \n")

    if options.get('-p') and not(options.get('-c') == '' or options.get('-h') == ''):
        die("option -p must be used on reads in fasta format or with option -h\n")

    if options.get('-q') == '' and not(options.get('p')):
        die("option -q must be given with option -p\n")

    if options.get('-s') and not(options.get('-h') == '' or options.get('-i') == '' or options.get('-j') == '' or options.get('-k') or options.get('-l') or options.get('-m') == '' or options.get('-p')):
        die("at least one processing step (-h, -i, -j, -k, -l, -m or -p) must be designated if processed file should be output (-s)\n")

    if options.get('-t') and not(options.get('-p')):
        die("reads must be mapped (-p) if mappings are to be output (-t)\n")

    if options.get('-k') and re.search(r'^-', options.get('-k')):
        die("please make sure that the adapter sequence designated with the -k option is correct\n")

    if options.get('-l') and re.search(r'^-', options.get('-l')):
        die("please make sure that the int given with the -l option is correct\n")

    if options.get('-p') and re.search(r'ebwt$', options.get('-p')):
        die("please make sure that you are using the -p option correctly.\nThe argument given after -p must be the _prefix_ of the bowtie\nindexed files and should not contain 'ebwt'. For instance,\nif the first indexed file is called 'h_sapiens_37_asm.1.ebwt'\nthen the prefix is 'h_sapiens_37_asm'.\n")

    if options.get('-p') and re.search(r'^-', options.get('-p')):
        die("please make sure that the genome index designated with the -p option is correct\n")

    # added by SM to check if bowtie is installed when reads should be mapped
    # to genome
    if options.get('-p'):
        # TODO: my $binst=`bowtie --version 2>&1`;
        binst = os.system('bowtie --version 2>&1')
        if binst:
            printErr()
            die("Bowtie mapping tool not installed.\nPlease download from http://downloads.sourceforge.net/project/bowtie-bio/bowtie/ the latest version and install it.\n\n")

    if options.get('-s') and re.search(r'^-', options.get('-s')):
        die("please make sure that the output file designated with the -s option is correct\n")

    if options.get('-t') and re.search(r'^-', options.get('-t')):
        die("please make sure that the output file designated with the -t option is correct\n")


def make_dir_tmp(pref, MAP):
    _today = datetime.datetime.now()
    _time = _today.strftime('%d_%m_%y_%H_%M_%S')
    MAP.write('\ntimestamp:\t{}\n\n'.format(_time))
    num = random.random()
    chance = substr(str(num), 2, 10)

    _dir = 'dir_mapper{}_{}_{}'.format(pref, chance, _time)
    os.mkdir(_dir)
    return _dir


def cat_to(file_1, file_2):
    OUT = open_or_die(file_2, 'a', 'cannot print to {}'.format(file_2))
    IN = open_or_die(file_1, 'rb', 'cannot read from {}'.format(file_1))

    while True:
        line = IN.readline()
        if not line:
            break
        OUT.write(line)

    IN.close()
    OUT.close()


def process_reads(file_reads_latest, prefix, MAP):
    global _dir, orig_file_reads
    orig_file_reads = file_reads_latest
    m = re.search(r'([_\-.a-zA-Z0-9]+)$', file_reads_latest)
    if m:
        orig_file_reads = m.groups()[0]

    _dir = make_dir_tmp("_{}_{}".format(prefix, orig_file_reads), MAP)

    # parse solexa to fasta
    if options.get('-h') == '':
        if options.get('-e') == '':
            MAP.write('parsing fastq to fasta format\n')

            if options.get('-v') == '':
                print_stderr('parsing fastq to fasta format\n')

            cmd = 'fastq2fasta.py {} > {}/reads.fa\n'.format(
                file_reads_latest, _dir)
            MAP.write(cmd)
            ret_format = os.system(cmd)
            file_reads_latest = '{}/reads.fa'.format(_dir)
        else:
            MAP.write('parsing Solexa / Illumina output to fasta format\n')

            if options.get('-v') == '':
                print_stderr(
                    'parsing Solexa / Illumina output to fasta format\n')

            line = 'illumina_to_fasta.py {}'.format(file_reads_latest)

            if options.get('-b') == '':
                line += ' -a'

            cmd = '{} > {}/reads.fa\n'.format(line, _dir)
            MAP.write(cmd)

            ret_format = os.system(cmd)
            file_reads_latest = '{}/reads.fa'.format(_dir)

    # RNA to DNA
    if options.get('-i') == '':
        MAP.write('converting rna to dna alphabet\n')

        if options.get('-v') == '':
            print_stderr('converting rna to dna alphabet\n')

        ret_rna2dna = os.system(
            'rna2dna.py {} > {}/reads_dna.fa'.format(file_reads_latest, _dir))
        file_reads_latest = '{}/reads_dna.fa'.format(_dir)

    # discard entries that contain non-canonical letters
    if options.get('-j') == '':
        MAP.write('discarding sequences with non-canonical letters\n')
        if options.get('-v') == '':
            print_stderr('discarding sequences with non-canonical letters\n')

        cmd = 'fastaparse.py {} -b > {}/reads_letters.fa 2>{}/reads_discarded.fa\n'.format(
            file_reads_latest, _dir, _dir)
        MAP.write(cmd)
        ret_clip = os.system(cmd.strip())
        file_reads_latest = '{}/reads_letters.fa'.format(_dir)

    # clip 3' adapters
    if options.get('-k'):
        MAP.write("clipping 3' adapters\n")
        if options.get('-v') == '':
            print_stderr("clipping 3' adapters\n")

        cmd = 'clip_adapters.py {} {} > {}/reads_clip.fa\n'.format(
            file_reads_latest, options.get('-k'), _dir)
        MAP.write(cmd)
        ret_clip = os.system(cmd.strip())
        file_reads_latest = '{}/reads_clip.fa'.format(_dir)

    if options.get('-l'):
        MAP.write('discarding short reads\n')

        if options.get('-v') == '':
            print_stderr('discarding short reads\n')

        cmd = 'fastaparse.py {} -a {} > {}/reads_no_short.fa 2>{}/reads_too_short.fa\n'.format(
            file_reads_latest, options.get('-l'), _dir, _dir)
        MAP.write(cmd)
        ret_rem_short = os.system(cmd.strip())
        file_reads_latest = '{}/reads_no_short.fa'.format(_dir)

    # collapse reads
    if options.get('-m') == '':
        MAP.write('collapsing reads\n')

        if options.get('-v') == '':
            print_stderr('collapsing reads\n')

        cmd = 'collapse_reads_md.py {} {} > {}/reads_nr.fa\n'.format(
            file_reads_latest, prefix, _dir)
        MAP.write(cmd)
        ret_collapse = os.system(cmd)
        file_reads_latest = '{}/reads_nr.fa'.format(_dir)

    # printing reads
    if options.get('-s'):
        cat_to(file_reads_latest, options.get('-s'))

    return file_reads_latest


def map_reads(file_reads_latest, MAP, options):
    global mismatches_seed, threads, orig_file_reads
    # map reads to genome
    MAP.write('mapping reads to genome index\n')
    if options.get('-v') == '':
        print_stderr('mapping reads to genome index\n')

    file_genome_latest = options.get('-p')
    mapping_loc = 5
    if '-r' in options.keys():
        mapping_loc = options.get('-r')

    cmd = 'bowtie -p {} -f -n {} -e 80 -l 18 -a -m {} --best --strata {}  --al {}/{}_mapped --un {}/{}_not_mapped  {} {}/mappings.bwt 2>bowtie.log\n\n'.format(
        threads,
        mismatches_seed,
        mapping_loc,
        file_genome_latest,
        _dir,
        orig_file_reads,
        _dir,
        orig_file_reads,
        file_reads_latest,
        _dir
    )
    MAP.write(cmd)
    ret_mapping = os.system(cmd.strip())
    file_mapping_latest = '{}/mappings.bwt'.format(_dir)

    cmd = 'convert_bowtie_output.py {} > {}/mappings.arf\n'.format(
        file_mapping_latest, _dir)
    MAP.write(cmd)
    ret_parse_to_arf = os.system(cmd.strip())
    file_mapping_latest = '{}/mappings.arf'.format(_dir)

    # trim unmapped nts in the 3' end
    MAP.write("trimming unmapped nts in the 3' ends\n")
    if options.get('-v') == '':
        print_stderr("trimming unmapped nts in the 3' ends\n")

    cmd = 'parse_mappings.py {} -j > {}/mappings_trim.arf\n\n'.format(
        file_mapping_latest, _dir)
    MAP.write(cmd)
    ret_trim = os.system(cmd.strip())
    file_mapping_latest = '{}/mappings_trim.arf'.format(_dir)

    if options.get('-v') == '':
        cat_to(file_mapping_latest, options.get('-t'))

    return file_mapping_latest


def remove_dir_tmp(options, MAP):
    global _dir
    if options.get('-u') != '':
        MAP.write('remove tmp dir\nrmtree({})\n\n'.format(_dir))
        rmtree(_dir)


def handle_one_file(file_reads, prefix, MAP, options):
    file_reads_latest = process_reads(file_reads, prefix, MAP)

    if options.get('-p'):
        file_mapping_latest = map_reads(file_reads_latest, MAP, options)

    remove_dir_tmp(options, MAP)


def handle_config_file(file_reads, MAP, options):
    FILE = open_or_die(
        file_reads, 'rb', 'can not open {}\n'.format(file_reads))
    while True:
        l = FILE.readline()
        if not l:
            break

        m = re.match(r'(^\S+)\s+(\S+)\s*.*$', l)
        if m:
            m = m.groups()
            file_reads = m[0]
            prefix = m[1]

            if (len(file_reads) < len(prefix)):
                file_reads = m[1]
                prefix = m[0]

            test_prefix(prefix)

            MAP.write("\nhandling file '{}' with prefix '{}'\n".format(
                file_reads, prefix))

            # check if files in config file are in accordance with option
            # specified
            if options.get('-a') == '':
                check_file_format_and_option(file_reads, 'a')
            if options.get('-b') == '':
                check_file_format_and_option(file_reads, 'b')
            if options.get('-c') == '':
                check_file_format_and_option(file_reads, 'c')
            if options.get('-e') == '':
                check_file_format_and_option(file_reads, 'e')

            if options.get('-v') == '':
                print_stderr("\nhandling file '{}' with prefix '{}'\n".format(
                    file_reads, prefix))

            handle_one_file(file_reads, prefix, MAP, options)

    FILE.close()


def read_stats(options):
    _hash = {}
    count = 0
    k2 = {}
    IN = open_or_die(options.get('-s'), 'rb',
                     'No reads file in fasta format given\n')
    while True:
        line = IN.readline()
        if not line:
            break

        m = re.match(r'^>*((\S\S\S)\S+_x(\d+))', line)
        if m:
            m = m.groups()

            try:
                if _hash[m[0]]:
                    continue
            except KeyError:
                pass

            # ATTR: Performance issue below, use logic above
            # if m[0] in _hash.keys() and _hash[m[0]]:
            #     continue

            _hash[m[0]] = 1
            count += int(m[2])

            if m[1] not in k2.keys():
                k2[m[1]] = 0
            k2[m[1]] += int(m[2])
    IN.close()

    _hash2 = {}
    count2 = 0
    k22 = {}

    print_stderr('Mapping statistics\n')
    IN = open_or_die(options.get('-t'), 'rb', 'No mapping file given\n')
    while True:
        line = IN.readline()
        if not line:
            break
        m = re.match(r'^>*((\S\S\S)\S+_x(\d+))', line)
        if m:
            m = m.groups()
            if m[0] in _hash2.keys() and _hash2[m[0]]:
                continue
            _hash2[m[0]] = 1
            count2 += int(m[2])

            if m[1] not in k22.keys():
                k22[m[1]] = 0

            k22[m[1]] += int(m[2])
    IN.close()

    print_stderr('\n#desc\ttotal\tmapped\tunmapped\t%mapped\t%unmapped\n')
    print_stderr("total: {}\t{}\t{}\t".format(count, count2, count - count2))
    print_stderr("{0:.3f}\t{1:.3f}\n".format(
        count2 / float(count), 1 - (count2 / float(count))))

    for k in k2.keys():
        print_stderr('{}: {}\t{}\t{}\t'.format(
            k, k2[k], k22[k], k2[k] - k22[k]))
        print_stderr('{0:.3f}\t{1:.3f}\n'.format(
            float(k22[k]) / k2[k], 1 - (float(k22[k]) / k2[k])))

if __name__ == '__main__':
    # create a log file for the mapper.py
    # the latest run of mapper will be on top of the log file
    #
    #
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(-1)

    if os.path.exists('mapper.log'):
        os.system('mv mapper.log mapper.log_bak')
    else:
        os.system('touch mapper.log_bak')

    MAP = open_or_die('mapper.log_tmp', 'w+',
                      'could not create mapper.log_tmp\n')
    # cdir = os.path.dirname(os.path.realpath(__file__))
    cdir = os.path.abspath('.')

    MAP.write('current dir:\t{}\n'.format(cdir))
    MAP.write('mapper command:\t{} {}\n'.format(
        sys.argv[0], ' '.join(sys.argv[1:])))

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('input_file_reads', help='input file')
    args = parser.parse_args(sys.argv[1:2])

    file_reads = args.input_file_reads

    if not os.path.exists(file_reads):
        print('No config or reads file could be found')
        print(usage)
        sys.exit(-1)

    line = '\t'.join(sys.argv[1:])
    check_line(line)

    opt, argss = getopt.getopt(sys.argv[2:], 'abcdeg:hijk:l:mp:qs:t:uvnr:o:')
    options = dict(opt)

    if options.get('-s') and os.path.exists(options.get('-s')) and options.get('-n') == '':
        os.system('rm {}'.format(options.get('-s')))
    if options.get('-t') and os.path.exists(options.get('-t')) and options.get('-n') == '':
        os.system('rm {}'.format(options.get('-t')))

    if not options.get('-l'):
        options['-l'] = 18

    check_options(options, file_reads)

    if '-o' in options.keys():
        threads = options.get('-o')

    cores = os.popen('grep -ic ^processor /proc/cpuinfo').read()
    if not re.search(r'^\d+$', cores):
        cores = os.popen('sysctl -n hw.physicalcpu').read()
        if not re.search(r'^\d+$', cores):
            cores = os.popen('sysctl -n hw.logicalcpu').read()

    if not re.search(r'^\d+$', cores):
        cores = 1

    if threads > cores:
        print_stderr(
            'More threads specified than cores on the system. Reducing the number of threads to {}\n'.format(cores))
        threads = cores

    if options.get('-q') == '':
        mismatches_seed = 1

    prefix_global = 'seq'

    if options.get('-g'):
        prefix_global = options.get('-g')

    if options.get('-d'):
        handle_config_file(file_reads)
    else:
        handle_one_file(file_reads, prefix_global, MAP, options)

    MAP.write('#' * 60)
    MAP.write('\n\n')

    MAP.close()

    os.system('cat mapper.log_tmp mapper.log_bak > mapper.log')
    os.system('rm mapper.log_tmp mapper.log_bak')

    # print_stderr('read statistics\n')
    if options.get('-s') and options.get('-t'):
        read_stats(options)
