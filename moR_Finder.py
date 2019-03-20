#!/usr/bin/env python
from __future__ import print_function

import os
import sys

from port import die, print_stderr


'''
mRNA analysis with moRNA Finder

reads_mr.fa                     mrna reads added to read file end
mrna_places.bed                 coordinates of mrna containing hairpin in the chromosome
cel_cluster.fa                  the genome while which is already in tutorial_dir folder
mature_ref_this_species.fa      fasta file with the reference miRBase mature miRNAs for the species, is already in tutorial_dir folder

1. bowtie-build cel_cluster.fa cel_cluster

2. bedtools getfasta -fi cel_cluster.fa -bed mrna_places.bed -name -fo mrna_hp.fa  * 2

3. mapper.py reads_mr.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT -l 18 -m -p cel_cluster -s reads_collapsed_mr.fa -t reads_collapsed_mr_vs_genome.arf -v

4. quantifier.py -p mrna_hp.fa -m mature_ref_this_species.fa -r reads_collapsed_mr.fa -t cel -y 27_10

'''

usage = '''
Usage:

    {} reads_mr.fa mrna_places.bed reads_collapsed_mr.fa reads_collapsed_mr_vs_genome.arf clip_seq timestamp

    reads_mr.fa                         mrna reads added to read file end
    mrna_places.bed                     coordinates of mrna containing hairpin in the chromosome
    reads_collapsed_mr.fa               print processed reads to this file
    reads_collapsed_mr_vs_genome.arf    print read mappings to this file
    clip_seq                            clip 3' adapter sequence
    timestamp                           The timestamp for this run

    eg.

    {} reads_mr.fa mrna_places.bed reads_collapsed_mr.fa reads_collapsed_mr_vs_genome.arf TCGTATGCCGTCTTCTGCTTGT timestamp

'''.format(sys.argv[0], sys.argv[0])


def test_bedtools():
    ret = os.system('which bedtools > /dev/null 2>&1')
    if ret:
        die('''bedtools is not installed in your environment.
For ubuntu/debian, run `apt-get install bedtools`.
For Redhat/CentOS, run `yum install BEDTools`.
For MacOS, run `brew install bedtools`.

If you need more information about how to install bedtools, please visit http://bedtools.readthedoc.org.
''')


def run_bowtie_cmd(fafile, fa_prefix):
    cmd = 'bowtie-build {} {}'.format(fafile, fa_prefix)
    print_stderr(cmd, '\n')
    ret = os.system(cmd)
    if ret:
        die('Run bowtie failed.\n')


def run_bedtools_cmd(run_times, fafile, mrna_places, mrna_hp):
    for i in range(run_times):
        cmd = 'bedtools getfasta -fi {} -bed {} -name -fo {}'.format(
            fafile, mrna_places, mrna_hp)
        print_stderr(cmd, '\n')
        ret = os.system(cmd)
        if ret:
            die('Run bedtools failed.\n')


def run_mapper(reads_file, clip_seq, fa_prefix, collapsed_file, mapping_file, ):
    cmd = 'mapper.py {} -c -j -k {} -l 18 -m -p {} -s {} -t {} -v'.format(
        reads_file, clip_seq, fa_prefix, collapsed_file, mapping_file)
    print_stderr(cmd, '\n')
    ret = os.system(cmd)
    if ret:
        die('Run mapper.py failed.\n')


def run_quantifier(mrna_hp, mature_this_file, collapsed_file, timestamp):
    cmd = 'quantifier.py -p {} -m {} -r {} -t cel -y {}'.format(
        mrna_hp, mature_this_file, collapsed_file, timestamp)
    print_stderr(cmd, '\n')
    ret = os.system(cmd)
    if ret:
        die('Run quantifier.py failed.\n')


if __name__ == '__main__':
    if len(sys.argv) < 7:
        die(usage)

    reads_file = sys.argv[1]
    mrna_places = sys.argv[2]
    collapsed_file = sys.argv[3]
    mapping_file = sys.argv[4]
    clip_seq = sys.argv[5]
    timestamp = sys.argv[6]

    fa_prefix = 'cel_cluster'
    fafile = '{}.fa'.format(fa_prefix)
    mrna_hp = 'mrna_hp.fa'
    mature_this_file = 'mature_ref_this_species.fa'

    # Test bedtools installed or not
    test_bedtools()

    # run bowtie
    run_bowtie_cmd(fafile, fa_prefix)

    # Run bedtools
    run_bedtools_cmd(2, fafile, mrna_places, mrna_hp)

    # run mapper
    run_mapper(reads_file, clip_seq, fa_prefix, collapsed_file, mapping_file)

    # run quantifier
    run_quantifier(mrna_hp, mature_this_file, collapsed_file, timestamp)
