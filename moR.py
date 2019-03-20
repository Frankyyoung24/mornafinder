#!/usr/bin/env python
from __future__ import print_function

import argparse
import getopt
import math
import os
import re
import shutil
import sys
import time

from port import (create_hash_key_chain, die, esplit, file_s, fileparse,
                  hash_defined, hash_sort_key, hash_value_true, localtime,
                  max2, min2, open_or_die, open_or_die2, pprint, print_stderr,
                  print_stderr_red, pround, revcom, ssplit, substr, tr)

version = '1.3.1'

print('''



#####################################
#                                   #
# moRNA Finder {}                   #
#                                   #
# last change: 10/09/2016           #
#                                   #
#####################################

'''.format(version))

usage = '''

{} reads genome mappings miRNAs_ref/none miRNAs_other/none precursors/none 2>report.log

This script enacts the moRNA Finder pipeline. The input files are:

reads         deep sequences in fasta format. The identifier should contain a prefix, a running
              number and a '_x' to indicate the number of reads that have this sequence.
              There should be no redundancy in the sequences.
genome        genome contigs in fasta format. The identifiers should be unique.
mappings      file_reads mapped against file_genome. The mappings should be in arf format.
              For details on the format see the documentation.
miRNAs_ref    miRBase miRNA sequences in fasta format. These should be the known mature
              sequences for the species being analyzed.
miRNAs_other  miRBase miRNA sequences in fasta format. These should be the pooled known
              mature sequences for 1-5 species closely related to the species being
              analyzed.

precursors    miRBase miRNA precursor sequences in fasta format. These should be the known precursor
              sequences for the species being analyzed.

The output files produced are:

result.html   a html table giving an overview of novel and known miRNAs detected in the
              data. The table is hyperlinked to pdfs showing the signature and structure of
              each hairpin.
result.csv    spread-sheet format of results.html
survey.csv    spread-sheet of prediction accuracy for all score-cutoffs between -10 and 10.
output.mrd    text output of the reported hairpins.

Options:

-a int        minimum read stack height that triggers analysis. Using this option disables
              automatic estimation of the optimal value and all detected precursors are analyzed

-g int        maximum number of precursors to analyze when automatic excision gearing is used.
              default=50.000, if set to -1 all precursors will be analyzed

-b int        minimum score cut-off for predicted novel miRNAs to be displayed in the overview
              table. This score cut-off is by default 0.

-c            disable randfold analysis

-d            disable pdf generation

-t species    species being analyzed - this is used to link to the appropriate UCSC browser entry

-u            output list of UCSC browser species that are supported and exit

-v            remove directory with temporary files

-o            do not sort aligned reads in pdf files by sample, only used if multiple samples were used as input (see Readme for mor information)

-s file       File with known miRBase star sequences

-r string     Prefix for output file names

-z tag        Additional tag appended to current time stamp

-P            use this switch if mature_ref_miRNAs contain miRBase v18 identifiers (5p and 3p) instead of previous ids from v17

Example of use:

{} reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs -t Mouse 2>report.log

'''.format(sys.argv[0], sys.argv[0])

warning = '''\n**********\nThe first three arguments to {} must be files while arguments 4-6 can be files or must be designated as 'none'. Examples:\n
{} reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs

or

{} reads.fa genome.fa reads_vs_genome.arf none none none


Please check if the supplied files exist and the command call to {} is correct
*******************\n'''.format(sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

# -q file  miRBase.mrd file from quantifier module to show miRBase miRNAs in data that were not scored by moRNA Finder

# hardcoded variables
read_align_mismatches = 1

_dir = None
dir_tmp = None
max_pres = 50000
# minimal precursor length, used for precheck of precursor file
minpreslen = 40
stack_height_min = None
ctime = None
ltime = None
scripts = None
sTimeG = None
stimeG = None
command_line = None
file_reads = None
file_genome = None
file_reads_vs_genome = None
file_mature_ref_this_species = None
file_mature_ref_other_species = None
file_precursors = None
parsed_arf = None
options = {}


def myTime(ctime):
    (sec, minute, hour, day, month, year, temp1, temp2, temp3) = localtime(ctime)
    year += 1900
    month += 1
    ret = "{:02}_{:02}_{:02}_t_{:02}_{:02}_{:02}".format(
        day, month, year, hour, minute, sec)
    return ret


def test_first_argument():
    if len(sys.argv) < 2:
        die(usage)

    if sys.argv[1] == '-u':
        os.system('make_html.py -u -y 1')
        sys.exit(0)

    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        die(usage)


def make_dir_tmp():
    global _dir, ltime, dir_tmp
    # make temporary directory
    if not os.path.isdir('moR_runs'):
        os.mkdir('moR_runs')

    _dir = "moR_runs/run_{}".format(ltime)

    print_stderr("mkdir {}\n\n".format(_dir))
    os.mkdir(_dir)

    dir_tmp = "{}/tmp".format(_dir)

    os.mkdir(dir_tmp)


def checkBIN(a, b):
    global _dir
    e = os.system(
        "{} 1>{}/tmp/binaries 2>{}/tmp/binaries2".format(a, _dir, _dir))

    IN = open_or_die("{}/tmp/binaries".format(_dir), 'rb',
                     'can not open {}/tmp/binaries'.format(_dir))
    found = 1

    while True:
        line = IN.readline()
        if not line:
            break

        if re.search(b, line):
            found = 0

    IN.close()

    if found:
        IN = open_or_die('{}/tmp/binaries2'.format(_dir), 'rb',
                         'can not open {}/tmp/binaries2'.format(_dir))
        while True:
            line = IN.readline()
            if not line:
                break

            if re.search(b, line):
                found = 0

    IN.close()
    return found


def test_installed_binaries():
    global scripts
    stdm = 'If you used the install.py script make sure that you started a complete new shell window after installation.\nIf this did not help please restart youer workstation.\n\n'

    ret = None

    ret = checkBIN("bowtie --version", "version")
    if ret:
        die("Error: \tbowtie not found\nCheck if bowtie is correctly installed and all Pathes were set correctly.\n", stdm)

    ret = checkBIN("RNAfold -h", "gamma")
    if ret:
        die("Error: \tRNAfold not found\nCheck if RNAfold is correctly installed and all Pathes were set correctly.\n", stdm)

    ret = checkBIN("randfold", "let7")
    if ret:
        die("Error: \trandfold not found\nCheck if randfold is correctly installed and all Pathes were set correctly.\n", stdm)

    # TODO, perl PDF lib requirement
    # ret = checkBIN("perl -e \'use PDF::API2; print \"installed\";\'","installed")
    # if ret:
    # die "Error: \tPerl PDF::API2 package not found\nCheck if the perl
    # PDF::API2 package is correctly installed and all Pathes were set
    # correctly.\n$stdm" if($ret);

    if not os.path.isfile('{}/Rfam_for_moR.fa'.format(scripts)):
        die("Error:\t Rfam_for_moR.fa not found in your moRNA Finder scripts directory\nPlease copy this file from the moRNA Finder archive to your moRNA Finder scripts directory\n\n")

    return 0


def cwd():
    '''
    Get current dir path
    '''
    return os.path.realpath('.')


def printUsedParameters(options):
    global _dir, ltime, command_line, file_reads, file_genome, file_reads_vs_genome
    fname = '{}/run_{}_parameters'.format(_dir, ltime)
    OUT = open_or_die(fname, 'w+', 'can not open {}\n'.format(fname))
    OUT.write("Start: {}\n".format(ltime))
    OUT.write("Script\t{}\n".format(sys.argv[0]))
    OUT.write("args {}\n".format(command_line))
    OUT.write("dir_with_tmp_files\tdir_moR_{}\n".format(ltime))

    d = cwd()

    OUT.write("dir\t{}\n".format(d))
    OUT.write("file_reads\t{}\n".format(file_reads))
    OUT.write("file_genome\t{}\n".format(file_genome))
    OUT.write("file_reads_vs_genome\t{}\n".format(file_reads_vs_genome))
    OUT.write("file_mature_ref_this_species\t{}\n".format(
        file_mature_ref_this_species))
    OUT.write("file_mature_ref_other_species\t{}\n".format(
        file_mature_ref_other_species))

    if options.get('-a'):
        OUT.write("option -a =\t{}\n".format(options.get('-a')))
    if options.get('-b'):
        OUT.write("option -b =\t{}\n".format(options.get('-b')))
    if options.get('-c') == '':
        OUT.write("option -c =\t{}\n".format(options.get('-c')))
    if options.get('-t'):
        OUT.write("option -t =\t{}\n".format(options.get('-t')))
    if options.get('-v') == '':
        OUT.write("option -v =\tused\n")
#    if($options{'q'}){print OUT "option{q} =\t$options{'q'}\n";}

    OUT.close()


def printErr():
    print_stderr_red("Error: ")


def test_input_files():
    global file_reads, file_reads_vs_genome, file_genome, file_precursors, minpreslen, file_mature_ref_other_species, file_mature_ref_this_species
    IN = open_or_die2(file_reads, 'rb')
    line = IN.readline().strip()
    if not re.search(r'^>\S+', line):
        printErr()
        die("The first line of file $file_reads does not start with '>identifier'\nReads file {} is not a valid fasta file\n\n".format(file_reads))

    if re.search(r'\s', line):
        printErr()
        die('File {} contains not allowed characters in sequences\nAllowed characters are ACGTUN\nReads file {} is not a fasta file\n\n'.format(
            file_reads, file_reads))

    line = IN.readline()
    if not re.search(r'^[ACGTUNacgtun]*$', line):
        printErr()
        die('File {} contains not allowed characters in sequences\nAllowed characters are ACGTUN\nReads file {} is not a fasta file\n\n'.format(
            file_reads, file_reads))

    IN.close()

    IN = open_or_die2(file_genome, 'rb')
    line = IN.readline().strip()
    if not re.search(r'>\S+', line):
        printErr()
        die("The first line of file {} does not start with '>identifier'\nGenome file {} is not a fasta file\n\n".format(
            file_genome, file_genome))

    if re.search(r'\s', line):
        printErr()
        die('Genome file {} has not allowed whitespaces in its first identifier\n\n'.format(
            file_genome))

    # get genome ids
    tmps = os.popen('grep ">" {}'.format(file_genome)).read().strip()
    genomeids = dict(map(lambda x: (x, 1), re.split("\n", tmps)))

    line = IN.readline()
    if not re.search(r'^[ACGTUNacgtun]*$', line):
        printErr()
        die('File {} contains not allowed characters in sequences\nAllowed characters are ACGTUN\nGenome file {} is not a fasta file\n\n'.format(
            file_genome, file_genome))

    IN.close()

    IN = open_or_die2(file_reads_vs_genome, 'rb')
    line = IN.readline()
    if not re.search(r'^(\S+_x\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+(\d+)\s*([mDIM]*)$', line):
        printErr()
        die('Mapping file {} is not in arf format\n\nEach line of the mapping file must consist of the following fields\nreadID_wo_whitespaces  length  start  end read_sequence genomicID_wo_whitspaces length  start   end     genomic_sequence  strand  #mismatches editstring\nThe editstring is optional and must not be contained\nThe readID must end with _xNumber and is not allowed to contain whitespaces.\nThe genomeID is not allowed to contain whitespaces.'.format(file_reads_vs_genome))

    IN.close()

    # get ids from arf file and compare them with ids from the genome file
    tmps = os.popen(
        'cut -f6 {}|sort -u'.format(file_reads_vs_genome)).read().strip()
    for s in re.split("\n", tmps):
        if not genomeids.get(">{}".format(s)):
            die("The mapped reference id {} from file {} is not an id of the genome file {}\n\n".format(
                s, file_reads_vs_genome, file_genome))

    if not re.search('none', file_mature_ref_this_species):
        IN = open_or_die2(file_mature_ref_this_species, 'rb')
        line = IN.readline().strip()
        if not re.search(r'>\S+', line):
            printErr()
            die("The first line of file {} does not start with '>identifier'\nmiRNA reference this species file {} is not a fasta file\n\n".format(
                file_mature_ref_this_species, file_mature_ref_this_species))

        if re.search(r'\s', line):
            printErr()
            die("miRNA reference this species file {} has not allowed whitespaces in its first identifier\n\n".format(
                file_mature_ref_this_species))

        line = IN.readline()
        if not re.search(r'^[ACGTUNacgtun]*$', line):
            printErr()
            die("File {} contains not allowed characters in sequences\nAllowed characters are ACGTUN\nmiRNA reference this species file {} is not a fasta file\n\n".format(
                file_mature_ref_this_species, file_mature_ref_this_species))

        IN.close()

    if not re.search('none', file_mature_ref_other_species):
        IN = open_or_die2(file_mature_ref_other_species, 'rb')
        line = IN.readline().strip()
        if not re.search(r'>\S+', line):
            printErr()
            die("The first line of file {} does not start with '>identifier'\nmiRNA reference this species file {} is not a fasta file\n\n".format(
                file_mature_ref_other_species, file_mature_ref_other_species))

        if re.search(r'\s', line):
            printErr()
            die("miRNA reference this species file {} has not allowed whitespaces in its first identifier\n\n".format(
                file_mature_ref_other_species))

        line = IN.readline()
        if not re.search(r'^[ACGTUNacgtun]*$', line):
            printErr()
            die("File {} contains not allowed characters in sequences\nAllowed characters are ACGTUN\nmiRNA reference this species file {} is not a fasta file\n\n".format(
                file_mature_ref_other_species, file_mature_ref_other_species))

        IN.close()

    if not re.search('none', file_precursors):
        IN = open_or_die2(file_precursors, 'rb')
        line = IN.readline().strip()
        if not re.search(r'>\S+', line):
            printErr()
            die("The first line of file {} does not start with '>identifier'\nmiRNA reference this species file {} is not a fasta file\n\n".format(
                file_precursors, file_precursors))

        if re.search(r'\s', line):
            printErr()
            die("precursor file {} has not allowed whitespaces in its first identifier\n\n".format(
                file_precursors))

        line = IN.readline()
        if not re.search(r'^[ACGTUNacgtun]*$', line):
            printErr()
            die("File {} contains not allowed characters in sequences\nAllowed characters are ACGTUN\nmiRNA reference this species file {} is not a fasta file\n\n".format(
                file_precursors, file_precursors))

        if len(line) < minpreslen:
            printErr()
            die("The precursor file {} does not contain sequences of at least {} nt\nPlease make sure that you provided the correct file and the correct parameter ordering when calling {}\nIf you have precursors with less than {} please use option -p <int> to specify this length\n".format(
                file_precursors,
                minpreslen,
                sys.argv[0],
                minpreslen
            ))

        IN.close()

    # #################################################
    # precheck finished
    # #################################################

    # do stringent testing of all input files
    pprint("#testing input files\n")
    print_stderr("#testing input files\n")

    if not re.search('none', file_mature_ref_this_species):
        start()
        cmd = "sanity_check_mature_ref.py {} 2>&1\n\n".format(
            file_mature_ref_this_species)
        print_stderr(cmd)
        ret_file_mature_ref_this_species = os.popen(cmd).read().strip()

        if ret_file_mature_ref_this_species:
            printErr()
            die("problem with {} {}\n".format(
                file_mature_ref_this_species, ret_file_mature_ref_this_species))
        end()

    if not re.search(r'none', file_mature_ref_other_species):
        start()

        cmd = "sanity_check_mature_ref.py {} 2>&1\n\n".format(
            file_mature_ref_other_species)
        print_stderr(cmd)
        ret_file_mature_ref_other_species = os.popen(cmd).read().strip()

        if ret_file_mature_ref_other_species:
            printErr()
            die("problem with {} {}\n".format(
                file_mature_ref_other_species, ret_file_mature_ref_other_species))
        end()

    cmd = "sanity_check_reads_ready_file.py {} 2>&1\n\n".format(file_reads)
    print_stderr(cmd)
    start()
    ret_test_file_reads = os.popen(cmd).read().strip()

    if ret_test_file_reads:
        printErr()
        die("problem with {} {}\n".format(file_reads, ret_test_file_reads))

    end()

    start()
    cmd = "sanity_check_genome.py {} 2>&1;\n\n".format(file_genome)
    print_stderr(cmd)
    ret_test_file_genome = os.popen(cmd).read().strip()

    if ret_test_file_genome:
        printErr()
        die("problem with {} {}\n".format(file_genome, ret_test_file_genome))

    end()
    start()

    cmd = "sanity_check_mapping_file.py {} 2>&1".format(file_reads_vs_genome)
    print_stderr(cmd)
    ret_test_file_reads_genome = os.popen(cmd).read().strip()

    if ret_test_file_reads_genome:
        printErr()
        die("problem with {} {}\n".format(
            file_reads_vs_genome, ret_test_file_reads_genome))

    end()

    if not re.search('none', file_precursors):
        start()

        cmd = "sanity_check_mature_ref.py {} 2>&1".format(file_precursors)
        print_stderr(cmd)
        ret_file_precursors = os.popen(cmd).read().strip()

        if ret_file_precursors:
            printErr()
            die("problem with {} {}\n".format(
                file_precursors, ret_file_precursors))

        end()

        start()
        if not re.search('none', file_mature_ref_this_species, re.IGNORECASE):
            print_stderr("Quantitation of expressed miRNAs in data\n\n\n")

            species = ''
            if options.get('-t'):
                species = "-t {}".format(options.get('-t'))

            file_star = ''
            if options.get('-s'):
                if file_s(options.get('-s')):
                    file_star = "-s {}".format(options.get('-s'))
                else:
                    print_stderr(
                        "File {} specified by option -s is empty or not found\n".format(options.get('-s')))
                    options['-s'] = 0

            print("#Quantitation of known miRNAs in data\n")
            dopt = ""
            Popt = ""
            if options.get('-d') == '':
                dopt = "-d"
            if options.get('-P') == '':
                Popt = "-P"

            quant = "quantifier.py -p {} -m {} -r {} {} {} -y {} -k {} {}".format(
                file_precursors, file_mature_ref_this_species, file_reads, file_star, species, ltime, dopt, Popt)
            print_stderr(quant, "\n")
            os.system(quant)
            options[
                '-q'] = "expression_analyses/expression_analyses_{}/miRBase.mrd".format(ltime)

            end()
        else:
            print_stderr(
                "Pre-quantitation is skipped caused by missing file with known miRNAs\n\n\n")

    else:
        print_stderr(
            "Pre-quantitation is skipped caused by missing file with known precursor miRNAs\n\n\n")


def start():
    global sTime, stime
    # measuring times
    (second, minute, hour, dayOfMonth, month, yearOffset,
     dayOfWeek, dayOfYear, daylightSavings) = localtime()
    if re.search(r'^\d$', str(second)):
        second = "0{}".format(second)
    sTime = "{}:{}:{}".format(hour, minute, second)
    stime = int(time.time())
    print_stderr("started: {}\n".format(sTime))


def end():
    global sTime, stime
    etime = int(time.time()) - stime
    (second, minute, hour, dayOfMonth, month, yearOffset,
     dayOfWeek, dayOfYear, daylightSavings) = localtime()
    if re.search(r'^\d$', str(second)):
        second = "0{}".format(second)
    eTime = "{}:{}:{}".format(hour, minute, second)

    print_stderr("\nended: {}\ntotal:".format(eTime, int(etime / 3600)),
                 " h:", int((etime % 3600) / 60), "m:", int(etime % 60), "s\n\n")


def rna2dna():
    global dir_tmp, file_mature_ref_other_species, file_mature_ref_this_species, file_precursors
    # process_input mirna files
    if not re.search('none', file_mature_ref_this_species, re.IGNORECASE):
        start()
        # copy file
        (file_mature_ref_this_species_tmp, path0, extension0) = fileparse(
            file_mature_ref_this_species, '\..*')
        cmd = "rna2dna.py {} > {}/{}{}\n\n".format(
            file_mature_ref_this_species,
            dir_tmp,
            file_mature_ref_this_species_tmp,
            extension0
        )
        print_stderr(cmd)
        ret_parse_mature_ref_this_species = os.popen(cmd).read()
        # rename orig file
        file_mature_ref_this_species = '{}{}'.format(
            file_mature_ref_this_species_tmp, extension0)

    if not re.search('none', file_mature_ref_other_species, re.IGNORECASE):
        # copy file
        (file_mature_ref_other_species_tmp, path0, extension0) = fileparse(
            file_mature_ref_other_species, '\..*')
        cmd = "rna2dna.py {} > {}/{}{}\n\n".format(
            file_mature_ref_other_species,
            dir_tmp,
            file_mature_ref_other_species_tmp,
            extension0
        )
        print_stderr(cmd)

        # here give file name
        ret_parse_mature_ref_other_species = os.popen(cmd).read()
        # rename orig file
        file_mature_ref_other_species = '{}{}'.format(
            file_mature_ref_other_species_tmp, extension0)
        end()

    if not re.search('none', file_precursors, re.IGNORECASE):
        # copy file
        (file_precursors_tmp, path0, extension0) = fileparse(
            file_precursors, '\..*')
        cmd = "rna2dna.py {} > {}/{}{}\n\n".format(
            file_precursors,
            dir_tmp,
            file_precursors_tmp,
            extension0
        )
        print_stderr(cmd)
        # here give file name
        ret_parse_precursors = os.popen(cmd).read()
        # rename orig file
        file_precursors = '{}{}'.format(file_precursors_tmp, extension0)
        end()

    return 0


def parse_mappings():
    global file_reads_vs_genome, parsed_arf, dir_tmp
    # parse mappings to retain only perfect mappings of reads 18 nt <= length
    # <= 25 nt that map perfectly to five loci or less
    pprint("#parsing genome mappings\n")
    print_stderr("#parsing genome mappings\n")

    cmd = "parse_mappings.py {} -a 0 -b 18 -c 25 -i 5 > {}/{}_parsed.arf\n\n".format(
        file_reads_vs_genome,
        dir_tmp,
        parsed_arf
    )
    print_stderr(cmd)

    start()
    ret_parse_mappings = os.popen(cmd).read()
    end()

    return 0


def prepare_signature():
    '''
    prepare signature file
    '''
    global file_reads, dir_tmp, read_align_mismatches, file_mature_ref_this_species, ltime
    pprint("#preparing signature\n")
    print_stderr("#preparing signature\n")

    if not re.search('none', file_mature_ref_this_species, re.IGNORECASE):
        cmd = "prepare_signature.py {} {}/precursors.fa {} -a {}/{} -o {}/signature.arf 2>>error_{}.log\n\n".format(
            file_reads,
            dir_tmp,
            read_align_mismatches,
            dir_tmp,
            file_mature_ref_this_species,
            dir_tmp,
            ltime
        )
        print_stderr(cmd)
        start()
        ret_prepare_signature = os.popen(cmd).read()
        end()
    else:
        cmd = "prepare_signature.py {} {}/precursors.fa {} -o {}/signature.arf 2>>error_{}.log\n\n".format(
            file_reads,
            dir_tmp,
            read_align_mismatches,
            dir_tmp,
            ltime
        )
        start()
        ret_prepare_signature = os.popen(cmd).read()
        end()

    return 0


def excise_precursors():
    global file_genome, parsed_arf, dir_tmp, stack_height_min, dir_tmp, max_pres
    # excise precursors from the genome
    pprint("#excising precursors\n")
    print_stderr("#excising precursors\n")

    start()
    ret_excise_precursors = None

    if options.get('-a'):
        cmd = "excise_precursors.py {} {}/{}_parsed.arf {}/precursors.coords -a {} > {}/precursors.fa\n\n".format(
            file_genome,
            dir_tmp,
            parsed_arf,
            dir_tmp,
            stack_height_min,
            dir_tmp
        )
        print_stderr(cmd)
        ret_excise_precursors = os.popen(cmd).read()
    else:
        cmd = "excise_precursors_iterative_final.py {} {}/{}_parsed.arf {}/precursors.fa {}/precursors.coords {}\n".format(
            file_genome,
            dir_tmp,
            parsed_arf,
            dir_tmp,
            dir_tmp,
            max_pres
        )
        print_stderr(cmd)
        ret_excise_precursors = os.popen(cmd).read()

        fname = '{}/precursors.fa_stack'.format(dir_tmp)
        OSS = open_or_die2(fname, 'rb')
        stack_height_min = OSS.readline().strip()
        OSS.close()

    end()

    fname = '{}/precursors.fa'.format(dir_tmp)
    # if (-z "$dir_tmp/precursors.fa" or not -f "$dir_tmp/precursors.fa"):
    if not file_s(fname) or not os.path.isfile(fname):  # empty or not a regular plain file
        die("No precursors excised\n")

    return 0


def fold_precursors():
    '''
    predicting RNA secondary structures with RNAfold
    '''
    global dir_tmp, ltime
    pprint("#folding precursors\n")
    print_stderr("#folding precursors\n")
    print_stderr(
        "RNAfold < {}/precursors.fa -noPS > {}/precursors.str\n\n".format(dir_tmp, dir_tmp))
    start()
    ret_fold_precursors = os.system(
        "RNAfold < {}/precursors.fa -noPS > {}/precursors.str 2>>error_{}.log".format(dir_tmp, dir_tmp, ltime))
    if ret_fold_precursors:
        ret_fold_precursors = os.system(
            "RNAfold < {}/precursors.fa --noPS > {}/precursors.str".format(dir_tmp, dir_tmp))
        if ret_fold_precursors:
            die("Some RNAfold error occurred. Error {}\n".format(ret_fold_precursors))

    end()


def compute_randfold():
    global options, dir_tmp
    if options.get('-c') == '':
        return

    # compute randfold p-values for the subset of precursors which are
    # plausible Dicer substrates

    pprint("#computing randfold p-values\n")
    print_stderr("#computing randfold p-values\n")
    cmd = "select_for_randfold.py {}/signature.arf {}/precursors.str > {}/precursors_for_randfold.ids\n\n".format(
        dir_tmp,
        dir_tmp,
        dir_tmp
    )
    print_stderr(cmd)
    start()
    ret_select_for_randfold = os.system(cmd)
    end()

    start()
    cmd = "fastaselect.py {}/precursors.fa {}/precursors_for_randfold.ids > {}/precursors_for_randfold.fa\n\n".format(
        dir_tmp,
        dir_tmp,
        dir_tmp
    )
    print_stderr(cmd)
    ret_fasta_select = os.system(cmd)
    end()

    start()
    cmd = "randfold -s {}/precursors_for_randfold.fa 99 > {}/precursors_for_randfold.rand\n\n".format(
        dir_tmp,
        dir_tmp
    )
    print_stderr(cmd)
    ret_randfold = os.system(cmd)
    end()


def get_longest_id(f):
    l = 0
    IN = open_or_die(f, 'rb', 'No file given for checking\n')
    while True:
        line = IN.readline()
        if not line:
            break
        m = re.findall(r'>(\S+)', line)
        if m:
            if len(m[0]) > l:
                l = len(m[0])

    IN.close()
    return l


def core_algorithm():
    '''
    run moRNA Finder core algorithm
    '''
    global _dir, dir_tmp, file_mature_ref_other_species, ltime
    pprint("#running moRNA Finder core algorithm\n")
    print_stderr("#running moRNA Finder core algorithm\n")
    line = None

    longest_id = 40
    if not re.search('none', file_mature_ref_this_species, re.IGNORECASE):
        longest_id = get_longest_id(
            "{}/{}".format(dir_tmp, file_mature_ref_this_species))

    start()

    if not re.search('none', file_mature_ref_other_species, re.IGNORECASE):
        line = "core_algorithm.py {}/signature.arf {}/precursors.str -s {}/{} -v -50 -l {}".format(
            dir_tmp,
            dir_tmp,
            dir_tmp,
            file_mature_ref_other_species,
            longest_id
        )
    else:
        line = "core_algorithm.py {}/signature.arf {}/precursors.str -v -50 -l {}".format(
            dir_tmp,
            dir_tmp,
            longest_id
        )

    if not options.get('-c') == '':
        line += " -y {}/precursors_for_randfold.rand".format(dir_tmp)

    cmd = "{} > {}/output.mrd\n".format(line, _dir)
    print_stderr(cmd)
    ret_mor_core = os.system(cmd)
    if options.get('-E'):
        ret_mor_core = os.system(
            '{} -t > {}/error.output.mrd'.format(line, _dir))

    end()

    # check if file is empty
    fname = "{}/output.mrd".format(_dir)
    if not file_s(fname):
        print_stderr("Error:\n\tFile {} is empty\n\n".format(fname))
        print_stderr(
            "Now running core_algorithm.py with option -t to see why all precursors were discarded\n")
        ret_mor_core = os.system(
            '{} -t > error.output.mrd_{}'.format(line, ltime))
        print_stderr(
            "The debug file is called error.output.mrd_{}\n".format(ltime))
        die("\nExiting now\n\n")


def perform_controls():
    global dir_tmp, _dir, file_mature_ref_other_species, ltime
    # run permuted controls:
    pprint("#running permuted controls\n")
    print_stderr("#running permuted controls\n")
    start()
    line = None

    if not re.search('none', file_mature_ref_other_species, re.IGNORECASE):
        line = "core_algorithm.py {}/signature.arf {}/precursors.str -s {}/{} -v -50".format(
            dir_tmp,
            dir_tmp,
            dir_tmp,
            file_mature_ref_other_species
        )
    else:
        line = "core_algorithm.py {}/signature.arf {}/precursors.str -v -50".format(
            dir_tmp,
            dir_tmp
        )

    if not(options.get('-c') == ''):
        line += " -y {}/precursors_for_randfold.rand".format(dir_tmp)

    cmd = "echo '{} > {}/output.mrd' > {}/command_line\n\n".format(
        line,
        _dir,
        dir_tmp,
    )
    print_stderr(cmd)
    ret_command_line = os.system(cmd)
    cmd = "perform_controls.py {}/command_line {}/precursors.str 100 -a > {}/output_permuted.mrd 2>>error_{}.log\n\n".format(
        dir_tmp,
        dir_tmp,
        dir_tmp,
        ltime
    )
    print_stderr(cmd)
    ret_perform_controls = os.system(cmd)
    end()


def make_survey():
    # get overview of the output:
    global _dir, dir_tmp, file_mature_ref_this_species, stack_height_min
    pprint("#doing survey of accuracy\n")
    print_stderr("#doing survey of accuracy\n")

    if not re.search('none', file_mature_ref_this_species, re.IGNORECASE):
        cmd = "survey.py {}/output.mrd -a {}/output_permuted.mrd -b {}/{} -c {}/signature.arf -d {} > {}/survey.csv\n\n".format(
            _dir,
            dir_tmp,
            dir_tmp,
            file_mature_ref_this_species,
            dir_tmp,
            stack_height_min,
            _dir
        )
        print_stderr(cmd)
        start()
        ret_survey = os.system(cmd)
        end()

    else:

        cmd = "survey.py {}/output.mrd -a {}/output_permuted.mrd -d {} > {}/survey.csv\n\n".format(
            _dir,
            dir_tmp,
            stack_height_min,
            _dir
        )
        print_stderr(cmd)
        start()
        ret_survey = os.system(cmd)
        end()


def output_results():
    '''
    making final results html file:
    '''
    global options, dir_tmp, ltime, version, scripts, file_mature_ref_this_species
    pprint("#producing graphic results\n")
    print_stderr("#producing graphic results\n")
    start()

    # sort aligned reads in pdf not by sample if option -o is given
    sort_by_sample = '-o'
    if options.get('-o'):
        sort_by_sample = ''

    line = None

    # choose file to use for counting miRNAs in data
    xopt = "{}/signature.arf".format(dir_tmp)
    if os.path.isfile("expression_analyses/expression_analyses_{}/miRNA_expressed.csv".format(ltime)):
        xopt = "expression_analyses/expression_analyses_{}/miRNA_expressed.csv".format(
            ltime)

    sc = 0
    if options.get('-b'):
        sc = options.get('-b')

    OE = ""
    if options.get('-E') == '':
        OE = " -E"

    if not re.search('none', file_mature_ref_this_species, re.IGNORECASE):
        if options.get('-q'):
            line = "make_html.py -f {}/output.mrd -k {}/{} -p {}/precursors.coords -s {}/survey.csv -c -e -q {} -x {} -r {}Rfam_for_moR.fa -v {} -y {} {} {}".format(
                _dir,
                dir_tmp,
                file_mature_ref_this_species,
                dir_tmp,
                _dir,
                options.get('-q'),
                xopt,
                scripts,
                sc,
                ltime,
                sort_by_sample,
                OE
            )
        else:
            line = "make_html.py -f {}/output.mrd -k {}/{} -p {}/precursors.coords -s {}/survey.csv -c -e -r {}Rfam_for_moR.fa -v {} -y {}  {} {}".format(
                _dir,
                dir_tmp,
                file_mature_ref_this_species,
                dir_tmp,
                _dir,
                scripts,
                sc,
                ltime,
                sort_by_sample,
                OE
            )
    else:
        if options.get('-q'):
            line = "make_html.py -f {}/output.mrd -p {}/precursors.coords -s {}/survey.csv -c -e -q {}  -x {} -r {}Rfam_for_moR.fa -v {} -y {} {} {}".format(
                _dir,
                dir_tmp,
                _dir,
                options.get('-q'),
                xopt,
                scripts,
                sc,
                ltime,
                sort_by_sample,
                OE
            )
        else:
            line = "make_html.py -f {}/output.mrd -p {}/precursors.coords -v {} -s {}/survey.csv -c -e -r {}Rfam_for_moR.fa -y {} {} {}".format(
                _dir,
                dir_tmp,
                sc,
                _dir,
                scripts,
                ltime,
                sort_by_sample,
                OE
            )

    if options.get('-t'):
        line += " -t {}".format(options.get('-t'))

    dopt = ""
    if options.get('-d'):
        dopt = "-d"

    cmd = '{} -V {} {}\n\n'.format(line, version, dopt)
    print_stderr(cmd)
    ret_make_html = os.system(cmd)

    end()


def remove_dir_tmp():
    global dir_tmp, options
    # remove temporary directory
    if options.get('-v') == '':
        print_stderr("rmtree({})\n\n".format(dir_tmp))
        shutil.rmtree(dir_tmp)


def make_bed():
    global ltime
    cmd = 'morbed.py result_{}.csv > result_{}.bed'.format(ltime, ltime)
    print_stderr(cmd, '\n')
    res = os.popen(cmd).read()
    if not res:
        return 0
    else:
        print_stderr(res, "\n")
        return 1


def extract_sequences_from_results():
    global ltime
    od = "mirna_results_{}".format(ltime)

    cmd = 'get_precursors.py -r result_{}.csv -p -d -o {}'.format(
        ltime, od)
    print_stderr(cmd, '\n')
    res = os.system(cmd)
    cmd = 'get_precursors.py -r result_{}.csv -m -p -d -o {}'.format(
        ltime, od)
    print_stderr(cmd, '\n')
    res1 = os.system(cmd)
    cmd = 'get_precursors.py -r result_{}.csv -k -p -d -o {}'.format(
        ltime, od)
    print_stderr(cmd, '\n')
    res2 = os.system(cmd)

    if not res and not res1 and not res2:
        print_stderr(
            "fasta and bed files have been created in subfolder {}\n".format(od))
        return 0
    else:
        print_stderr(res, "\n")
        print_stderr(res1, "\n")
        print_stderr(res2, "\n")
        return 1


if __name__ == '__main__':
    # measuring times
    (sTime, eTime, stime, etime, sTimeG, eTimeG, stimeG, etimeG) = (
        None, None, None, None, None, None, None, None)
    (second, minute, hour, dayOfMonth, month, yearOffset, dayOfWeek, dayOfYear,
     daylightSavings) = (None, None, None, None, None, None, None, None, None)
    (second, minute, hour, dayOfMonth, month, yearOffset,
     dayOfWeek, dayOfYear, daylightSavings) = localtime()

    if re.search(r'^\d$', str(second)):
        second = "0{}".format(second)

    sTimeG = "{}:{}:{}".format(hour, minute, second)
    stimeG = int(time.time())

    pprint("moRNA Finder started at {}\n\n\n".format(sTimeG))

    ctime = int(time.time())
    ltime = myTime(ctime)

    scripts = os.popen('which moR.py').read()
    # scripts = re.sub(r'moR.py', '', scripts, count=1)
    # scripts = re.sub(r'\s+', '', scripts)
    scripts = os.path.dirname(scripts) + '/'

    pprint('#Starting moRNA Finder\n')
    print_stderr('#Starting moRNA Finder\n{} {}\n\n'.format(
        sys.argv[0], ' '.join(sys.argv[1:])))
    print_stderr("moRNA Finder started at {}\n\n\n".format(sTimeG))

    test_first_argument()

    command_line = "{} {}\n".format(sys.argv[0], ' '.join(sys.argv[1:]))

    if len(sys.argv) < 7:
        die(usage)

    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument('reads', help='reads')
    parser.add_argument('genome', help='genome')
    parser.add_argument('mappings', help='mappings')
    parser.add_argument('miRNAs_ref', help='miRNAs_ref')
    parser.add_argument('miRNAs_other', help='miRNAs_other')
    parser.add_argument('precursors', help='precursors')
    args = parser.parse_args(sys.argv[1:7])
    file_reads = args.reads
    file_genome = args.genome
    file_reads_vs_genome = args.mappings
    file_mature_ref_this_species = args.miRNAs_ref
    file_mature_ref_other_species = args.miRNAs_other
    file_precursors = args.precursors

    parsed_arf = None
    m = re.findall(r'\/*([a-zA-Z_0-9\.]+)$', file_reads_vs_genome)
    if m:
        parsed_arf = m[0]
    else:
        die("could not match arf file\n")

    if os.path.exists(file_mature_ref_this_species) or file_mature_ref_this_species == 'none':
        pass
    else:
        die(usage, '\n\nError: no file containing miRNAs of investigating species specified\n\nEither specify a valid fasta file or say none\n\n', warning)

    if os.path.exists(file_mature_ref_other_species) or file_mature_ref_other_species == 'none':
        pass
    else:
        die(usage, '\n\nError: no file containing miRNAs of other species specified\n\nEither specify a valid fasta file or say none\n\n', warning)

    if os.path.exists(file_precursors) or file_precursors == 'none':
        pass
    else:
        die(usage, '\n\nError: no file containing miRNAs precursors of investigating species specified\n\nEither specify a valid fasta file or say none\n\n', warning)

    opts, argss = getopt.getopt(sys.argv[7:], 'a:b:cdt:uvq:s:z:r:p:g:EP')
    options = dict(opts)

    max_pres = 50000
    if options.get('-g'):
        max_pres = options.get('-g')

    # minimal precursor length, used for precheck of precursor file
    minpreslen = 40
    if options.get('-p'):
        minpreslen = options.get('-p')

    stack_height_min = None
    if options.get('-a'):
        stack_height_min = options.get('-a')

    if options.get('-z'):
        time += options.get('-z')

    make_dir_tmp()

    test_installed_binaries()

    printUsedParameters(options)

    test_input_files()

    rna2dna()  # this makes u to t, lc to uc in sequence and skips everything behind first whitespace in fasta identifier

    parse_mappings()

    excise_precursors()

    prepare_signature()

    fold_precursors()

    compute_randfold()

    core_algorithm()

    perform_controls()

    make_survey()

    output_results()

    make_bed()

    extract_sequences_from_results()

    remove_dir_tmp()

    etimeG = int(time.time()) - stimeG
    (second, minute, hour, dayOfMonth, month, yearOffset,
     dayOfWeek, dayOfYear, daylightSavings) = localtime()
    if re.search(r'^\d$', str(second)):
        second = "0{}".format(second)
    eTimeG = "{}:{}:{}".format(hour, minute, second)

    message = '\nmoRNA Finder runtime: \n\nstarted: {}\nended: {}\ntotal:{}h:{}m:{}s\n\n'.format(
        sTimeG,
        eTimeG,
        int(etimeG / 3600),
        int((etimeG % 3600) / 60),
        int(etimeG % 60), "s\n\n"
    )
    pprint(message)
    print_stderr(message)

    fname = '{}/run_{}_parameters'.format(_dir, ltime)
    OUT = open_or_die(
        fname, 'a', 'parameter file {} not found\n'.format(fname))
    OUT.write(message)
    OUT.close()

    sys.exit(0)
