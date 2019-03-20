#!/usr/bin/env python
from __future__ import print_function

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
from reportlab.lib.colors import (black, darkviolet, green, grey, lightskyblue,
                                  orange, red)
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

__DATA__ = '''Human
Chimp
Orangutan
Rhesus
Marmoset
Mouse
Rat
Guinea Pig
Cat
Dog
Horse
Cow
Opossum
Platypus
Chicken
Zebra finch
Lizard
X. tropicalis
Zebrafish
Tetraodon
Fugu
Stickleback
Medaka
Lamprey
Lancelet
C. intestinalis
S. purpuratus
C. elegans
C. brenneri
C. briggsae
C. remanei
C. japonica
P. pacificus
D. melanogaster
D. simulans
D. sechellia
D. yakuba
D. erecta
D. ananassae
D. pseudoobscura
D. persimilis
D. virilis
D. mojavensis
D. grimshawi
A. gambiae
A. mellifera
S. cerevisiae'''

options = {}

n_count = 0
k_count = 0
e_count = 0
sig = 0
infile = None   # moRNA Finder output file
pdfs = 0        # force pdf creation
threshold = 0    # hairpins have to have score above threshold in order to be reported
csv = 0
xposshift = 40
lstruct_multi = None
sb_obs = None
known = None
organisms = {}
# read in moRNA Finder output file
_id = None
_hash = {}   # takes up all entries from the moRNA Finder module
seen = {}
created = 1
_in = None
counter = 0
hstruct = {}  # counts how often a nt is covered reads
i = None
offset = 0
me = 0     # mature end coordinate
desc = []
lflank1 = 0  # length(string of left flank)
fl1 = None
lflank2 = None  # length string of right flank
fl2b = 0        # right flank begin
lloop = None    # string of loop
lb = 0          # starting position of loop
lstar = None    # string of star sequence
sb = 0          # starting
lmature = None  # string of mature
mb = 0          # mature begin
sstruct = ''    # structure string
pri_seq = None  # pri-cursor sequence
lenstr = 0
pdf = None            # pdf descriptor
page = None           # page descriptor
gfx = None            # graphic variable
trb = None            # fontvariable
text = None
text2 = None
aligned = None                        # reads reading
hash2 = {}                            # begin of read in precursor
hash2c = {}                           # number of reads per read sequence
hash2key = {}
hash2mm = {}                          # number of mismatches
hash2order = {}                       # output order saved
hash2seq = {}
hash2sample = {}
# stores begin coordinates of fl1,m,l,s,fl2 sequences
order = {}
multiplier = 3.6  # 4.825                # minimal distance between two letters
# calculate predefined pdf loci for alignment characters
position_hash = {}
counter = 0
yorig = 500  # 500
downy = 50
dline = None                              # line graphic handler
first = 1
lastx = None
lasty = None
final = None                            # final output string of a read
pseq = []                              # precursor sequence
rseq = []                              # read sequence
totalreads = 0
# color assigned to position where letter is drawn
assign_str = {}
assign_str_exp = {}
bpo1 = -10                              # left nt pos in first bp
bpo2 = -10                              # right nt pos in first bp
bpo1r = -10                             # left nt pos in second bp
bpo2r = -10                             # right nt pos in second bp
ffe = 0                                 # first flank end position
ff2b = 0                                # second flank begin position
# array that stores sorted order of fl1,m,l,s,fl2
_sorted = []
y = yorig                               # y coordinate

# min and max x,y coordinates of rna sequence
minx, miny, maxx, maxy = (None, None, None, None)
rna = []                  # rna sequence
rna_d = []
xc = {}                   # holds x cooridnate of each nt
yc = {}                   # holds y coordinate of each nt
sid = ""

pres_coords = {}

# pdf histogram colors
col_star_exp = lightskyblue
col_star_obs = darkviolet
col_mature = red
col_loop = orange
hm = {}
hs = {}
hp = {}
weighted = {}
files_mirnaex = []
ltime = None

# determine pdf path when running on a cluster
pdf_path = None
cwd = None
org = None

# some quantifier variables
mirbase = 0
mature2hairpin = {}
# some hairpins have more than 1 mature assigned, circumvent this problem
hairpin2mature = {}
hash_q = {}           # takes up all entries from the quantifier module
blast = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY="
blast_query = "&db=nucleotide&QUERY_FROM=&QUERY_TO=&QUERYFILE=&GENETIC_CODE=1&SUBJECTS=&stype=nucleotide&SUBJECTS_FROM=&SUBJECTS_TO=&SUBJECTFILE=&DBTYPE=gc&DATABASE=nr&EQ_MENU=&NUM_ORG=1&EQ_TEXT=&BLAST_PROGRAMS=blastn&PHI_PATTERN=&MAX_NUM_SEQ=100&SHORT_QUERY_ADJUST=on&EXPECT=10&WORD_SIZE=7&MATRIX_NAME=PAM30&MATCH_SCORES=2,-3&GAPCOSTS=5+2&COMPOSITION_BASED_STATISTICS=0&FILTER=L&REPEATS=repeat_9606&FILTER=m&TEMPLATE_LENGTH=0&TEMPLATE_TYPE=0&PSSM=&I_THRESH=&SHOW_OVERVIEW=true&SHOW_LINKOUT=true&GET_SEQUENCE=auauauaauauauauauauuauaa&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&ALIGNMENT_VIEW=Pairwise&MASK_CHAR=2&MASK_COLOR=1&DESCRIPTIONS=100&ALIGNMENTS=100&NEW_VIEW=true&OLD_BLAST=false&NCBI_GI=false&SHOW_CDS_FEATURE=false&NUM_OVERVIEW=100&FORMAT_EQ_TEXT=&FORMAT_ORGANISM=&EXPECT_LOW=&EXPECT_HIGH=&QUERY_INDEX=&CLIENT=web&SERVICE=plain&CMD=request&PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=&RUN_PSIBLAST=&TWO_HITS=&DEFAULT_PROG=megaBlast&WWW_BLAST_TYPE=&DB_ABBR=&SAVED_PSSM=&SELECTED_PROG_TYPE=blastn&SAVED_SEARCH=true&BLAST_SPEC=&QUERY_BELIEVE_DEFLINE=&DB_DIR_PREFIX=&USER_DATABASE=&USER_WORD_SIZE=&USER_MATCH_SCORES=&USER_FORMAT_DEFAULTS=&NO_COMMON=&NUM_DIFFS=2&NUM_OPTS_DIFFS=1&UNIQ_DEFAULTS_NAME=A_SearchDefaults_1Mn7ZD_2Sq4_1Z58HQ5Jb_23tpbD_167y9p&PAGE_TYPE=BlastSearch&USER_DEFAULT_PROG_TYPE=blastn&USER_DEFAULT_MATCH_SCORES=3."

# get mature positions if options a is set
mature_pos_hash = {}
confident = {}

knownones = {}
starknownones = {}
seen = {}

CSV = None
HTML = None


def current_dir():
    '''
    Get current dir path
    '''
    return os.path.realpath('.')


def Usage():
    print_stderr(
        "\n\n\n\n[usage]\n\python make_html.py -f moR_outfile [options]\n\n")
    print_stderr("[options]\n-v 2\t only output hairpins with score above 2\n")
    print_stderr("-c\t also create overview in excel format.\n")
    print_stderr("-k file\t supply file with known miRNAs\n")
    print_stderr("-s file\t supply survey file if score cutoff is used to get information about how big is the confidence of resulting reads\n\n\n")
    print_stderr("-e \t report complete survey file\n")
    print_stderr("-g \t report survey for current score cutoff\n")
#    print_stderr("-w project_folder\t automatically used when running webinterface, otherwise don't use it\n")
    print_stderr(
        "-r file\t Rfam file to check for already reported small RNA sequences\n\n")
    print_stderr("-q file\t miRBase.mrd file produced by quantifier module\n")
    print_stderr(
        "-x file\t signature.arf file with mapped reads to precursors\n")
    print_stderr(
        "-t spec\t specify the organism from which your sequencing data was obtained\n")
    print_stderr("-u \t print all available UCSC input organisms\n")
    print_stderr("-y \ttimestamp of this run\n")
    print_stderr(
        "-o \tsort signature by sample in pdf file, default is by beginning position\n")
    print_stderr("-d \t do not generate pdfs\n\n\n")
    print_stderr(
        "-a \tprint genomic coordinates of mature sequence (still testing)\n")
    print_stderr("-b \tsupply confidence file\n")
    print_stderr("-V \tmoRNA Finder version used\n\n")
#   print_stderr("-E \t not fully implemented yet: prints error messages for known precursors that have not been scored by moRNA Finder\n")

    sys.exit(0)


def CreateHTML(HTML):
    data = '''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
    "http://www.w3.org/TR/html4/strict.dtd">
    <html>
    <head>
    <title>moRNA Finder</title>
    <!-- CSS code -->
    <style type="text/css">
    body{
      font: .75em Arial,sans-serif;
      background: #FFFFFF;
      color: #333333
    }
    div#
    container
    {
      width: 1000px;
      margin:0 auto
      }
    h1{
      color: #F60;
      margin: 1em 0 0;
        letter-spacing: -2px;
    }
    p{
      margin: 0 0 1.7em;
    }

    a.tooltip{
        text-decoration:none;
        color:black;
            font-weight:bold;
    }

    a.tooltip:hover
    {    position: relative;
         background: transparent;

    }

    a.tooltip span
    {    position: absolute;
         visibility: hidden;

         width: 20em;
         top: 1.5em;
         background: #ffffdd;
    }

    a.tooltip:hover span
    {    visibility: visible;   }


    a.tooltip2{
        text-decoration:none;
        color:black;
            font-weight:bold;
    }

    a.tooltip2:hover
    {    position: relative;
         background: transparent;

    }

    a.tooltip2 span
    {    position: absolute;
         visibility: hidden;

         width: 20em;
         top: 1.5em;left:-10em;
         background: #ffffdd;
    }

    a.tooltip2:hover span
    {    visibility: visible;   }


    a.tooltip3{
        text-decoration:none;
        color:black;
            font-weight:bold;
    }

    a.tooltip3:hover
    {    position: relative;
         background: transparent;

    }

    a.tooltip3 span
    {    position: absolute;
         visibility: hidden;

         width: 20em;
         top: 3em;
         background: #ffffdd;
    }

    a.tooltip3:hover span
    {    visibility: visible;   }


    a.tooltip4{
        text-decoration:none;
        color:black;
            font-weight:bold;
    }

    a.tooltip4:hover
    {    position: relative;
         background: transparent;

    }

    a.tooltip4 span
    {    position: absolute;
         visibility: hidden;

         width: 20em;
         top:3em;left:-10em;
         background: #ffffdd;
    }

    a.tooltip4:hover span
    {    visibility: visible;   }



    </style>

        </head>
        <body>
        <table border="0" width="100%">
        <colgroup>
        <col width="5*">
        <col width="5*">
        </colgroup>
        <tr height="200" valign="top">
        <td><font face="Times New Roman" size="8">
        <b><a href="http://www.mirbase.org" target="_blank" style="color:Black;text-decoration:none" title="moRNA Finder homepage" >moRNA Finder</a></b>
        <br>
        <br>
        </font></td>
        <td> <a href="http://www.mirbase.org" target="_blank" >
        <img src="http://www.mirbase.org/images/mirbase-logo-blue-web.png" style="border-style: none" name="precursor miRNA" title="moRNA Finder homepage" align=right alt="moRNA Finder home"/></a></td>
        </tr>
        </table>'''
    HTML.write(data)


def get_precursor_pos():
    global options
    IN = open_or_die(options.get('-p'), 'rb',
                     "no precursor.coords file found with argument -p\n")
    _id = None
    while True:
        l = IN.readline()
        if not l:
            break

        m = re.match(r'^>*((\S+)_\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s*$', l)
        if m:
            m = m.groups()
            _id = tr(m[0], '|', '_')
            create_hash_key_chain(pres_coords, '', _id, 'chr')
            pres_coords[_id]['chr'] = m[1]
            create_hash_key_chain(pres_coords, '', _id, 'strand')
            pres_coords[_id]['strand'] = m[2]
            create_hash_key_chain(pres_coords, '', _id, 's')
            pres_coords[_id]['s'] = int(m[3])
            create_hash_key_chain(pres_coords, '', _id, 'e')
            pres_coords[_id]['e'] = int(m[4])

    IN.close()


def get_mature_pos():
    '''
    This function is defined in make_html.py
    '''
    global options, ltime
    if not options.get('-f'):
        die('no input file given for option -f\n')

    IN = open_or_die2(options.get('-f'), 'rb')

    mhash = {}
    mid = None
    start = -1
    seq = None
    (pos, rpos, pos_last, end, read) = (None, None, None, None, None)

    while True:
        line = IN.readline().strip()
        if not line:
            break

        m = re.search(r'>(\S+)')
        if m:
            mid = m.groups()[0]

        m = re.search(r'exp\s+(\S+)', line)
        if m:
            seq = m.groups()[0]
            start = seq.find('M')
            end = seq.rfind('M')

        m = re.match(r'^(\S+_x\d+)\s+(\S+)', line)
        if m:
            m = m.groups()
            read = m[0]
            seq = m[1]
            seq = tr(seq, 'acgtuACGTU', 'AAAAAAAAAA')
            pos = seq.find('A')
            rpos = seq.rfind('A')
            if pos == start and rpos == end and not hash_defined(mhash[read]):
                mhash[mid] = read

    IN.close()

    fname = "moR_runs/run_{}/tmp/signature.arf".format(ltime)
    IN = open_or_die(fname, 'rb', 'no signature file given for parsing\n')
    line = []

    while True:
        l = IN.readline().strip()
        if not l:
            break

        line = esplit(l)
        if hash_value_true(mhash[line[5]]) and mhash[line[5]] == line[0]:
            mature_pos_hash[line[5]]['s'] = line[7]
            mature_pos_hash[line[5]]['e'] = line[8]
            mature_pos_hash[line[5]]['strand'] = line[10]

    IN.close()


def PrintKnownnotfound():
    global CSV, csv, options, ltime, mature2hairpin, hairpin2mature, hash_q, counter, reads, _hash, seen
    not_seen = {}
    signature = {}
    error = {}
    (mature, path0, extension0) = fileparse(options.get('-k'), '\..*')
    # change .fa suffix to _mapped.arf
    mature += "{}_mapped.arf".format(extension0)
    exprs = {}

    IN = open_or_die2(
        'expression_analyses/expression_analyses_{}/miRNA_expressed.csv'.format(ltime), 'rb')
    while True:
        l = IN.readline().strip()
        if not l:
            break

        if re.search(r'precursor', l):
            continue

        line = re.split(r'\t', l)

        # here comes up the trouble when two mature map to same precursor
        mature2hairpin[line[0]] = line[2]
        hairpin2mature[line[2]] = line[0]
        exprs[line[0]] = line[1]

    IN.close()

    # read in error messages of moRNA Finder run
    if options.get('-E'):
        infile = options.get('-s')
        infile = re.sub('survey.csv', 'error.output.mrd', infile, count=1)
        try:
            IN = open(infile, 'rb')
            _in = 0
            _id = None
            mess = ''
            while True:
                l = IN.readline()
                if not l:
                    break

                m = re.search(r'>(\S+)', l)
                if m:
                    _in = 1
                    _id = m.groups()[0]
                    mess = ''
                    continue

                m = re.match(r'^exp', l)
                if m:
                    _in = 0
                    continue

                if _in:
                    mess += l
                else:
                    if re.search(r'pri_seq', l):
                        continue
                    if re.search(r'pri_struct', l):
                        continue
                    m = re.match(r'^(\S+)', l)
                    if m:
                        m = m.groups()
                        if hash_value_true(knownones, m[0]) and not hash_value_true(seen, m[0]):
                            create_hash_key_chain(error, m[0], _id)
                            error[m[0]][_id] += mess

            IN.close()
        except IOError:
            print_stderr(infile, ' not found\n')

    # only performed when option z is NOT given
    if not options.get('-z') == '':
        # this are all reads of known miRBase miRNAs seen in output.mrd

        # put all not seen in output.mrd miRNAs in new hash %not_seen
        for key in knownones.keys():   # all miRNAs in mature_ref_this_species listed things
            if not(key in seen.keys() and seen[key]):
                not_seen[key] = 1

        # read in all micrornas in data
        if options.get('-x'):
            IN = open_or_die2(options.get('-x'), 'rb')
            mmu = []
            while True:
                l = IN.readline()
                if not l:
                    break

                if re.search(r'ID\s+read count', l):
                    continue

                mmu = re.split(r'\t', l)
                if mmu[0] in knownones.keys() and knownones[mmu[0]] or (mmu[0] in mature2hairpin.keys() and mature2hairpin[mmu[0]] in knownones.keys() and knownones[mature2hairpin[mmu[0]]]):
                    signature[mmu[0]] = 1

            IN.close()

    # open miRBase.mrd from quantifier module ## precursor ids as hash
    IN = open_or_die2(options.get('-q'), 'rb')
    while True:
        l = IN.readline()
        if not l:
            break

        m = re.match(r'^\>(\S+)', l)
        if m:
            m = m.groups()
            _id = m[0]
            oid = m[0]
            _id = tr(_id, '|', '_')

            create_hash_key_chain(hash_q, '', _id, 'oid')
            hash_q[_id]["oid"] = oid
            create_hash_key_chain(hash_q, '', _id, 'id')
            hash_q[_id]["id"] = _id
            counter = 0
        else:
            m = re.match(r'^remaining read count\s*(\d*)', l)
            if m:
                m = m.groups()
                create_hash_key_chain(hash_q, '', _id, 'remaining_loop')
                hash_q[_id]["remaining_loop"] = m[0]
            else:
                m = re.match(r'^total read count\s*(\d*)', l)
                if m:
                    m = m.groups()
                    create_hash_key_chain(hash_q, '', _id, 'freq_total')
                    hash_q[_id]["freq_total"] = m[0]
                else:
                    m = re.match(r'^(\S+) read count\s*(\d*)', l)
                    if m:
                        hash_q[_id][m.groups()[0]] = m.groups()[1]
                    else:
                        m = re.match(r'^pri_seq\s+(\S+)', l)
                        if m:
                            m = m.groups()
                            create_hash_key_chain(hash_q, '', _id, 'pri_seq')
                            hash_q[_id]["pri_seq"] = m[0]
                            d = []
                            if hash_value_true(hash_q, _id, 'obs'):
                                d = ssplit(hash_q[_id]['obs'])
                            else:
                                d = ssplit(hash_q[_id]['exp'])

                            # put precursor sequence to array
                            s = ssplit(hash_q[_id]["pri_seq"])
                            mseq = ""
                            sseq = ""

                            # now set labels for mature and star seq in
                            # explanation string
                            for i in range(0, len(m[0])):
                                if d[i] != "f":
                                    create_hash_key_chain(
                                        hash_q, '', _id, 'ucsc_seq')
                                    hash_q[_id]["ucsc_seq"] += s[i]
                                # accoutn for everything else
                                if d[i] == "M" or d[i] == "5" or d[i] == "3":
                                    mseq += s[i]
                                elif d[i] == "S":
                                    sseq += s[i]

                            # if there is an observed star sequence use this
                            sseq_obs = ""

                            defined = False
                            try:
                                hash_q[_id]['obs']
                                defined = True
                            except KeyError:
                                pass

                            if defined and hash_q[_id]['obs']:
                                d = ssplit(hash_q[_id]['obs'])
                                for i in range(0, len(m[0])):
                                    if d[i] == "S":
                                        sseq_obs += s[i]

                            create_hash_key_chain(hash_q, '', _id, 'mat_seq')
                            hash_q[_id]["mat_seq"] = mseq
                            create_hash_key_chain(hash_q, '', _id, 'star_seq')
                            hash_q[_id]["star_seq"] = sseq
                            create_hash_key_chain(
                                hash_q, '', _id, 'star_seq_obs')
                            hash_q[_id]["star_seq_obs"] = sseq_obs
                        else:
                            m = re.match(r'^exp\s+(\S+)', l)
                            if m:
                                m = m.groups()
                                create_hash_key_chain(hash_q, '', _id, 'exp')
                                hash_q[_id]['exp'] = m[0]
                            else:
                                m = re.match(r'^obs\s+(\S+)', l)
                                if m:
                                    m = m.groups()
                                    create_hash_key_chain(
                                        hash_q, '', _id, 'obs')
                                    hash_q[_id]['obs'] = m[0]
                                else:
                                    m = re.match(r'^pri_struct\s+(\S+)', l)
                                    if m:
                                        m = m.groups()
                                        hash_q[_id]["pri_struct"] = m[0]
                                        reads = 1
                                        counter = 0
                                        continue
                                    else:
                                        m = re.match(
                                            r'^(\S\S\S)(\S+)\s+(\S+)\s+(\S+)$', l)
                                        if m:
                                            m = m.groups()
                                            counter += 1
                                            create_hash_key_chain(
                                                hash_q, '', _id, 'reads', m[0], counter, 'rid')
                                            hash_q[_id]["reads"][m[0]][counter][
                                                "rid"] = "{}{}".format(m[0], m[1])
                                            create_hash_key_chain(
                                                hash_q, '', _id, 'reads', m[0], counter, 'seq')
                                            hash_q[_id]["reads"][m[0]][
                                                counter]["seq"] = m[2]
                                            create_hash_key_chain(
                                                hash_q, '', _id, 'reads', m[0], counter, 'mm')
                                            hash_q[_id]["reads"][m[0]][
                                                counter]["mm"] = m[3]

    IN.close()
    print_stderr('parsing miRBase.mrd file finished\n')

    # print out not scored but in data miRNAs by moRNA Finder to html.
    _id = None
    sf = None
    mf = None

    # print to csv these miRNAs as well
    if csv:
        fname = '{}/result_{}.csv'.format(cwd, ltime)
        CSV = open_or_die2(fname, 'a+')
        CSV.write('\n#miRBase miRNAs not detected by moRNA Finder\n')

    # if make_html was called by quantifier.py then output all precursors that
    # are expressed
    if re.search('miRBase.mrd', options.get('-f')):
        pass
    else:
        # if processing only the output of a normal moRNA Finder run
        exhtml = {}
        start = 0
        mirc = None

        rt = 0

        fname = 'expression_analyses/expression_analyses_{}/expression_{}.html'.format(
            options.get('-y'), options.get('-y'))
        IN = open_or_die2(fname, 'rb')
        while True:
            l = IN.readline()
            if not l:
                break

            if re.search(r'miRBase precursor id', l):
                start = 1

            if not start:
                continue

            # read in header
            if re.search(r'nowrap', l):
                start = 2

            if re.search(r'^\s+<\/tr>\s*$', l):
                start = 3

            if start == 1:
                create_hash_key_chain(exhtml, '', 'header')
                exhtml['header'] += l
            elif start == 2:
                if re.search(r'pdf', l):
                    m = re.search(r'([a-zA-Z0-9-]+)<\/a><\/td>', l)
                    if m:
                        m = m.groups()
                        mirc = m[0]
                        create_hash_key_chain(exhtml, '', mirc)
                        exhtml[mirc] += l
                else:
                    create_hash_key_chain(exhtml, '', mirc)
                    exhtml[mirc] += l

        IN.close()

        HTML.write('<br><br><h2>mature miRBase miRNAs not detected by moRNA Finder</h2><br><font face=\"Times New Roman\" size=\"2\">\n <table border=\"1\"> {}'.format(
            exhtml['header']
        ))

        csvout = []
        m = re.findall(r'tooltip\d+\">(.+)<span>', exhtml['header'])
        for mm in m:
            csvout.append(mm)

        if csv:
            CSV.write('\t'.join(csvout))
            CSV.write('\n')

        skip = {}
        for oid in hash_sort_key(exprs, lambda x: x[1] * -1):
            if oid in mature2hairpin.keys() and mature2hairpin[oid] in skip.keys() and skip[mature2hairpin[oid]]:
                continue

            skip[mature2hairpin[oid]] = 1
            if ((hash_value_true(not_seen, oid) or hash_value_true(not_seen, mature2hairpin[oid])) and hash_value_true(signature, oid)) or options.get('-z') == '':
                create_hash_key_chain(exhtml, '', mature2hairpin[oid])
                exhtml[mature2hairpin[oid]] = re.sub(
                    r'<td>\s*<\/td>', '<td>-<\/td>', exhtml[mature2hairpin[oid]])
                csvout = re.split('\n', exhtml[mature2hairpin[oid]])

                count = None
                icount = None
                _sum = 0
                _in = None
                notp, line = (None, None)

                for cl in csvout:
                    notp = 0
                    line = cl

                    if re.search(r'>\S+<\/a><\/td>', cl):
                        m = re.search(r'pdf\s*<\/div>(\S+)<\/a><\/td>', cl)
                        if m:
                            if csv:
                                CSV.write(m.groups()[0])
                                CSV.write('\t')

                        if re.search(r'nobr', cl):
                            icount = 0
                            _in = 1
                            m = re.findall(r'<nobr>(\d+)', cl)
                            for mm in m:
                                notp = 1
                                _in += 1
                                if _in % 2 == 1:
                                    continue
                                _sum += int(mm)
                                icount += 1
                                if count == icount:
                                    break

                            m = re.match(r'^(.+<\/a>)<\/td><td><nobr>', line)
                            if m:
                                line = m.groups()[0]
                    elif icount > 0 and re.search(r'<\/table>', cl):
                        if csv:
                            CSV.write(str(_sum))
                            CSV.write('\t')

                        HTML.write(str(_sum))
                        HTML.write('</td>')
                        icount = 0
                        _sum = 0
                    elif re.search('norm', cl):
                        count = 0
                        m = re.findall('norm', cl)
                        for mm in m:
                            count += 1

                        line = substr(line, 0, 7)
                    else:
                        m = re.search(r'>(\d+\.*\d*)<\/td>', cl)
                        if m:
                            if csv:
                                CSV.write(m.groups()[0])
                                CSV.write('\t')
                        else:
                            if re.search(r'<td>\s*([acgtACGTuU]+)\s*<\/td>', cl):
                                m = re.findall(
                                    r'<td>\s*([acgtACGTuU]+)\s*<\/td>', cl)
                                for mm in m:
                                    if csv:
                                        CSV.write(mm)
                                        CSV.write('\t')
                            else:
                                m = re.search(r'"nowrap">(\d+)<\/td>', cl)
                                if m:
                                    if csv:
                                        CSV.write(m.groups()[0])
                                        CSV.write('\t')

                    if notp:
                        pass
                    else:
                        if re.search(r'^\s*-\s*$', cl) and csv:
                            CSV.write('-\t')
                        elif re.search(r'<td>\s*-\s*<\/td>', cl) and csv:
                            CSV.write('-\t')

                        if re.search(r'blast.ncbi', cl) and csv:
                            CSV.write('-\t')

                        HTML.write(line)

                if csv:
                    CSV.write('\n')

    if csv:
        CSV.close()

    return


def check_Rfam(hsh):
    global ltime, _hash, options
    err = None
    if options.get('-q'):
        TMP = open_or_die2(
            'expression_analyses/expression_analyses_{}/identified_precursors.fa'.format(ltime), 'w+')
    else:
        TMP = open_or_die2(
            'moR_runs/run_{}/identified_precursors.fa'.format(ltime), 'w+')

    seqm, seql, seqs = None, None, None
    for key in hsh.keys():
        k = key.strip()
        seqm = hsh[k]['mat_seq']
        seqm = tr(seqm, 'uUacg', 'TTACG')

        seql = hsh[k]['loop_seq']
        seql = tr(seql, 'uUacg', 'TTACG')

        seqs = hsh[k]['star_seq']
        seqs = tr(seqs, 'uUacg', 'TTACG')

        TMP.write('>{}\n{}\n'.format(k, hsh[k]['ucsc_seq']))
        TMP.write(">M:{}\n{}\n".format(k, seqm))

        if len(seql) > 15:
            TMP.write(">L:{}\n{}\n".format(k, seql))

        TMP.write(">S:{}\n{}\n".format(k, seqs))

    TMP.close()

    # get script directory
    scripts = os.popen('which moR.py').read()
    # scripts = re.sub('moR.py', '', scripts)
    # scripts = re.sub('\s+', '', scripts)
    scripts = os.path.dirname(scripts) + '/'

    # bowtie index is placed in folder indexes in folder that holds the moRNA Finder
    # scripts
    tmp = None
    print_stderr("Build bowtie index of Rfam entries\n\n")
    if not os.path.isdir("{}indexes".format(scripts)):
        tmp = os.popen('mkdir "{}indexes"'.format(scripts)).read()
        if tmp:
            print_stderr(tmp, "\n")

    if not os.path.isfile("{}indexes/Rfam_index.1.ebwt".format(scripts)):
        err = os.system(
            'bowtie-build {} {}indexes/Rfam_index'.format(options.get('-r'), scripts))

        if not os.path.isfile("{}indexes/Rfam_index.1.ebwt".format(scripts)):
            print_stderr("{}\n\nRfam index could not be created in {}indexes/\nPlease check if moRNA Finder is allowed to create the directory {}indexes/ and write to it\nThe Rfam analysis will be skipped as long as this is not possible\n\n".format(
                err,
                scripts,
                scripts,
            ))
            return 256

    print_stderr("Mapping mature,star and loop sequences against index\n")
    IN = None

    # I think 0 MM would be too conservative
    if options.get('-q'):
        err = os.system('bowtie -f -v 1 -a --best --strata --norc {}indexes/Rfam_index expression_analyses/expression_analyses_{}/identified_precursors.fa expression_analyses/expression_analyses_{}/rfam_vs_precursor.bwt'.format(scripts, ltime, ltime))
        IN = open_or_die2(
            "expression_analyses/expression_analyses_{}/rfam_vs_precursor.bwt".format(ltime), 'rb')
    else:
        err = os.system('bowtie -f -v 1 -a --best --strata --norc {}indexes/Rfam_index moR_runs/run_{}/identified_precursors.fa moR_runs/run_{}/rfam_vs_precursor.bwt'.format(scripts, ltime, ltime))
        IN = open_or_die2(
            "moR_runs/run_{}/rfam_vs_precursor.bwt".format(ltime), 'rb')

    line = []
    ids = []
    while True:
        l = IN.readline()
        if not l:
            break

        line = re.split(r'\t', l)
        ids = re.split(r':', line[0])

        if re.search('rRNA', line[2], re.IGNORECASE):
            create_hash_key_chain(_hash, 0, ids[1], ids[0], 'rRNA', 'c')
            _hash[ids[1]][ids[0]]['rRNA']['c'] += 1
        elif re.search('tRNA', line[2], re.IGNORECASE):
            create_hash_key_chain(_hash, 0, ids[1], ids[0], 'tRNA', 'c')
            _hash[ids[1]][ids[0]]['tRNA']['c'] += 1
        else:
            _hash[ids[1]][ids[0]]['rRNA']['c'] = 0
            _hash[ids[1]][ids[0]]['tRNA']['c'] = 0
            _hash[ids[1]]['rfam'] = ''

    IN.close()

    _str = ('M', 'L', 'S')
    count_o = 0

    for k in _hash.keys():
        count_o = 0

        for i in range(0, 3):
            try:
                if _hash[k][_str[i]]['rRNA']['c'] > 0:
                    count_o += 1
            except KeyError:
                pass

        if count_o > 1:
            _hash[k]['rfam'] = 'rRNA'
            continue

        count_o = 0
        for i in range(0, 3):
            try:
                if _hash[k][_str[i]]['rRNA']['c'] > 0 or _hash[k][_str[i]]['tRNA']['c'] > 0:
                    count_o += 1
            except KeyError:
                pass

        if count_o > 1:
            _hash[k]['rfam'] = 'rRNA/tRNA'
            continue

        count_o = 0
        for i in range(0, 3):
            try:
                if _hash[k][_str[i]]['tRNA']['c'] > 0:
                    count_o += 1
            except KeyError:
                pass

        if count_o > 1:
            _hash[k]['rfam'] = 'tRNA'
            continue

    return 0


def getEvalue(v, dig):
    if v == 0:
        return 0

    sign = ''
    if v < 0:
        sign = '-'
        v = abs(v)

    count = 0

    while v >= 10:
        count += 1
        v /= 10

    while v < 0.1:
        count -= 1
        v *= 10
        if v >= 1:
            break

    v = str(v)
    if not re.search(r'\.', v):
        v += "."

    if len(v) > (dig + 2):
        v = substr(v, 0, (dig + 2))
    else:
        v += "0" * ((dig + 2) - len(v))

    if not count:
        return "{}{}".format(sign, v)
    else:
        if count > 0:
            return "{}{}e+{}\n".format(sign, v, count)
        else:
            return "{}{}e{}\n".format(sign, v, count)


def CloseHTML(HTML):
    HTML.write('''</table>
        </body>
    </html>''')


def PrintHtmlTableHeader(hl='', csv=''):
    global HTML
    h = {}

    # divide string by linebreaks every x characters
    p1 = '<th><a href="" class="tooltip3">'
    p11 = '<th><a href="" class="tooltip4">'
    p2 = '<span>'
    q = '</span></a></th>'

    if re.search('novel', hl, re.IGNORECASE):
        for i in range(1, 18):
            h[i] = {}

        h[1][1] = 'provisional id'
        h[1][2] = 'this is a provisional miRNA name assigned by moRNA Finder. The first part of the id designates the chromosome or genome contig on which the miRNA gene is located. The second part is a running number that is added to avoid identical ids. The running number is incremented by one for each potential miRNA precursor that is excised from the genome. Clicking this field will display a pdf of the structure, read signature and score breakdown of the reported miRNA.'
        h[2][1] = 'moRNA Finder score'
        h[2][2] = 'the log-odds score assigned to the hairpin by moRNA Finder'
        h[3][1] = 'estimated probability that the miRNA candidate is a true positive'
        h[3][2] = 'the estimated probability that a predicted novel miRNA with a score of this or higher is a true positive. To see exactly how this probability is estimated, mouse over the \'novel miRNAs, true positives\' in the table at the top of the webpage.'
        h[4][1] = 'rfam alert'
        h[4][
            2] = 'this field indicates if the predicted miRNA hairpin has sequence similarity to reference rRNAs or tRNAs. Warnings in this field should overrule the estimated probability that a reported miRNA is a true positive (previous field).'
        h[5][1] = 'total read count'
        h[5][2] = 'this is the sum of read counts for the predicted mature, loop and star miRNAs.'
        h[6][1] = 'mature read count'
        h[6][2] = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted mature miRNA, including 2 nts upstream and 5 nts downstream.'
        h[7][1] = 'loop read count'
        h[7][2] = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted miRNA loop, including 2 nts upstream and 5 nts downstream.'
        h[8][1] = 'star read count'
        h[8][2] = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted star miRNA, including 2 nts upstream and 5 nts downstream.'
        h[9][1] = 'significant randfold p-value'
        h[9][
            2] = 'this field indicates if the estimated randfold p-value of the excised potential miRNA hairpin is equal to or lower than 0.05 (see Bonnet et al., Bioinformatics, 2004).'
        h[10][1] = 'miRBase miRNA'
        h[10][
            2] = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA map to the miRNA hairpin, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs for the species is displayed.'
        h[11][1] = 'example miRBase miRNA with the same seed'
        h[11][2] = 'this field displays the ids of any reference mature miRNAs from related species that have a seed sequence identical to that of the reported mature miRNA. The seed is here defined as nucleotides 2-8 from the 5\' end of the mature miRNA. If more than one reference mature miRNA have identical seed, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs from related species is displayed.'
        h[12][1] = 'UCSC browser'
        h[12][2] = 'if a species name was input to moRNA Finder, then clicking this field will initiate a UCSC blat search of the consensus precursor sequence against the reference genome.'
        h[13][1] = 'NCBI blastn'
        h[13][
            2] = 'clicking this field will initiate a NCBI blastn search of the consensus precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).'
        h[14][1] = 'consensus mature sequence'
        h[14][2] = 'this is the consensus mature miRNA sequence as inferred from the deep sequencing reads.'
        h[15][1] = 'consensus star sequence'
        h[15][2] = 'this is the consensus star miRNA sequence as inferred from the deep sequencing reads.'
        h[16][1] = 'consensus precursor sequence'
        h[16][2] = 'this is the consensus precursor miRNA sequence as inferred from the deep sequencing reads. Note that this is the inferred Drosha hairpin product, and therefore does not include substantial flanking genomic sequence as does most miRBase precursors.'
        h[17][1] = 'precursor coordinate'
        h[17][2] = 'The given precursor coordinates refer do absolute position in the mapped reference sequence'

    elif re.search(r'miRBase miRNAs detected by moRNA Finder', hl, re.IGNORECASE):
        for i in range(1, 18):
            h[i] = {}

        h[1][1] = 'tag id'
        h[1][2] = 'this is a tag id assigned by moRNA Finder. The first part of the id designates the chromosome or genome contig on which the miRNA gene is located. The second part is a running number that is added to avoid identical ids. The running number is incremented by one for each potential miRNA precursor that is excised from the genome. Clicking this field will display a pdf of the structure, read signature and score breakdown of the miRNA.'
        h[2][1] = 'moRNA Finder score'
        h[2][2] = 'the log-odds score assigned to the hairpin by moRNA Finder'
        h[3][1] = 'estimated probability that the miRNA is a true positive'
        h[3][2] = 'the estimated probability that a predicted miRNA with a score of this or higher is a true positive. To see exactly how this probability is estimated, mouse over the \'novel miRNAs, true positives\' in the table at the top of the webpage. For miRBase miRNAs, this reflects the support that the data at hand lends to the miRNA.'
        h[4][1] = 'rfam alert'
        h[4][
            2] = 'this field indicates if the miRNA hairpin has sequence similarity to reference rRNAs or tRNAs. Warnings in this field should overrule the estimated probability that a reported miRNA is a true positive (previous field).'
        h[5][1] = 'total read count'
        h[5][2] = 'this is the sum of read counts for the mature, loop and star miRNAs.'
        h[6][1] = 'mature read count'
        h[6][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus mature miRNA, including 2 nts upstream and 5 nts downstream.'
        h[7][1] = 'loop read count'
        h[7][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus miRNA loop, including 2 nts upstream and 5 nts downstream.'
        h[8][1] = 'star read count'
        h[8][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus star miRNA, including 2 nts upstream and 5 nts downstream.'
        h[9][1] = 'significant randfold p-value'
        h[9][
            2] = 'this field indicates if the estimated randfold p-value of the miRNA hairpin is equal to or lower than 0.05 (see Bonnet et al., Bioinformatics, 2004).'
        h[10][1] = 'mature miRBase miRNA'
        h[10][
            2] = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA maps to the miRNA hairpin, then only the id of the reference miRBase miRNA that matches the predicted mature sequence is output.'
        h[11][1] = 'example miRBase miRNA with the same seed'
        h[11][2] = 'this field displays the ids of any reference mature miRNAs from related species that have a seed sequence identical to that of the reported mature miRNA. The seed is here defined as nucleotides 2-8 from the 5\' end of the mature miRNA. If more than one reference mature miRNA have identical seed, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs from related species is displayed.'
        h[12][1] = 'UCSC browser'
        h[12][2] = 'if a species name was input to moRNA Finder, then clicking this field will initiate a UCSC blat search of the consensus precursor sequence against the reference genome.'
        h[13][1] = 'NCBI blastn'
        h[13][
            2] = 'clicking this field will initiate a NCBI blastn search of the consensus precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).'
        h[14][1] = 'consensus mature sequence'
        h[14][2] = 'this is the consensus mature miRNA sequence as inferred from the deep sequencing reads.'
        h[15][1] = 'consensus star sequence'
        h[15][2] = 'this is the consensus star miRNA sequence as inferred from the deep sequencing reads.'
        h[16][1] = 'consensus precursor sequence'
        h[16][2] = 'this is the consensus precursor miRNA sequence as inferred from the deep sequencing reads. Note that this is the inferred Drosha hairpin product, and therefore does not include substantial flanking genomic sequence as does most miRBase precursors.'
        h[17][1] = 'precursor coordinate'
        h[17][2] = 'The given precursor coordinates refer do absolute position in the mapped reference sequence'
        h[4][3] = 'predicted mature seq. in accordance with miRBase mature seq.'
        h[4][4] = 'If the predicted moRNA Finder sequence overlaps with the miRBase annotated mature sequence than this is indicated by \'TRUE\'. If the predicted moRNA Finder star sequence overlaps with the miRBase annotated mature sequence this is inidicated by \'STAR\'.'
    else:
        for i in range(1, 18):
            h[i] = {}

        h[1][1] = 'miRBase precursor id'
        h[1][2] = 'Clicking this field will display a pdf of the structure and read signature of the miRNA.'
        h[2][1] = '-'
        h[2][1] += "&#160" * 5
        h[2][2] = '-'
        h[3][1] = '-'
        h[3][1] += "&#160" * 5
        h[3][2] = '-'
        h[4][1] = '-'
        h[4][1] += "&#160" * 5
        h[4][2] = '-'
        h[5][1] = 'total read count'
        h[5][2] = 'this is the sum of read counts for the mature and star miRNAs.'
        h[6][1] = 'mature read count(s)'
        h[6][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the mature miRNA, including 2 nts upstream and 5 nts downstream. If more than one mature sequence is given this will be a comman separated list'
        h[7][1] = '-    '
        h[7][2] = '-    '
        h[8][1] = 'star read count'
        h[8][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the star miRNA, including 2 nts upstream and 5 nts downstream. This field is empty unless a reference star miRNA was given as input to quantifier.py. If more than one mature sequence is given this will be a comman separated list'
        h[9][1] = 'remaining reads'
        h[9][2] = 'this is the number of reads that did not map to any of the mature and star sequences'
        h[10][1] = '-'   # 'miRBase mature id'
        h[10][2] = '-'   # 'Clicking this field will link to miRBase.'
        h[11][1] = '-'
        h[11][2] = '-'
        h[12][1] = 'UCSC browser'
        h[12][2] = 'if a species name was input to moRNA Finder, then clicking this field will initiate a UCSC blat search of the miRNA precursor sequence against the reference genome.'
        h[13][1] = 'NCBI blastn'
        h[13][
            2] = 'clicking this field will initiate a NCBI blastn search of the miRNA precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).'
        h[14][1] = 'miRBase mature sequence(s)'
        h[14][
            2] = 'this is/are the mature miRNA sequence(s) input to quantifier.py.'
        h[15][1] = 'miRBase star sequence(s)'
        h[15][
            2] = 'this is/are the star miRNA sequence(s) input to quantifier.py. This field is empty unless a reference star miRNA was given as input to quantifier.py.'
        h[16][1] = 'miRBase precursor sequence'
        h[16][2] = 'this is the precursor miRNA sequence input to quantifier.py.'
        h[17][1] = 'precursor coordinate'
        h[17][2] = 'The given precursor coordinates refer do absolute position in the mapped reference sequence'

    if csv:
        f = 1
        CSV.write('{}\n'.format(hl))
        for k in hash_sort_key(h, lambda x: x[0]):
            if f:
                f = 0
                CSV.write(h[k][1])
            else:
                if re.search('miRBase miRNAs detected', hl, re.IGNORECASE) and k == 5:
                    try:
                        CSV.write(h['mirbase'][1])
                    except KeyError:
                        CSV.write('')

                CSV.write('\t{}'.format(h[k][1]))
        CSV.write('\n')
        return

    HTML.write('''<br>
    <br>
    <br><h2>{}</h2><br>
    <font face="Times New Roman" size="2">
    <table border="1">'''.format(hl))

    for k in sorted(h.keys()):
        # if k == 2 or k == 3 or k == 4 or k == 7 or k == 10 or k == 11:
        #     continue

        if k != 16:
            if k == 9 and not re.search(r'not', hl):
                HTML.write('<th><a href=\"http://www.ncbi.nlm.nih.gov/entrez/utils/fref.fcgi?PrId=3051&itool=AbstractPlus-def&uid=15217813&nlmid=9808944&db=pubmed&url=http://bioinformatics.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=15217813\" target=\"blank\"class=\"tooltip3\">{}{}{}{}\n\n'.format(
                    h[k][1],
                    p2,
                    h[k][2],
                    q
                ))
            else:
                HTML.write('{}{}{}{}{}\n\n'.format(
                    p1,
                    h[k][1],
                    p2,
                    h[k][2],
                    q
                ))
                if re.search('miRBase miRNAs detected', hl) and k == 4:
                    HTML.write("{}{}{}{}{}\n\n".format(
                        p1,
                        h[k][3],
                        p2,
                        h[k][4],
                        q
                    ))
        else:
            HTML.write('{}{}{}{}{}\n\n'.format(
                p11,
                h[k][1],
                p2,
                h[k][2],
                q
            ))


def ReadInParameters():
    global options, HTML, ltime

    if os.path.isfile("moR_runs/run_{}/run_{}_parameters".format(ltime, ltime)):
        HTML.write("<h2>Parameters used</h2>\n")
        HTML.write(" <table border=\"0\">\n")
        HTML.write(
            "<tr><td>moRNA Finder version</td><td>{}</td></tr><br>\n".format(options.get('-V')))

        PAR = open_or_die2(
            'moR_runs/run_{}/run_{}_parameters'.format(ltime, ltime), 'rb')

        while True:
            l = PAR.readline()
            if not l:
                break

            m = re.search(r'args\s+(\S.+)', l)
            if m:
                HTML.write(
                    "<tr><td>Program call</td><td>{}</td><tr>\n".format(m.groups()[0]))
                continue

            m = re.match(r'^dir with tmp files\s+(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Working directory</td><td>{}<td></tr><br>\n".format(m.groups()[0]))
                continue

            m = re.match(r'^file_reads\s+(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Reads</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.match(r'^file_genome\s+(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Genome</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.match(r'^file_reads_vs_genome\s+(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Mappings</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.match(r'^file_mature_ref_this_species\s+(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Reference mature miRNAs</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.match(r'^file_mature_ref_other_species\s+(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Other mature miRNAs</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'option -t\s*=\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Species</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'option -q\s*=\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Expression analysis raw file</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'option -a\s*=\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>minimum read stack height</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'option -b\s*=\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>minimum score for novel precursors shown in table</td><td>{}</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'option -c\s*=\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>randfold analysis disabled</td><td>yes</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'option -v\s*=\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>remove temporary files</td><td>yes</td></tr>\n".format(m.groups()[0]))
                continue

            m = re.search(r'started:\s*(\S+)', l)
            if m:
                HTML.write(
                    "<tr><td>Start</td><td>{}</td></tr>\n".format(m.groups()[0]))
                tmp = PAR.readline()
                m = re.search(r'ended:\s*(\S+)', tmp)
                if m:
                    HTML.write(
                        "<tr><td>End</td><td>{}</td></tr>\n".format(m.groups()[0]))

                tmp = PAR.readline()
                m = re.search(r'total:\s*(\S+)', tmp)
                if m:
                    HTML.write(
                        "<tr><td>Total</td><td>{}</td></tr>\n".format(m.groups()[0]))

                continue

        PAR.close()
        HTML.write("</table><br><br>")
    else:
        print_stderr(
            "File moR_runs/run_{}/run_{}_parameters not found\n".format(ltime, ltime))


def ClosePDF(drawCanvas):
    '''
    Save the canvas to file
    '''
    drawCanvas.save()


def Line(drawCanvas, x1, y1, x2, y2, color, width):
    '''
    Draw a line
    '''
    drawCanvas.setStrokeColor(color)
    drawCanvas.setLineWidth(width)
    drawCanvas.line(x1, y1, x2, y2)


def Base(drawCanvas, x1, y1, base, color, size):
    '''
    Draw a letter
    '''
    drawCanvas.setFont('Courier-Bold', size)
    drawCanvas.setFillColor(color)
    drawCanvas.drawString(x1, y1, base)


def Shifting():
    global xc, yc, minx, miny, maxx, maxy
    minx = min(xc.values())
    maxx = max(xc.values())
    miny = min(yc.values())
    maxy = max(yc.values())
    shiftx = abs(minx) + 10
    shifty = abs(miny) + 10

    if minx < 0:
        for i in range(0, len(rna_d) - 1):
            xc[i] += shiftx

    if miny < 0:
        for i in range(0, len(rna_d) - 1):
            yc[i] += shifty

    if maxx > 600:
        for i in range(0, len(rna_d) - 1):
            xc[i] -= minx + 10

    if maxy > 800:
        for i in range(0, len(rna_d) - 1):
            yc[i] -= miny + 10


def drawLabel(drawCanvas, x, y, fontName, fontSize, string, color):
    '''
    Draw a label
    '''
    c = drawCanvas
    c.setFont(fontName, fontSize)
    c.setFillColor(color)
    c.drawString(x, y, str(string))


def CreatePDF(_hash, filename):
    global trb, rna_d, xc, yc, bpo2r, bpo1r, bpo1, bpo2, y, sid, hairpin2mature
    fname = '{}/pdfs_{}/{}.pdf'.format(cwd, ltime, filename)
    c = canvas.Canvas(fname, pagesize=A4)
    spacer = len(sid)
    # c.setFont('Times-Roman', 20)
    trb = 'Times-Roman'
    madd = 0

    if mirbase:
        madd = 60

    drawLabel(c, xposshift + 20, y + 300 + downy,
              trb, 8, "Provisional ID", black)
    drawLabel(c, xposshift + 100, y + 300 + downy, trb, 8, ": " + sid, black)

    # only print for discovered miRNAs
    if not mirbase:
        spaces = " " * (spacer - len(_hash[sid]['score']))
        drawLabel(c, xposshift + 20, y + 290 + downy,
                  trb, 8, "Score total", black)
        drawLabel(c, xposshift + 100, y + 290 + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['score']), black)

        spaces = " " * (spacer - len(_hash[sid]['score_star']))
        drawLabel(c, xposshift + 20, y + 280 + downy,
                  trb, 8, "Score for star read(s)", black)
        drawLabel(c, xposshift + 100, y + 280 + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['score_star']), black)

        spaces = " " * (spacer - len(_hash[sid]['score_read']))
        drawLabel(c, xposshift + 20, y + 270 + downy,
                  trb, 8, "Score for read counts", black)
        drawLabel(c, xposshift + 100, y + 270 + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['score_read']), black)

        spaces = " " * (spacer - len(_hash[sid]['score_mfe']))
        drawLabel(c, xposshift + 20, y + 260 + downy,
                  trb, 8, "Score for mfe", black)
        drawLabel(c, xposshift + 100, y + 260 + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['score_mfe']), black)

        spaces = " " * (spacer - len(_hash[sid]['rand']))
        drawLabel(c, xposshift + 20, y + 250 + downy,
                  trb, 8, "Score for randfold", black)
        drawLabel(c, xposshift + 100, y + 250 + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['rand']), black)

        spaces = " " * (spacer - len(_hash[sid].get('score_cons', '')))
        drawLabel(c, xposshift + 20, y + 240 + downy,
                  trb, 8, "Score for cons. seed", black)
        drawLabel(c, xposshift + 100, y + 240 + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid].get('score_cons', '')), black)

    spaces = " " * (spacer - len(_hash[sid]["freq_total"]))
    drawLabel(c, xposshift + 20, y + 230 + madd +
              downy, trb, 8, "Total read count", black)
    drawLabel(c, xposshift + 100, y + 230 + madd + downy, trb, 8,
              ": {}{}".format(spaces, _hash[sid]['freq_total']), black)

    spaces = " " * (spacer - len(_hash[sid]["freq_mature"]))
    drawLabel(c, xposshift + 20, y + 220 + madd +
              downy, trb, 8, "Mature read count", black)
    drawLabel(c, xposshift + 100, y + 220 + madd + downy, trb, 8,
              ": {}{}".format(spaces, _hash[sid]['freq_mature']), black)

    if not mirbase:
        spaces = " " * (spacer - len(_hash[sid]["freq_loop"]))
        drawLabel(c, xposshift + 20, y + 210 + madd +
                  downy, trb, 8, "Loop read count", black)
        drawLabel(c, xposshift + 100, y + 210 + madd + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['freq_loop']), black)

        spaces = " " * (spacer - len(_hash[sid]["freq_star"]))
        drawLabel(c, xposshift + 20, y + 200 + madd +
                  downy, trb, 8, "Star read count", black)
        drawLabel(c, xposshift + 100, y + 200 + madd + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['freq_star']), black)
    else:
        spaces = " " * (spacer - len(_hash[sid]["freq_star"]))
        drawLabel(c, xposshift + 20, y + 210 + madd +
                  downy, trb, 8, "Star read count", black)
        drawLabel(c, xposshift + 100, y + 210 + madd + downy, trb, 8,
                  ": {}{}".format(spaces, _hash[sid]['freq_star']), black)

    trb = 'Courier-Bold'
    return c


def CreateHistogram(drawCanvas, _hash):
    global y, order, hstruct, xposshift, lstruct_multi, totalreads
    trb = 'Courier-Bold'
    c = drawCanvas
    # c.setLineWidth(2)
    # c.setStrokeColor(black)
    y = yorig - downy
    Line(c, xposshift + 20, y + 160, xposshift + 20, y + 50 - 1, black, 2)
    Line(c, xposshift + 20, y + 50, xposshift +
         20 + lstruct_multi, y + 50, grey, 2)

    c.setLineWidth(2)
    c.setStrokeColor(black)
    p = c.beginPath()
    p.moveTo(xposshift + 17 + lstruct_multi, y + 53)
    p.lineTo(xposshift + 20 + lstruct_multi, y + 50)
    p.lineTo(xposshift + 17 + lstruct_multi, y + 47)
    c.drawPath(p)

    p = c.beginPath()
    p.moveTo(xposshift + 17, y + 157)
    p.lineTo(xposshift + 20, y + 160)
    p.lineTo(xposshift + 23, y + 157)
    c.drawPath(p)

    p = c.beginPath()
    p.moveTo(xposshift + 17, y + 150)
    p.lineTo(xposshift + 23, y + 150)
    c.drawPath(p)

    drawLabel(c, xposshift + 12, y + 165, trb, 6, 'freq.', black)
    drawLabel(c, xposshift + lstruct_multi, y + 40, trb, 8, 'length', black)
    drawLabel(c, xposshift + 10, y + 148, trb, 6, '1', black)

    # 0.75
    p = c.beginPath()
    p.moveTo(xposshift + 17, y + 125)
    p.lineTo(xposshift + 23, y + 125)
    c.drawPath(p)
    drawLabel(c, xposshift + 2, y + 122, trb, 6, '0.75', black)

    # 0.5
    p = c.beginPath()
    p.moveTo(xposshift + 17, y + 100)
    p.lineTo(xposshift + 23, y + 100)
    c.drawPath(p)
    drawLabel(c, xposshift + 6, y + 98, trb, 6, '0.5', black)

    # 0.25
    p = c.beginPath()
    p.moveTo(xposshift + 17, y + 75)
    p.lineTo(xposshift + 23, y + 75)
    c.drawPath(p)
    drawLabel(c, xposshift + 2, y + 73, trb, 6, '0.25', black)

    # 0
    drawLabel(c, xposshift + 12, y + 48, trb, 6, '0', black)

    # draw flank1
    p = c.beginPath()
    p.moveTo(xposshift + 20, y + 50)

    # example case for mmu-mir-497
    # 0..12 means to 13th char in string ## print one char further to have the transition
    # for i in range(0, lflank1 + 2):
    i = 0
    while i <= lflank1:
        p.lineTo(position_hash[
                 i + 1], ((float(hstruct[i]) / totalreads) * 100) + y + 50)    # .25
        i += 1

    c.drawPath(p)
    i -= 1  # $i is 13 now where mature starts

    lastx = position_hash[i + 1]
    lasty = ((float(hstruct[i]) / totalreads) * 100) + y + 50

    c.setLineWidth(2)
    _sorted = hash_sort_key(order, lambda x: x[1])

    # import pdb
    # pdb.set_trace()

    for k in _sorted:
        if k == "m":
            c.setStrokeColor(col_mature)
            p = c.beginPath()
            p.moveTo(lastx, lasty)

            # graphical output corrections
            if mb > sb:
                # for i in range(mb + 1, mb + lmature):
                i = mb + 1
                while i < mb + lmature:
                    p.lineTo(position_hash[
                             i + 1], ((float(hstruct[i]) / totalreads) * 100) + y + 50)   # .25
                    i += 1
            else:
                # for i in range(mb, mb + lmature):
                i = mb
                while i < mb + lmature:
                    p.lineTo(position_hash[
                             i + 1], ((float(hstruct[i]) / totalreads) * 100) + y + 50)   # .25
                    i += 1

            i -= 1    # $i is 34 now
            c.drawPath(p)

            lastx = position_hash[i + 1]
            lasty = ((float(hstruct[i]) / totalreads) * 100) + y + 50
            c.setStrokeColor(black)

        elif k == "s":

            defined = False
            try:
                _hash[sid]['obs']
                defined = True
            except KeyError:
                pass

            if defined and _hash[sid]['obs']:
                c.setStrokeColor(col_star_obs)
            else:
                c.setStrokeColor(col_star_exp)

            p = c.beginPath()
            p.moveTo(lastx, lasty)

            # graphical output corrections
            if sb > mb:
                # for i in range(sb + 1, sb + lstar):
                i = sb + 1
                while i < sb + lstar:
                    p.lineTo(position_hash[
                             i + 1], ((float(hstruct[i]) / totalreads) * 100) + y + 50)  # .25
                    i += 1
            else:
                # for i in range(sb, sb + lstar):
                i = sb
                while i < sb + lstar:
                    p.lineTo(position_hash[
                             i + 1], ((float(hstruct[i]) / totalreads) * 100) + y + 50)  # .25
                    i += 1

            i -= 1
            c.drawPath(p)
            lastx = position_hash[i + 1]
            lasty = ((float(hstruct[i]) / totalreads) * 100) + y + 50
            c.setStrokeColor(black)

        elif k == "l":
            c.setStrokeColor(orange)
            p = c.beginPath()
            p.moveTo(lastx, lasty)

            # for i in range(lb + 1, lb + lloop + 1):
            i = lb + 1
            while i <= lb + lloop:
                p.lineTo(position_hash[
                         i + 1], ((float(hstruct[i]) / totalreads) * 100) + y + 50)   # .25
                i += 1

            i -= 1
            c.drawPath(p)

            lastx = position_hash[i + 1]
            lasty = ((float(hstruct[i]) / totalreads) * 100) + y + 50
            c.setStrokeColor(black)
        elif k == "f":
            c.setStrokeColor(black)
            p = c.beginPath()
            p.moveTo(lastx, lasty)
            # for i in range(fl2b + 1, fl2b + lflank2 + 1):
            i = fl2b + 1
            while i <= fl2b + lflank2:
                p.lineTo(position_hash[
                         i + 1], ((float(hstruct.get(i, 0)) / totalreads) * 100) + y + 50)   # .25
                i += 1

            i -= 1
            c.drawPath(p)
            lastx = position_hash[i + 1]
            lasty = ((float(hstruct.get(i, 0)) / totalreads) * 100) + y + 50
            c.setStrokeColor(black)


def DrawStructure(drawCanvas, filename):
    global minx, miny, maxx, maxy, star_exp_hit_pos, rna_d, xc, yc, bpo2r, bpo1r, bpo1, bpo2, y, sid, hairpin2mature
    c = drawCanvas
    os.system("RNAfold -d 0 < {}/pdfs_{}/{}.tmp > {}/pdfs_{}/tmp".format(cwd,
                                                                         ltime, filename, cwd, ltime))

    in_pos = 0
    in_pairs = 0
    count = 0

    bp = {}
    PS = open_or_die2('{}/pdfs_{}/rna.ps'.format(cwd, ltime), 'rb')
    minx, miny = 10000, 10000
    maxx, maxy = 0, 0

    centering_x = 0
    twisted = 0
    sums = 0
    while True:
        l = PS.readline()
        if not l:
            break

        if re.search(r'\/sequence', l):
            line = PS.readline().strip()
            rna_d = ssplit(line)
            continue

        if re.search(r'\/coor\s*\[', l):
            in_pos = 1
            continue

        m = re.search(r'\[(\S+)\s+(\S+)\]', l)
        if in_pos and m:
            xc[count] = float(m.groups()[0])  # x coordinate
            yc[count] = float(m.groups()[1])  # y coordinate
            count += 1
            continue

        if re.search(r'\/pairs', l):   # read in base pairs
            count = 0
            in_pos = 0
            in_pairs = 1
            continue

        # mature begin is after loop begin
        if mb > lb:
            twisted = 1

        m = re.search(r'\[(\S+)\s+(\S+)\]', l)
        if in_pairs and m:
            if twisted:
                bpo2r = int(m.groups()[0])
                bpo1r = int(m.groups()[1])
            else:
                bpo2r = int(m.groups()[1])
                bpo1r = int(m.groups()[0])

            # determine two subsequence in mature having subsequent paired
            # bases
            if bpo1r >= mb - offset and bpo1r < me - offset and centering_x == 0:
                if twisted:
                    sums = -1
                else:
                    sums = 1

                if (bpo1r - bpo1) == sums and (bpo2 - bpo2r) == sums:
                    if twisted:
                        centering_x = bpo1r
                    else:
                        centering_x = bpo1r - 2

            bpo1 = bpo1r
            bpo2 = bpo2r
            bp[bpo1r - 1] = bpo2r - 1  # saving nt pairts in hash bp
            continue

        if in_pairs and re.search(r'\]\s*def', l):
            in_pairs = 0
            break

    PS.close()

    # print(xc)
    # print(yc)

    Shifting()

    minx = min(xc.values())
    maxx = max(xc.values())
    miny = min(yc.values())
    maxy = max(yc.values())

    # determine if mirror or not
    # mirror sequence so that loop is on the right hand side
    mir = 0

    if twisted:
        mir = 1

    yshift = 0
    cshift = 3

    if mir:
        for i in range(0, len(rna_d) - 1):
            xc[i] *= -1

    minx = min(xc.values())
    maxx = max(xc.values())
    miny = min(yc.values())
    maxy = max(yc.values())

    if mir:
        for i in range(0, len(rna_d) - 1):
            xc[i] += abs(minx) + 10

    minx = min(xc.values())
    maxx = max(xc.values())
    miny = min(yc.values())
    maxy = max(yc.values())

    ax = xc[centering_x]
    ay = yc[centering_x]

    # point relative to center
    bx = 0
    by = 0

    if twisted:
        bx = xc.get(centering_x - 1, 0)
        by = yc.get(centering_x - 1, 0)
    else:
        bx = xc.get(centering_x + 1, 0)
        by = yc.get(centering_x + 1, 0)

    gk = by - ay
    ak = bx - ax

    r = math.sqrt((ak**2 + gk**2))
    phi = math.asin(gk / r)

    if bx < ax and by > ay:
        phi = math.pi - phi

    if bx <= ax and by <= ay:
        phi *= -1
        phi += math.pi

    # print('cen:', centering_x, 't:', twisted, 'bx:', bx, 'by:', by, 'gk:', gk, 'ak:', ak, 'phi:', phi, 'ay:', ay, 'ax:', ax)

    alpha = None
    do_rot = 1
    if do_rot:
        last = xc[0]

        # rotate every point in a designated angle of phi
        for i in range(0, len(rna_d) - 1):
            if i == centering_x:
                continue

            bx = xc[i]
            by = yc[i]

            gk = by - ay
            ak = bx - ax
            r = math.sqrt(ak**2 + gk**2)
            alpha = math.asin(gk / r)

            if bx < ax and by > ay:
                alpha = math.pi - alpha

            if bx <= ax and by <= ay:
                alpha *= -1
                alpha += math.pi

            alpha -= phi
            xc[i] = ax + r * math.cos(alpha)
            yc[i] = ay + r * math.sin(alpha)

            dif = xc[i] - last
            last = xc[i]

    _reduce = 0
    # import pdb
    # pdb.set_trace()
    # red_dist = abs(xc[mb + cshift] - xc[mb - 1 + cshift])

    Shifting()

    if not twisted:
        bpkeys = sorted(bp.keys())
        maxy = max(yc.values())
        if yc[bpkeys[0]] < yc[bp[bpkeys[0]]]:
            for i in range(0, len(rna_d) - 1):
                yc[i] = - yc[i] + maxy + 10

    minx = min(xc.values())
    maxx = max(xc.values())
    miny = min(yc.values())
    maxy = max(yc.values())

    y = yorig + 300
    x = 550
    rx = 300

    if maxx < x:
        x = x - maxx
    else:
        x = -abs(x - maxx)

    if maxy < y:
        y = y - maxy
    else:
        y = -abs(y - maxy)

    scx, scy = 250, 250
    tx, ty = 0, 0
    scfactor = rx / float(maxx - minx)

    if scfactor < 1:
        for i in range(0, len(rna_d) - 1):
            tx = xc[i] - scx
            ty = yc[i] - scy
            xc[i] = tx * scfactor + scx
            yc[i] = ty * scfactor + scx

    minx = min(xc.values())
    maxx = max(xc.values())
    miny = min(yc.values())
    maxy = max(yc.values())

    y = yorig + 280
    x = 550

    if maxx < x:
        x = x - maxx
    else:
        x = - abs(x - maxx)

    if maxy < y:
        y = y - maxy
    else:
        y = - abs(y - maxy)

    if minx + x < x - rx:
        for k in xc.keys():
            xc[k] += (x - rx) - (minx + x)

    for i in range(1, len(rna_d) - 1):
        # print(xc[i - 1] + x + 2.5, yc[i - 1] + 2 + y)
        Line(c, xc[i - 1] + x + 2.5, yc[i - 1] + 2 + y,
             xc[i] + 2.5 + x, yc[i] + 2 + y, grey, 0.5)

    Base(c, xc[0] + x, yc[0] + y, "{}'".format(rna_d[0]), black, 8)
    for i in range(1, len(rna_d) - 2):
        # ATTR : star_exp_hit_pos is always empty
        # if star_exp_hit_pos[i - 1 + offset]:
        #     Base(c, xc[i] + x, yc[i] + y, rna_d[i], col_star_exp, 8)
        # else:
        defined = False
        try:
            star_exp_hit_pos[i - 1 + offset]
            defined = True
        except KeyError:
            pass

        if defined and star_exp_hit_pos[i - 1 + offset]:
            Base(c, xc[i] + x, yc[i] + y, rna_d[i], col_star_exp, 8)
        else:
            Base(c, xc[i] + x, yc[i] + y, rna_d[i],
                 assign_str[i - 1 + offset], 8)

    Base(c, xc[len(rna_d) - 2] + x, yc[len(rna_d) - 2] +
         y, "{}'".format(rna_d[len(rna_d) - 2]), black, 8)

    # drawing the bp lines
    scfactorl = 0.4
    fx = 0                     # from x coordinate
    tox = 0                    # from y coordinate
    fy = 0                     # from y cooridinate
    toy = 0                    # to y coordinate
    dx = 0                     # xlength
    dy = 0                     # y length
    dx1 = 0   # difference between orig x length and scaled x length
    dy1 = 0   # difference between orig y length and scaled y length

    for k in bp.keys():
        dx = abs(xc[k] - xc[bp[k]])
        dy = abs(yc[k] - yc[bp[k]])
        dx1 = (dx - scfactorl * dx) / float(2)
        dy1 = (dy - scfactorl * dy) / float(2)

        if xc[k] > xc[bp[k]]:
            fx = xc[k] - dx1
            tox = xc[bp[k]] + dx1
        else:
            fx = xc[k] + dx1
            tox = xc[bp[k]] - dx1

        if yc[k] > yc[bp[k]]:
            fy = yc[k] - dy1
            toy = yc[bp[k]] + dy1
        else:
            fy = yc[k] + dy1
            toy = yc[bp[k]] - dy1

        Line(c, fx + 2.5 + x, fy + 2 + y, tox +
             2.5 + x, toy + 2 + y, black, 0.5)


def CreateAlignment(drawCanvas, _hash):
    global y
    c = drawCanvas
    # draw left flank of sequence
    drawLabel(c, xposshift + (18 + ((mb + 1) * multiplier)),
              y + 40, trb, 6, mb - lflank1 + 1, black)
    drawLabel(c, xposshift + (20 + ((mb + 1) * multiplier)),
              y + 10, trb, 12, 'Mature', col_mature)
    drawLabel(c, xposshift + (18 + ((lb + 1) * multiplier)),
              y + 40, trb, 6, lb - lflank1 + 1, black)
    drawLabel(c, xposshift + (18 + ((sb + 1) * multiplier)),
              y + 40, trb, 6, sb - lflank1 + 1, black)

    defined = False
    try:
        _hash[sid]['obs']
        defined = True
    except KeyError:
        pass

    if defined and re.search(r'S', _hash[sid]['obs']):
        drawLabel(c, xposshift + (20 + ((sb_obs + 1) * multiplier)),
                  y + 10, trb, 12, 'Star', col_star_obs)
    elif re.search(r'S', _hash[sid]['exp']):
        drawLabel(c, xposshift + (20 + ((sb + 1) * multiplier)),
                  y + 10, trb, 12, 'Star', col_star_exp)

    drawLabel(c, xposshift + (18 + ((fl2b + 1) * multiplier)),
              y + 40, trb, 6, fl2b - lflank1 + 1, black)

    for i in range(0, len(rna)):
        drawLabel(c, position_hash[i], y, trb, 6, rna[i], assign_str[i])

    drawLabel(c, xposshift + 25 + lstruct_multi, y, trb, 6, '-3\'', black)
    drawLabel(c, xposshift + 10, y, trb, 6, '5\'-', black)

    defined = False
    try:
        _hash[sid]['obs']
        defined = True
    except KeyError:
        pass

    if defined and _hash[sid]['obs']:
        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'obs', black)
    else:
        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'exp', black)

    defined = False
    try:
        _hash[sid]['obs']
        defined = True
    except KeyError:
        pass

    if defined and _hash[sid]['obs']:
        y -= 10
        for i in range(0, len(rna)):
            drawLabel(c, position_hash[i], y, trb,
                      6, rna[i], assign_str_exp[i])

        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'exp', black)

    defined = False
    try:
        _hash[sid]['known']
        defined = True
    except KeyError:
        pass

    if defined and _hash[sid]['known']:
        y -= 10
        mid = _hash[sid]["known"]
        mse = knownones[mid]

        # now check if mature of moRNA Finder is same as mirbase or not
        mbaseb = _hash[sid]['pri_seq'].find(mse)
        if mbaseb >= 0:
            # draw flank first
            drawLabel(c, xposshift + 20 + multiplier, y, trb, 6,
                      substr(_hash[sid]['pri_seq'], 0, mbaseb), black)
            # draw mature part
            drawLabel(c, xposshift + 20 + ((mbaseb + 1) *
                                           multiplier), y, trb, 6, mse, green)
            # draw right flank last
            drawLabel(c, xposshift + 20 + ((len(mse) + mbaseb + 1) * multiplier),
                      y, trb, 6, substr(_hash[sid]['pri_seq'], mbaseb + len(mse)), black)
        else:
            print_stderr("mature sequence not found \n")

        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'known', green)

    y -= 10
    structx = ssplit(sstruct)
    sadd = 0

    drawLabel(c, position_hash[0], y, trb, 6, sstruct, black)

    drawLabel(c, xposshift + 30 + lstruct_multi, y, trb, 6, 'reads', black)
    drawLabel(c, xposshift + 70 + lstruct_multi, y, trb, 6, 'mm', black)
    drawLabel(c, xposshift + 110 + lstruct_multi, y, trb, 6, 'sample', black)
    y -= 10
    if options.get('-o') == '':
        for tag in _hash[sid]['reads'].keys():
            for k in hash_sort_key(hash2order, lambda x: x[1]):
                if hash2sample[k] != tag:
                    continue

                drawLabel(c, position_hash[0], y, trb, 6, hash2seq[k], black)
                # matches and read numbers
                drawLabel(c, xposshift + 30 + lstruct_multi, y, trb,
                          6, str(hash2c[tag][hash2key[k]]), black)
                drawLabel(c, xposshift + 70 + lstruct_multi,
                          y, trb, 6, str(hash2mm[k]), black)
                drawLabel(c, xposshift + 110 + lstruct_multi,
                          y, trb, 6, str(hash2sample[k]), black)

                y -= 10
                if y < 100:
                    # add a new page
                    c.showPage()
                    y = 800
                    for i in range(0, len(rna)):
                        drawLabel(c, position_hash[i], y, trb, 6, rna[
                                  i], assign_str[i])

                    drawLabel(c, xposshift + (20 + ((mb + 1) * multiplier)),
                              y + 10, trb, 12, 'Mature', col_mature)

                    defined = False
                    try:
                        _hash[sid]['obs']
                        defined = True
                    except KeyError:
                        pass

                    if defined and re.search(r'S', _hash[sid]['obs']):
                        drawLabel(c, xposshift + (20 + ((sb_obs + 1) * multiplier)),
                                  y + 10, trb, 12, 'Star', col_star_obs)
                    elif re.search(r'S', _hash[sid]['exp']):
                        drawLabel(c, xposshift + (20 + ((sb + 1) * multiplier)),
                                  y + 10, trb, 12, 'Star', col_star_exp)

                    y -= 10

            y -= 10

    else:
        for k in hash_sort_key(hash2order, lambda x: x[1]):
            tag = hash2sample[k]
            drawLabel(c, position_hash[0], y, trb, 6, hash2seq[k], black)

            # matches and read numbers
            drawLabel(c, xposshift + 30 + lstruct_multi, y,
                      trb, 6, hash2c[tag][hash2key[k]], black)
            drawLabel(c, xposshift + 70 + lstruct_multi,
                      y, trb, 6, hash2mm[k], black)
            drawLabel(c, xposshift + 110 + lstruct_multi,
                      y, trb, 6, hash2sample[k], black)

            y -= 10
            if y < 100:
                # Add a new page
                c.showPage()
                y = 800
                for i in range(0, len(rna)):
                    drawLabel(c, position_hash[i], y,
                              trb, 6, rna[i], assign_str[i])

                drawLabel(c, xposshift + (20 + ((mb + 1) * multiplier)),
                          y + 10, trb, 12, 'Mature', col_mature)

                if re.search(r'S', _hash[sid]['obs']):
                    drawLabel(c, xposshift + (20 + ((sb_obs + 1) * multiplier)),
                              y + 10, trb, 12, 'Star', col_star_obs)
                elif re.search(r'S', _hash[sid]['exp']):
                    drawLabel(c, xposshift + (20 + ((sb + 1) * multiplier)),
                              y + 10, trb, 12, 'Star', col_star_exp)

                y -= 10


def CreateStructurePDF(_hash):
    global sb, fl2b, offset, star_exp_hit_pos, sb_obs, lmature, lloop, order, lflank1, lflank2, hash2order, rna, hstruct, sstruct, y, yorig, sid, assign_str, assign_str_exp, lstruct_multi, totalreads, mb, lb, me, hash2c, hash2key, hash2order, hash2mm, hash2seq, hash2sample, lstar

    filename = None
    print_stderr('Creating PDF files\n')
    for k in hash_sort_key(_hash, lambda x: _hash[x[0]]['freq_total'] * -1):

        defined = False
        try:
            _hash[k]['pdf']
            defined = True
        except KeyError:
            pass

        if not (defined and _hash[k]['pdf']):
            continue
        if not _hash[k]['freq_total']:
            continue

        # print_stderr(k, ' pdf: ', _hash[k]['pdf'], " freq: ", _hash[k]['freq_total'], '\n')

        sid = k
        sid = tr(sid, '|', '_')
        star_exp_hit_pos = {}
        filename = sid

        # if seen[filename]:
        #     continue

        if os.path.isfile("{}/pdfs_{}/{}.pdf".format(cwd, ltime, filename)):
            continue

        # # reinit variables
        i = 0
        offset = 0
        me = 0       # mature end coordinate
        desc = []
        lflank1 = 0  # length(string of left flank)
        fl1 = 0      #
        lflank2 = 0  # length string of right flank
        fl2b = -1   # right flank begin
        lloop = 0   # string of loop
        lb = -1     # starting position of loop
        lstar = 0   # string of star sequence
        sb = -1     # starting
        lmature = 0  # string of mature
        mb = -1     # mature begin
        sstruct = ''     # structure string
        pri_seq = ""     # pri-cursor sequence
        lenstr = 0
        pdf = None     # pdf descriptor
        page = None    # page descriptor
        gfx = None     # graphic variable
        trb = None     # fontvariable
        hash2 = {}
        hash2c = {}
        hash2mm = {}
        hash2order = {}
        order = {}

        yorig = 500
        downy = 50

        dline = None                     # line graphic handler

        first = 1
        lastx = 0
        lasty = 0

        final = ""                       # final output string of a read
        pseq = []                        # precursor sequence
        rseq = []                        # read sequence

        totalreads = 0

        assign_str = {}
        assign_str_exp = {}
        hstruct = {}
        bpo1 = -10                        # left nt pos in first bp
        bpo2 = -10                        # right nt pos in first bp
        bpo1r = -10                       # left nt pos in second bp
        bpo2r = -10                       # right nt pos in second bp
        ffe = 0                           # first flank end position
        ff2b = 0                          # second flank begin position
        _sorted = []                      # array that stores sorted order of fl1,m,l,s,fl2
        y = yorig                         # y coordinate
        # min and max x,y coordinates of rna sequence
        (minx, miny, maxx, maxy) = None, None, None, None
        rna = []                         # rna sequence

        xc = {}
        yc = {}

        # main program;
        pri_seq = _hash[sid]["pri_seq"]
        rna = ssplit(pri_seq)
        pri_seq = pri_seq.strip()

        defined = False
        try:
            _hash[sid]['obs']
            defined = True
        except KeyError:
            pass

        if defined and _hash[sid]['obs']:
            desc2 = ssplit(_hash[sid]['exp'])
            for i in range(0, len(desc2)):
                if (desc2[i] == "f"):   # # assign_str now starts at 0 not at one
                    assign_str_exp[i] = black
                elif (desc2[i] == "M"):
                    assign_str_exp[i] = col_mature
                elif (desc2[i] == "l"):
                    assign_str_exp[i] = col_loop
                elif (desc2[i] == "S"):
                    assign_str_exp[i] = col_star_exp
                else:
                    print_stderr(
                        '"something went wrong while parsing alignment in output.mrd file\n"')

        # if observed star strand then use this one otherwise use the expected
        # one
        desc_tmp = []

        defined = False
        try:
            _hash[sid]['obs']
            defined = True
        except KeyError:
            pass

        if defined and _hash[sid]['obs']:
            desc = ssplit(_hash[sid]['obs'])
            desc_tmp = ssplit(_hash[sid]['exp'])
        else:
            desc = ssplit(_hash[sid]['exp'])

        run_var_loop = 0
        for i in range(0, len(desc)):
            if desc[i] == "f":   # assign_str now starts at 0 not at one
                assign_str[i] = "black"

                if mb == -1:           # if in first flank
                    if desc[i + 1] == "f":
                        ffe = i
                    lflank1 += 1
                else:
                    # second flank   ## dependent on exp/obs both or not
                    if star_exp_hit and i < len(desc_tmp) and desc_tmp[i] == "f":
                        if fl2b == -1:
                            fl2b = i

                        try:
                            if not order["f"]:
                                order["f"] = fl2b
                        except KeyError:
                            order["f"] = fl2b

                        lflank2 += 1
                    elif star_exp_hit and i < len(desc_tmp) and desc_tmp[i] == "S":
                        star_exp_hit_pos[i] = i      # star hit expression
                        continue
                    else:
                        if fl2b == -1:
                            fl2b = i

                        try:
                            if not order["f"]:
                                order["f"] = fl2b
                        except KeyError:
                            order["f"] = fl2b

                        lflank2 += 1
            elif desc[i] == 'M':
                if mb == -1:
                    mb = i
                if desc[i + 1] != "M":
                    me = i

                try:
                    if not order["m"]:
                        order["m"] = i
                except KeyError:
                    order["m"] = i

                assign_str[i] = col_mature
                lmature += 1
            elif desc[i] == 'l':
                if lb == -1:
                    lb = i

                try:
                    if not order["l"]:
                        order["l"] = i
                except KeyError:
                    order["l"] = i

                assign_str[i] = col_loop
                lloop += 1

                if star_exp_hit and i < len(desc_tmp) and desc_tmp[i] == "S":
                    star_exp_hit_pos[i] = i       # star hit expression
            elif desc[i] == 'S':
                if sb == -1:
                    sb = i
                    sb_obs = i
                    run = 1
                    while 1:
                        if (i - run) < len(desc_tmp) and desc_tmp[i - run] == "S":
                            star_exp_hit_pos[i - run] = i - \
                                run        # star hit expression
                            sb = i - run
                            run += 1
                        else:
                            break

                try:
                    if not order["s"]:
                        order["s"] = i
                except KeyError:
                    order["s"] = i

                defined = False
                try:
                    _hash[sid]['obs']
                    defined = True
                except KeyError:
                    pass

                if defined and _hash[sid]['obs']:
                    assign_str[i] = col_star_obs
                else:
                    assign_str[i] = col_star_exp

                lstar += 1
            else:
                print_stderr(
                    "something went wrong while parsing alignment in output.mrd file \n")

        pdir = '{}/pdfs_{}'.format(cwd, ltime)
        if not os.path.isdir(pdir):
            os.mkdir(pdir)

        FOLD = open_or_die2('{}/{}.tmp'.format(pdir, filename), 'wb')
        sstruct = _hash[sid]["pri_struct"]
        lstruct_multi = ((len(sstruct) + 2) * multiplier)

        # only for quantitation module
        if mirbase:
            FOLD.write("5{}3\n".format(pri_seq))
        else:    # this is for the moRNA Finder module, where only the assumed precursor without flanks is folded
            if mb < lb:
                err = substr(pri_seq, mb, fl2b - mb)
                FOLD.write("5{}3\n".format(err))
                offset = mb
            else:
                err = substr(pri_seq, sb, fl2b - sb)
                FOLD.write("5{}3\n".format(err))
                offset = sb

        FOLD.close()

        for i in range(0, len(sstruct)):
            hstruct[i] = 0

        for tag in _hash[sid]['reads'].keys():
            for read in sorted(_hash[sid]['reads'][tag].keys()):
                m = re.match(r'^(\.*)(\w+)\.*$',
                             _hash[sid]['reads'][tag][read]['seq'])
                if m:
                    v1 = m.groups()[0]
                    v2 = m.groups()[1]

                    hash2[read] = len(v1)    # begin of read in precursor
                    m = re.search(
                        r'_x(\d+)', _hash[sid]["reads"][tag][read]["rid"])
                    if m:
                        dc = int(m.groups()[0])
                        totalreads += dc

                        create_hash_key_chain(hash2c, 0, tag, v2)
                        hash2c[tag][v2] += dc
                        hash2key[read] = v2
                        hash2order[read] = read
                        # number of mismatches with precursor sequenc
                        hash2mm[read] = _hash[sid]["reads"][tag][read]["mm"]
                        hash2seq[read] = _hash[sid]["reads"][tag][read]["seq"]
                        hash2sample[read] = tag

                        for i in range(len(v1), len(v1) + len(v2)):
                            # saves how often a nt in precursor is covered by a
                            # read
                            hstruct[i] += dc

        y = yorig - downy
        os.chdir('./pdfs_{}'.format(ltime))

        drawCanvas = CreatePDF(_hash, filename)

        DrawStructure(drawCanvas, filename)

        if totalreads != '0':
            CreateHistogram(drawCanvas, _hash)

        CreateAlignment(drawCanvas, _hash)

        y -= 20

        ClosePDF(drawCanvas)

        # ATTR filename_ss.ps does not exist
        # os.unlink("{}/pdfs_{}/{}_ss.ps".format(cwd, ltime, filename))
        os.unlink("{}/pdfs_{}/{}.tmp".format(cwd, ltime, filename))
        os.unlink("{}/pdfs_{}/rna.ps".format(cwd, ltime))
        os.chdir("..")

        print_stderr('creating pdf for {} finished\n'.format(filename))

    try:
        os.unlink('{}/pdfs_{}/tmp'.format(cwd, ltime))
    except:
        pass


if __name__ == '__main__':

    if len(sys.argv) < 2 or re.search(r'-*-he*l*p*', sys.argv[1]):
        Usage()

    opts, args = getopt.getopt(
        sys.argv[1:], "ug:v:f:ck:os:t:er:q:dx:zy:ab:p:V:E")
    options = dict(opts)

    # read in available organisms at the end of the script
    for line in __DATA__.split('\n'):
        line = line.strip()
        if not line:
            continue

        tmp = line
        line = re.sub('\s+', '', line)
        organisms[line] = tmp

    for i in range(0, 200):
        position_hash[counter] = xposshift + multiplier + 20 + i * multiplier
        counter += 1

    if options.get('-u'):
        pprint("\n\nAvailable species organisms were:\n\n")
        for key in organisms.keys():
            pprint(key, "\n")

        pprint("\n\n\n")
        sys.exit(0)

    if options.get('-p'):
        get_precursor_pos()

    if not options.get('-y'):
        die("no timestamp given with parameter y\n")
    else:
        ltime = options.get('-y')

    if options.get('-x') and not options.get('-q'):
        die("\nError:\n\toption -x can only be used together with option -q\n\n")

    cwd = current_dir()

    # organism parameter
    org = organisms[options.get('-t')]

    if options.get('-a') == '' and os.path.isfile("moR_runs/run_{}/tmp/signature.arf".format(ltime)):
        get_mature_pos()

    if options.get('-b'):
        IN = open_or_die2(options.get('-b'), 'rb')
        while True:
            line = IN.readline()
            if not line:
                break

            if re.search(r'\#', line):
                continue

            line = line.strip()
            r = line
            r = re.sub(r'\|', '_', line)
            confident[r] = 1

        IN.close()

    # this option is only set if the quantifier module is used alone!!!
    if options.get('-z') == '':
        PrintKnownnotfound()
        CloseHTML()
        if not options.get('-d') == '':
            mirbase = 1
            CreateStructurePDF(hash_q)

        sys.exit(0)

    # scan program parameters
    if options.get('-f'):
        infile = options.get('-f')
    else:
        print_stderr(
            "Error: no moRNA Finder output file specified by -f option\n\n")
        Usage()

    if not options.get('-t'):
        print_stderr(
            "\nNo organism specified by switch -t [organism] so UCSC browser BLAT will not be used.\n\nTo get a list of possible species names type make_html.py -u \n\nSpecify organism Reads are coming from by  -t [organism] to get links in html file\n")

    if options.get('-v'):
        threshold = options.get('-v')

    # read in known miRNA identifier
    if options.get('-k'):
        KN = open_or_die2(options.get('-k'), 'rb')
        seq1 = None
        while True:
            line = KN.readline()
            if not line:
                break

            m = re.search(r'\>(\S+)', line)
            if m:
                m = m.groups()
                # knownones[m[0]]
                seq1 = KN.readline().lower().strip()
                seq1 = tr(seq1, 't', 'u')
                knownones[m[0]] = seq1

        KN.close()

    # print_stderr(knownones, '\n')

    # create csv file
    if options.get('-c') == '':
        csv = 1
        CSV = open_or_die2("{}/result_{}.csv".format(cwd, ltime), 'w+')

    reads = 0
    line = []

    HTML = open_or_die2('{}/result_{}.html'.format(cwd, ltime), 'w+')

    # Create HTML file
    CreateHTML(HTML)

    # print paramters and files used to html file
    ReadInParameters()

    spacer = None   # length of longest entry
    spaces = None   # string of spaces to fill up spacer

    # Specificity value hash for each score cutoff value
    SP_values = {}

    # read in survey file
    if options.get('-s'):
        pprint("making survey now\n")
        spacer = len("known hairpins (estimated false positives):       ")
        spaces = 0

        # is for extended output
        if options.get('-e') == '':

            p1 = '<th><a href="#" class="tooltip">'
            p11 = '<th><a href="#" class="tooltip2">'

            p2 = '<span>'
            q = '</span></a></th>'

            # hash for mouse overs
            ha = {}
            ii = 1
            for i in range(1, 10):
                ha[i] = {}
                ha[i][1] = None
                ha[i][2] = None

            ha[1][1] = 'moRNA Finder score'
            ha[1][2] = 'for details on how the log-odds score is calculated, see Friedlander et al., Nature Biotechnology, 2008.'
            ha[2][1] = 'predicted by moRNA Finder'
            ha[2][
                2] = 'novel miRNA hairpins are here defined by not having any of the reference mature miRNAs mapping perfectly (full length, no mismatches). The numbers show how many novel miRNA hairpins have a score equal to or exceeding the cut-off.'
            ha[3][1] = 'estimated false positives'
            ha[3][
                2] = 'number of false positive miRNA hairpins predicted at this cut-off, as estimated by the moRNA Finder controls (see Friedlander et al., Nature Biotechnology, 2008). Mean and standard deviation is estimated from 100 rounds of permuted controls.'
            ha[4][1] = 'estimated true positives'
            ha[4][2] = 'the number of true positive miRNA hairpins is estimated as t = total novel miRNAs - false positive novel miRNAs. The percentage of the predicted novel miRNAs that is estimated to be true positives is calculated as p = t / total novel miRNAs. The number of false positives is estimated from 100 rounds of permuted controls. In each of the 100 rounds, t and p are calculated, generating mean and standard deviation of t and p. The variable p can be used as an estimation of moRNA Finder positive predictive value at the score cut-off. '
            ha[5][1] = 'in species'
            ha[5][
                2] = 'number of reference mature miRNAs for that species given as input to moRNA Finder.'
            ha[6][1] = 'in data'
            ha[6][
                2] = 'number of reference mature miRNAs for that species that map perfectly (full length, no mismatches) to one or more of precursor candidates that have been excised from the genome by moRNA Finder.'
            ha[7][1] = 'detected by moRNA Finder'
            ha[7][2] = 'number of reference mature miRNAs for that species that map perfectly (full length, no mismatches) to one or more of predicted miRNA hairpins that have a score equal to or exceeding the cut-off. The percentage of reference mature miRNAs in data that is detected by moRNA Finder is calculated as s = reference mature miRNAs detected / reference mature miRNAs in data. s can be used as an estimation of moRNA Finder sensitivity at the score cut-off.'
            ha[8][1] = 'estimated signal-to-noise'
            ha[8][2] = 'for the given score cut-off, the signal-to-noise ratio is estimated as r = total miRNA hairpins reported / mean estimated false positive miRNA hairpins over 100 rounds of permuted controls.'
            ha[9][1] = 'excision gearing'
            ha[9][2] = 'this is the minimum read stack height required for excising a potential miRNA precursor from the genome in this analysis.'

            # hash finished
            SURVEY = open_or_die(options.get('-s'), 'rb',
                                 'options -s file not found')
            tmp = re.split('\t', SURVEY.readline())
            reduced = 0

            if len(tmp) < 4:
                reduced = 1

            SURVEY.close()

            SURVEY = open_or_die(options.get('-s'), 'rb',
                                 'options -s file not found')
            HTML.write('<font face=\"Times New Roman\" size=\"4\"><b>Survey of moRNA Finder performance for score cut-offs -10 to 10</b>\n<font face=\"Times New Roman\" size=\"3\"><table border=\"1\" >\n')
            i = 0
            sh = {}

            if not reduced:
                HTML.write(
                    "<tr><th></th><TH colspan=\"3\">novel miRNAs</th><th colspan=\"3\">known miRBase miRNAs</th><th></th><th></th></tr>\n")

                for k in sorted(ha.keys()):
                    if k == 1:
                        HTML.write("<th><a href=\"http://www.nature.com/nbt/journal/v26/n4/abs/nbt1394.html\" target=\"blank\" class=\"tooltip\">{}{}{}{}</a>\n\n".format(
                            ha[k][1],
                            p2,
                            ha[k][2],
                            q
                        ))
                    elif k != 9:
                        HTML.write("{}{}{}{}{}\n\n".format(
                            p1,
                            ha[k][1],
                            p2,
                            ha[k][2],
                            q
                        ))
                    else:
                        HTML.write("{}{}{}{}{}\n\n".format(
                            p11,
                            ha[k][1],
                            p2,
                            ha[k][2],
                            q
                        ))
            else:
                HTML.write("<th><a href=\"http://www.nature.com/nbt/journal/v26/n4/abs/nbt1394.html\" target=\"blank\" class=\"tooltip\">{}{}{}{}</a>\n\n".format(
                    ha[1][1],
                    p2,
                    ha[1][2],
                    q
                ))
                HTML.write("{}{}{}{}{}\n\n".format(
                    p1, ha[8][1], p2, ha[8][2], q))
                HTML.write("{}{}{}{}{}\n\n".format(
                    p1, ha[9][1], p2, ha[9][2], q))

            while True:
                l = SURVEY.readline()
                if not l:
                    break

                line = re.split('\t', l)
                if csv and not i:
                    CSV.write(l)

                i += 1
                if i > 1:
                    if csv:
                        CSV.write(l)

                    l = re.sub(r'\t', '<td>', l)
                    if int(line[0]) < 0:
                        sh[i] = "<tr bgcolor=silver><td>{}</td></tr>\n".format(
                            l)
                    else:
                        sh[i] = "<tr><td>{}</td></tr>\n".format(l)

            if csv:
                CSV.write('\n\n\n')

            for k in sorted(sh.keys()):
                sh[k] = re.sub(r'\+\/-', '&#177', sh[k])
                HTML.write(sh[k])

            HTML.write('</table><br><br>\n')
            SURVEY.close()

        SURVEY = open_or_die(options.get('-s'), 'rb',
                             'options -s file not found\n')
        header = re.split('\t', SURVEY.readline())
        survey = []
        std_survey = []

        # also report survey for score cutoff
        if options.get('-g'):
            HTML.write(
                "<font face=\"Times New Roman\" size=\"3\"><table border=\"1\" >\n<th>Survey of run</th>\n")

        while True:
            l = SURVEY.readline()
            if not l:
                break

            survey = re.split(r'\t', l)
            lcol = ''
            m = re.search(r'\((.+)\)', survey[3])
            if m:
                lcol = m.groups()[0]
                lcol = re.sub(r'\+\/\-', '&#177', lcol)

            if survey[0] == 0:
                std_survey = survey

            SP_values[survey[0]] = lcol

            # if line in survey file is reached that hits the threshold then
            # print out stuff
            if survey[0] == threshold:
                for i in range(0, len(survey)):
                    spaces = spacer - (len(header[i]) + 1)
                    if i == 0:
                        # print CSV "Survey of run\t\n$header[$i]
                        # used\t$survey[$i]\n" if($csv);
                        if options.get('-g'):
                            HTML.write(
                                "<tr><td>{} used</td><td>{}</td></tr>\n".fromat(header[i], survey[i]))
                    else:
                        # print CSV "$header[$i]\t$survey[$i]\n" if($csv);
                        if options.get('-g'):
                            HTML.write(
                                "<tr><td>{}</td><td>{}</td></tr>\n".format(header[i], survey[i]))

                expected = SP_values[survey[0]]
                novel = 0
                m = re.search(r'(\d+)\(\d+', survey[6])
                if m:
                    novel = int(m.groups()[0] + 0.5)

                # print CSV "novel hairpins, expected true\t$novel
                # ($expected%)\n" if($csv);
                if options.get('-g'):
                    HTML.write(
                        "<tr><td>novel hairpins, expected true</td><td>{} ({}%)</td></tr>\n".fromat(novel, expected))

        SURVEY.close()
        if threshold == "na":
            threshold = std_survey[0]
            for i in range(0, len(std_survey)):
                spaces = spacer - (len(header[i]) + 1)
                if options.get('-g'):
                    HTML.write(
                        "<tr><td>{}</td><td>{}</td></tr>\n".format(header[i], std_survey[i]))

        if options.get('-g'):
            HTML.write("</table></font>\n")
    else:
        die("Error: no survey file specified with option -s. Please also make sure that your survey file is not empty\n\n")

    # ########################################################
    # #
    # #
    # # here the output.mrd file given by switch -f is parsed
    # #
    # #
    # ########################################################
    IN = open_or_die2(infile, 'rb')
    oid = None
    star_exp_hit = 0
    star_exp_hit_pos = {}

    while True:
        l = IN.readline()
        if not l:
            break

        m = re.match('^\>\s*(\S+)', l)
        if m:
            m = m.groups()
            _id = m[0]
            oid = m[0]
            _id = tr(_id, '|', '_')
            create_hash_key_chain(_hash, oid, _id, 'oid')
            create_hash_key_chain(_hash, _id, _id, 'id')
            counter = 0
            continue

        m = re.match(r'^score_mfe\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'score_mfe')
            continue

        m = re.match(r'^score for star read.+\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'score_star')
            continue

        m = re.match(r'^score for read counts\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'score_read')
            continue

        m = re.match('^score for mfe\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'score_mfe')
            continue

        m = re.match(r'^score total\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'score')
            continue

        m = re.match(r'^total read count\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'freq_total')
            continue

        m = re.match(r'^mature read count\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'freq_mature')
            continue

        m = re.match(r'^loop read count\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'freq_loop')
            continue

        m = re.search('^star read count\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'freq_star')
            continue

        m = re.match(r'^pri_seq\s+(\S+)', l)
        if m:
            create_hash_key_chain(_hash, m.groups()[0], _id, 'pri_seq')
            d = []

            defined = False
            try:
                dummy = _hash[_id]['obs']
                defined = True
            except KeyError:
                pass

            if defined and _hash[_id]['obs']:
                d = ssplit(_hash[_id]['obs'])
            else:
                d = ssplit(_hash[_id]['exp'])

            s = ssplit(_hash[_id]['pri_seq'])
            mseq = ""
            sseq = ""
            lseq = ''

            for i in range(0, len(m.groups()[0])):

                if d[i] != "f":
                    create_hash_key_chain(_hash, '', _id, 'ucsc_seq')
                    _hash[_id]['ucsc_seq'] += s[i]

                if d[i] == "M":
                    mseq += s[i]
                elif d[i] == "S":
                    sseq += s[i]
                elif d[i] == "l":
                    lseq += s[i]

            sseq_obs = ""

            defined = False
            try:
                dummy = _hash[_id]['obs']
                defined = True
            except KeyError:
                pass

            if defined and _hash[_id]['obs']:
                d = ssplit(_hash[_id]['obs'])
                for i in range(0, len(m.groups()[0])):
                    if d[i] == "S":
                        sseq_obs += s[i]

            _hash[_id]['mat_seq'] = mseq
            _hash[_id]['loop_seq'] = lseq
            _hash[_id]["star_seq"] = sseq
            _hash[_id]["star_seq_obs"] = sseq_obs
            continue

        m = re.match(r'^score for randfold\s+(\S+)', l)
        if m:
            _hash[_id]['rand'] = m.groups()[0]
            if float(m.groups()[0]) > 0:
                _hash[_id]["randfold"] = "yes"
            else:
                _hash[_id]["randfold"] = "no"
            continue

        m = re.match(r'^score for cons.\s+seed\s+(\S+)', l)
        if m:
            _hash[_id]["score_cons"] = m.groups()[0]
            if float(m.groups()[0]) > 0:
                _hash[_id]["cons_seed"] = "yes"
            else:
                _hash[_id]["cons_seed"] = ""
            continue

        m = re.match(r'^miRNA with same seed\s+(\S+)', l)
        if m:
            _hash[_id]["cons_seed"] = m.groups()[0]
            continue

        m = re.match(r'^exp\s+(\S+)', l)
        if m:
            _hash[_id]['exp'] = m.groups()[0]
            continue

        m = re.match(r'^obs\s+(\S+)', l)
        if m:
            _hash[_id]['obs'] = m.groups()[0]
            continue

        m = re.match(r'^pri_struct\s+(\S+)', l)
        if m:
            _hash[_id]["pri_struct"] = m.groups()[0]
            reads = 1
            counter = 0
            continue

        m = re.match(r'^(\S\S\S)(\S+)\s+(\S+)\s+(\S+)$', l)
        if m and reads:    # do not read in known miRNAs
            m = m.groups()
            counter += 1

            mm = "{}{}".format(m[0], m[1])
            create_hash_key_chain(_hash, mm, _id, 'reads',
                                  m[0], counter, 'rid')
            create_hash_key_chain(
                _hash, m[2], _id, 'reads', m[0], counter, 'seq')
            create_hash_key_chain(
                _hash, m[3], _id, 'reads', m[0], counter, 'mm')

            _hash[_id]["reads"][m[0]][counter]["rid"] = mm
            _hash[_id]["reads"][m[0]][counter]["seq"] = m[2]
            _hash[_id]["reads"][m[0]][counter]["mm"] = m[3]

            defined = False
            try:
                dummy = knownones[mm]
                defined = True
            except KeyError:
                pass

            if defined and knownones[mm]:
                seen[mm] = 1

            # check here if expected sequence is hit partially by at least one
            # read ## maybe full coverage should be check => later
            if star_exp_hit == 0 and re.search(_hash[_id]["star_seq"], _hash[_id]["reads"][m[0]][counter]["seq"], re.IGNORECASE):
                star_exp_hit = 1

            # open MMU,">>recovered_miRNA";
            if defined and knownones[mm] and reads:

                defined = False
                try:
                    dummy = _hash[_id]['known']
                    defined = True
                except KeyError:
                    pass

                if not (defined and _hash[_id]['known']):
                    _hash[_id]["known"] = mm

                    # if _id == 'CHROMOSOME_I_544':
                    #     print_stderr(_hash[_id]['mat_seq'], '\t', l, '\n')

                    # now check if mature of moRNA Finder is same as mirbase or
                    # not
                    if knownones[mm].find(substr(_hash[_id]['mat_seq'], 4, 12)) >= 0:
                        _hash[_id]['mirbase'] = 'TRUE'
                        _hash[_id]['mirbase-out'] = mm
                        # if star sequence mirbase mature then output star
                    elif knownones[mm].find(substr(_hash[_id]['star_seq'], 4, 12)) >= 0:
                        _hash[_id]['mirbase'] = 'STAR'
                        _hash[_id]['mirbase-out'] = mm
                    else:   # this case should not happen
                        _hash[_id]['mirbase'] = 'NA'
                        _hash[_id]['mirbase-out'] = "NA"

                continue

            if re.match('^\s+$', l) and reads:
                reads = 0
                continue

    IN.close()
    print_stderr("parsing input file finished\n")

    # #########################################################
    #
    #
    # Check precursors for Rfam hits
    #
    #
    # #########################################################

    if options.get('-r'):
        print_stderr("checking Rfam for hits to precursors\n")
        check_Rfam(_hash)
    else:
        print_stderr("no Rfam file specified\tskipping\n")

    percentage = 0
    klist = sorted(SP_values.keys(), reverse=True)
    klen = len(klist)
    blat = None

    # count number of novel miRNAs
    novelc = 0
    for k in _hash.keys():

        defined = False
        try:
            dummy = _hash[k]['known']
            defined = True
        except KeyError:
            pass

        if defined and _hash[k]['known']:
            continue
        novelc += 1

    ##################################################################
    #
    #
    # Print to HTML file the novel miRNAs
    #
    #
    ##################################################################

    if not novelc:
        HTML.write(' <br><h2>no novel miRNAs detected</h2><br>')
    else:
        PrintHtmlTableHeader("novel miRNAs predicted by moRNA Finder")

        if csv:
            PrintHtmlTableHeader("novel miRNAs predicted by moRNA Finder", 1)

        for _id in hash_sort_key(_hash, lambda x: _hash[x[0]]['score'] * -1):

            defined = False
            try:
                dummy = _hash[_id]['known']
                defined = True
            except KeyError:
                pass

            if defined and _hash[_id]['known']:
                continue

            if options.get('-b'):
                if confident[_id] != 1:
                    continue

            if float(_hash[_id]["score"]) >= float(threshold):
                _hash[_id]["pdf"] = 1

                if _hash[_id]["score"] >= klist[0]:
                    percentage = SP_values[klist[0]]
                else:
                    for i in range(0, klen):
                        if _hash[_id]["score"] >= klist[i]:
                            percentage = SP_values[klist[i]]
                            break

                # print to CSV file
                if csv:
                    s_star = _hash[_id]['star_seq']
                    if _hash[_id]['star_seq_obs']:
                        s_star = _hash[_id]['star_seq_obs']
                    p = percentage
                    p = re.sub(r'&#177', '+/-', p, count=1)
                    CSV.write("{}\t{}\t{}\t".format(
                        _hash[_id]['oid'],
                        _hash[_id]['score'],
                        p
                    ))

                    defined = False
                    try:
                        _hash[_id]['rfam']
                        defined = True
                    except KeyError:
                        pass

                    if defined and _hash[_id]['rfam']:
                        CSV.write(_hash[_id]['rfam'])
                        CSV.write('\t')
                    else:
                        CSV.write("-\t")

                    CSV.write("{}\t{}\t{}\t{}\t".format(
                        _hash[_id]['freq_total'],
                        _hash[_id]['freq_mature'],
                        _hash[_id]['freq_loop'],
                        _hash[_id]['freq_star']
                    ))

                    try:
                        if _hash[_id]['randfold']:
                            CSV.write("{}\t-\t".format(_hash[_id]['randfold']))
                        else:
                            CSV.write("-\t-\t")
                    except KeyError:
                        CSV.write("-\t-\t")

                    try:
                        if _hash[_id]['cons_seed']:
                            CSV.write(
                                "{}\t-\t-\t".format(_hash[_id]['cons_seed']))
                        else:
                            CSV.write("-\t-\t-\t")
                    except KeyError:
                        CSV.write("-\t-\t-\t")

                    CSV.write("{}\t{}\t{}".format(
                        _hash[_id]['mat_seq'],
                        s_star,
                        _hash[_id]['ucsc_seq']
                    ))

                    if options.get('-p'):
                        offset = _hash[_id]['pri_seq'].find(
                            _hash[_id]['ucsc_seq'])
                        if pres_coords[_id]['strand'] == '-':
                            CSV.write("\t{}:{}..{}:{}".format(
                                pres_coords[_id]['chr'],
                                pres_coords[_id]['e'] - offset -
                                len(_hash[_id]['ucsc_seq']),
                                pres_coords[_id]['e'] - offset,
                                pres_coords[_id]['strand']
                            ))
                        else:
                            CSV.write("\t{}:{}..{}:{}".format(
                                pres_coords[_id]['chr'],
                                int(pres_coords[_id]['s']) + offset - 1,
                                int(pres_coords[_id]['s']) + offset +
                                len(_hash[_id]['ucsc_seq']) - 1,
                                pres_coords[_id]['strand']
                            ))

                    CSV.write("\n")

                if org == "":
                    blat = '<td></td>'
                else:
                    blat = '''<td><a href="http://genome.ucsc.edu/cgi-bin/hgBlat?org={}&type=BLAT's guess&userSeq={}" target="_blank">blat</a></td>'''.format(org, _hash[
                                                                                                                                                              _id]['ucsc_seq'])

                s_star = _hash[_id]['star_seq']
                if _hash[_id]['star_seq_obs']:
                    s_star = _hash[_id]['star_seq_obs']

                known = "<td nowrap=\"nowrap\"></td>"

                n_count += 1

                science = getEvalue(float(_hash[_id]["score"]), 1)

                m = re.match(r'^(\d+)\s+(\S+)\s+(\d+)\%$', percentage)
                if m:
                    m = m.groups()
                    one = str(float(m[0]) / 100)
                    mid = m[1]
                    two = str(float(m[2]) / 100)

                    if len(str(one)) < 4:
                        one += "0" * (4 - len(one))

                    if re.search(r'0000', one):
                        one = '0.00'

                    if len(two) < 4:
                        two += "0" * (4 - len(two))

                    if re.match(r'0000', two):
                        two = '0.00'

                    percentage = "{} {} {}".format(one, mid, two)

                HTML.write('''<tr><td><a href="pdfs_{}/{}.pdf">{}</a></td>
                <td nowrap="nowrap">{}</td>
                <td>{}</td>
                <td>{}</td>
                <td>{}</td>
                <td>{}</td>
                <td>{}</td>
                <td>{}</td>
                <td>{}</td>
                {}
                <td nowrap="nowrap">{}</td>
                {}
                <td><a href="{}{}&JOB_TITLE={}{}" target="_blank">blast</a></td>
                <td>{}</td>
                <td>{}</td>
                <td>{}</td>'''.format(
                    ltime,
                    _id,
                    _hash[_id]['oid'],
                    science,
                    percentage,
                    _hash[_id].get('rfam', ''),
                    _hash[_id]['freq_total'],
                    _hash[_id]['freq_mature'],
                    _hash[_id]['freq_loop'],
                    _hash[_id]['freq_star'],
                    _hash[_id]['randfold'],
                    known,
                    _hash[_id].get('cons_seed', ''),
                    blat,
                    blast,
                    _hash[_id]['ucsc_seq'],
                    _hash[_id]["oid"],
                    blast_query,
                    _hash[_id]["mat_seq"],
                    s_star,
                    _hash[_id]["ucsc_seq"],
                ))

                if options.get('-p'):
                    offset = _hash[_id]['pri_seq'].find(_hash[_id]['ucsc_seq'])
                    if pres_coords[_id]['strand'] == '-':
                        HTML.write("<td>{}:{}..{}:{}</td>\n".format(
                            pres_coords[_id]['chr'],
                            pres_coords[_id]['e'] - offset -
                            len(_hash[_id]['ucsc_seq']),
                            pres_coords[_id]['e'] - offset,
                            pres_coords[_id]['strand']
                        ))
                    else:
                        HTML.write("<td>{}:{}..{}:{}</td>\n".format(
                            pres_coords[_id]['chr'],
                            int(pres_coords[_id]['s']) + offset - 1,
                            int(pres_coords[_id]['s']) + offset +
                            len(_hash[_id]['ucsc_seq']) - 1,
                            pres_coords[_id]['strand']
                        ))

                HTML.write('</tr>')

        # now print the unconfident ones
        if options.get('-b'):
            HTML.write("<tr><td>--</td></tr>")

            for _id in hash_sort_key(_hash, lambda x: _hash[x[0]]['score'] * -1):
                if _hash[_id]['known']:
                    continue

                if options.get('-b'):
                    if confident[_id] != 1:
                        continue

                if float(_hash[_id]["score"]) >= float(threshold):
                    _hash[_id]["pdf"] = 1

                    if _hash[_id]["score"] >= klist[0]:
                        percentage = SP_values[klist[0]]
                    else:
                        for i in range(0, klen):
                            if _hash[_id]["score"] >= klist[i]:
                                percentage = SP_values[klist[i]]
                                break

                    # print to CSV file
                    if csv:
                        s_star = _hash[_id]['star_seq']
                        if _hash[_id]['star_seq_obs']:
                            s_star = _hash[_id]['star_seq_obs']
                        p = percentage
                        p = re.sub(r'&#177', '+/-', p, count=1)
                        CSV.write("{}\t{}\t{}\t".format(
                            _hash[_id]['oid'],
                            _hash[_id]['score'],
                            p
                        ))

                        if _hash[_id]['rfam']:
                            CSV.write(_hash[_id]['rfam'])
                            CSV.write('\t')
                        else:
                            CSV.write("-\t")

                        CSV.write("{}\t{}\t{}\t{}\t".format(
                            _hash[_id]['freq_total'],
                            _hash[_id]['freq_mature'],
                            _hash[_id]['freq_loop'],
                            _hash[_id]['freq_star']
                        ))

                        if _hash[_id]['randfold']:
                            CSV.write("{}\t-\t".format(_hash[_id]['randfold']))
                        else:
                            CSV.write("-\t-\t")
                        if _hash[_id]['cons_seed']:
                            CSV.write(
                                "{}\t-\t-\t".format(_hash[_id]['cons_seed']))
                        else:
                            CSV.write("-\t-\t-\t")

                        CSV.write("{}\t{}\t{}".format(
                            _hash[_id]['mat_seq'],
                            s_star,
                            _hash[_id]['ucsc_seq']
                        ))

                        if options.get('-p'):
                            offset = _hash[_id]['pri_seq'].find(
                                _hash[_id]['ucsc_seq'])
                            if pres_coords[_id]['strand'] == '-':
                                CSV.write("\t{}:{}..{}:{}".format(
                                    pres_coords[_id]['chr'],
                                    pres_coords[_id]['e'] - offset -
                                    len(_hash[_id]['ucsc_seq']),
                                    pres_coords[_id]['e'] - offset,
                                    pres_coords[_id]['strand']
                                ))
                            else:
                                CSV.write("\t{}:{}..{}:{}".format(
                                    pres_coords[_id]['chr'],
                                    int(pres_coords[_id]['s']) + offset - 1,
                                    int(pres_coords[_id][
                                        's']) + offset + len(_hash[_id]['ucsc_seq']) - 1,
                                    pres_coords[_id]['strand']
                                ))

                        CSV.write("\n")

                    if org == "":
                        blat = '<td></td>'
                    else:
                        blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org={}&type=BLAT's guess&userSeq={}\" target=\"_blank\">blat</a></td>".format(org, _hash[
                                                                                                                                                                  _id]['ucsc_seq'])

                    s_star = _hash[_id]['star_seq']
                    if _hash[_id]['star_seq_obs']:
                        s_star = _hash[_id]['star_seq_obs']

                    known = "<td nowrap=\"nowrap\"></td>"

                    n_count += 1

                    science = getEvalue(float(_hash[_id]["score"]), 1)

                    m = re.match(r'^(\d+)\s+(\S+)\s+(\d+)\%$', percentage)
                    if m:
                        m = m.groups()
                        one = str(float(m[0]) / 100)
                        mid = m[1]
                        two = str(float(m[2]) / 100)

                        if len(str(one)) < 4:
                            one += "0" * (4 - len(one))

                        if re.search(r'0000', one):
                            one = '0.00'

                        if len(two) < 4:
                            two += "0" * (4 - len(two))

                        if re.match(r'0000', two):
                            two = '0.00'

                        percentage = "{} {} {}".format(one, mid, two)

                    HTML.write('''<tr><td><a href="pdfs_{}/{}.pdf">{}</a></td>
                    <td nowrap="nowrap">{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    {}
                    <td nowrap="nowrap">{}</td>
                    {}
                    <td><a href="{}{}&JOB_TITLE={}{}" target="_blank">blast</a></td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>'''.format(
                        ltime,
                        _id,
                        _hash[_id]['oid'],
                        science,
                        percentage,
                        _hash[_id]['rfam'],
                        _hash[_id]['freq_total'],
                        _hash[_id]['freq_mature'],
                        _hash[_id]['freq_loop'],
                        _hash[_id]['freq_star'],
                        _hash[_id]['randfold'],
                        known,
                        _hash[_id]['cons_seed'],
                        blat,
                        blast,
                        _hash[_id]['ucsc_seq'],
                        _hash[_id]["oid"],
                        blast_query,
                        _hash[_id]["mat_seq"],
                        s_star,
                        _hash[_id]["ucsc_seq"],
                    ))

                    if options.get('-p'):
                        offset = _hash[_id]['pri_seq'].find(
                            _hash[_id]['ucsc_seq'])
                        if pres_coords[_id]['strand'] == '-':
                            HTML.write("<td>{}:{}..{}:{}</td>\n".format(
                                pres_coords[_id]['chr'],
                                pres_coords[_id]['s'],
                                pres_coords[_id]['e'],
                                pres_coords[_id]['strand']
                            ))
                        else:
                            HTML.write("<td>na</td>\n")

                    HTML.write('</tr>')

        HTML.write("</table></font>")

    # #################################################################
    #
    #
    # Print to HTML file the known miRNAs
    #
    #
    # #################################################################
    if options.get('-k'):
        PrintHtmlTableHeader("mature miRBase miRNAs detected by moRNA Finder")

        if csv:
            CSV.write('\n\n\n')
            PrintHtmlTableHeader(
                "mature miRBase miRNAs detected by moRNA Finder", 1)

        for _id in hash_sort_key(_hash, lambda x: _hash[x[0]]['score'] * -1):

            defined = False
            try:
                dummy = _hash[_id]['known']
                defined = True
            except KeyError:
                pass

            if not (defined and _hash[_id]['known']):
                continue

            _hash[_id]["pdf"] = 1

            if _hash[_id]["score"] >= klist[0]:
                percentage = SP_values[klist[0]]
            else:
                for i in range(0, klen):
                    if _hash[_id]["score"] >= klist[i]:
                        percentage = SP_values[klist[i]]
                        break

            if _hash[_id]["known"]:
                known = "<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms={} target=\"_blank\"> {}</a></td>".format(
                    _hash[_id]["known"],
                    _hash[_id]["known"]
                )
            else:
                known = "<td nowrap=\"nowrap\"></td>"

            if csv:
                s_star = _hash[_id]['star_seq']
                if _hash[_id]['star_seq_obs']:
                    s_star = _hash[_id]['star_seq_obs']
                p = percentage
                p = re.sub(r'&#177', '+/-', p, count=1)

                if not novelc:
                    CSV.write("{}\t{}\t-\t".format(
                        _hash[_id]['oid'],
                        _hash[_id]['score']
                    ))
                else:
                    CSV.write("{}\t{}\t{}\t".format(
                        _hash[_id]['oid'],
                        _hash[_id]['score'],
                        p
                    ))

                defined = False
                try:
                    _hash[_id]['rfam']
                    defined = True
                except KeyError:
                    pass

                if defined and _hash[_id]['rfam']:
                    CSV.write("{}\t".format(_hash[_id]['rfam']))
                else:
                    CSV.write("-\t")

                CSV.write('{}\t{}\t{}\t{}\t'.format(
                    _hash[_id]['freq_total'],
                    _hash[_id]['freq_mature'],
                    _hash[_id]['freq_loop'],
                    _hash[_id]['freq_star']
                ))

                try:
                    if _hash[_id]['randfold']:
                        CSV.write("{}\t".format(_hash[_id]['randfold']))
                    else:
                        CSV.write('-\t')
                except KeyError:
                    CSV.write('-\t')

                try:
                    if _hash[_id]['known']:
                        CSV.write("{}\t".format(_hash[_id]['known']))
                    else:
                        CSV.write('-\t')
                except KeyError:
                    CSV.write('-\t')

                try:
                    if _hash[_id]['cons_seed']:
                        CSV.write("{}\t-\t-\t".format(_hash[_id]['cons_seed']))
                    else:
                        CSV.write("-\t-\t-\t")
                except KeyError:
                    CSV.write("-\t-\t-\t")

                CSV.write("{}\t{}\t{}".format(
                    _hash[_id]['mat_seq'],
                    s_star,
                    _hash[_id]['ucsc_seq']
                ))

                if options.get('-p'):
                    offset = _hash[_id]['pri_seq'].find(_hash[_id]['ucsc_seq'])
                    if pres_coords[_id]['strand'] == '-':
                        CSV.write("\t{}:{}..{}:{}".format(
                            pres_coords[_id]['chr'],
                            pres_coords[_id]['e'] - offset -
                            len(_hash[_id]['ucsc_seq']),
                            pres_coords[_id]['e'] - offset,
                            pres_coords[_id]['strand']
                        ))
                    else:
                        CSV.write("\t{}:{}..{}:{}".format(
                            pres_coords[_id]['chr'],
                            int(pres_coords[_id]['s']) + offset - 1,
                            int(pres_coords[_id]['s']) + offset +
                            len(_hash[_id]['ucsc_seq']) - 1,
                            pres_coords[_id]['strand']
                        ))

                CSV.write('\n')

            if org == "":
                blat = '<td></td>'
            else:
                blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org={}&type=BLAT's guess&userSeq={}\" target=\"_blank\">blat</a></td>".format(org, _hash[
                                                                                                                                                          _id]['ucsc_seq'])

            s_star = _hash[_id]['star_seq']
            if _hash[_id]['star_seq_obs']:
                s_star = _hash[_id]['star_seq_obs']

            n_count += 1

            science = getEvalue(float(_hash[_id]["score"]), 1)

            if not novelc:
                percentage = '-'

            m = re.match(r'^(\d+)\s+(.+)\s+(\d+)\%$', percentage)
            if m:
                m = m.groups()
                one = str(float(m[0]) / 100)
                mid = m[1]
                two = str(float(m[2]) / 100)

                if len(str(one)) < 4:
                    one += "0" * (4 - len(one))

                if re.search(r'0000', one):
                    one = '0.00'

                if len(two) < 4:
                    two += "0" * (4 - len(two))

                if re.match(r'0000', two):
                    two = '0.00'

                percentage = "{} {} {}".format(one, mid, two)

            HTML.write('''<tr><td><a href="pdfs_{}/{}.pdf">{}</a></td>
            <td nowrap="nowrap">{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>
            {}
            <td nowrap="nowrap">{}</td>
            {}
            <td><a href="{}{}&JOB_TITLE={}{}" target="_blank">blast</a></td>
            <td>{}</td>
            <td>{}</td>
            <td>{}</td>'''.format(
                ltime,
                _id,
                _hash[_id]['oid'],
                science,
                percentage,
                _hash[_id].get('rfam', ''),
                _hash[_id]['mirbase'],
                _hash[_id]['freq_total'],
                _hash[_id]['freq_mature'],
                _hash[_id]['freq_loop'],
                _hash[_id]['freq_star'],
                _hash[_id]['randfold'],
                known,
                _hash[_id].get('cons_seed', ''),
                blat,
                blast,
                _hash[_id]['ucsc_seq'],
                _hash[_id]["oid"],
                blast_query,
                _hash[_id]["mat_seq"],
                s_star,
                _hash[_id]["ucsc_seq"],
            ))

            if options.get('-p'):
                offset = _hash[_id]['pri_seq'].find(_hash[_id]['ucsc_seq'])
                if pres_coords[_id]['strand'] == '-':
                    HTML.write("<td>{}:{}..{}:{}</td>\n".format(
                        pres_coords[_id]['chr'],
                        pres_coords[_id]['e'] - offset -
                        len(_hash[_id]['ucsc_seq']),
                        pres_coords[_id]['e'] - offset,
                        pres_coords[_id]['strand']
                    ))
                else:
                    HTML.write("<td>{}:{}..{}:{}</td>\n".format(
                        pres_coords[_id]['chr'],
                        int(pres_coords[_id]['s']) + offset - 1,
                        int(pres_coords[_id]['s']) + offset +
                        len(_hash[_id]['ucsc_seq']) - 1,
                        pres_coords[_id]['strand']
                    ))
            else:
                HTML.write('<td>na</td>')

        HTML.write('</tr>')

    if csv:
        CSV.close()

    IN.close()
    HTML.write("</table></font>")

    # #################################################################
    #
    #
    # Print to HTML file the miRNAs in data but not detected by moRNA Finder
    #
    #
    # #################################################################
    if options.get('-q'):
        PrintKnownnotfound()
    else:
        if csv:
            CSV.close()
        IN.close()

    CloseHTML(HTML)

    # #################################################################
    #
    #
    # Create PDF files for all entrys in the HTML file
    #
    #
    # #################################################################

    if not options.get('-d'):
        CreateStructurePDF(_hash)
        if options.get('-q'):
            mirbase = 1
            CreateStructurePDF(hash_q)

    sys.exit(0)
