#!/usr/bin/env python
from __future__ import print_function

import getopt
import math
import os
import re
import sys

from port import (create_hash_key_chain, die, esplit, file_s, fileparse,
                  hash_defined, hash_sort_key, hash_value_true, open_or_die,
                  open_or_die2, print_stderr, ssplit, tr)
from reportlab.lib.colors import (black, darkviolet, grey, lightskyblue,
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
mat_pre_arf = {}
lflank1 = None  # length(string of left flank)
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
star_exp_hit_pos = {}
blat = None
spacer = None     # length of longest entry
spaces = None     # string of spaces to fill up spacer
# stores begin coordinates of fl1,m,l,s,fl2 sequences
order = None
multiplier = 3.6  # 4.825                # minimal distance between two letters
# calculate predefined pdf loci for alignment characters
position_hash = {}
counter = 0
yorig = 500  # 500
downy = 50
dline = None                             # line graphic handler
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

HTML = None


def current_dir():
    '''
    Get current dir path
    '''
    return os.path.realpath('.')


def Usage():
    print_stderr(
        "\n\n\n\n[usage]\n\tpython make_html.py -f moR_outfile [options]\n\n")
    print_stderr("[options]\n-v 2\t only output hairpins with score above 2\n")
    print_stderr("-c\t also create overview in excel format.\n")
    print_stderr("-k file\t supply file with known miRNAs\n")
    print_stderr("-s file\t supply survey file if score cutoff is used to get information about how big is the confidence of resulting reads\n\n\n")
    print_stderr("-e \t report complete survey file\n")
    print_stderr("-g \t report survey for current score cutoff\n")
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
    print_stderr(
        "-l \tbe stringent when assigning miRNA-precursor connections like mmu-mir only is assigned to mmu-precursor\n")
    sys.exit(0)


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

        m = re.findall(r'>(\S+)')
        if m:
            mid = m[0]

        m = re.findall(r'exp\s+(\S+)', line)
        if m:
            seq = m[0]
            start = seq.find('M')
            end = seq.rfind('M')

        m = re.findall(r'^(\S+_x\d+)\s+(\S+)', line)
        if m:
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


def PrintHtmlTableHeader(hl='', csv=''):
    global HTML
    h = {}

    # divide string by linebreaks every x characters
    p1 = '<th><a href="" class="tooltip3">'
    p11 = '<th><a href="" class="tooltip4">'
    p2 = '<span>'
    q = '</span></a></th>'

    if re.search('novel', hl, re.IGNORECASE):
        if options.get('-a') == '':
            end = 18
        else:
            end = 17
        for i in range(1, end):
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
        if options.get('-a') == '':
            h[17][1] = 'genomic position'

    elif re.search(r'miRBase miRNAs in dataset', hl, re.IGNORECASE):
        if options.get('-a') == '':
            end = 18
        else:
            end = 17
        for i in range(1, end):
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
        h[10][1] = 'miRBase miRNA'
        h[10][2] = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA map to the miRNA hairpin, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs for the species is displayed. Clicking this field will link to miRBase.'
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

        if options.get('-a') == '':
            h[17][1] = 'genomic position'

    else:
        if options.get('-a') == '':
            end = 18
        else:
            end = 17
        for i in range(1, end):
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

        if options.get('-P'):
            h[5][2] = 'this is the sum of read counts for the 5p and 3p sequences.'
            h[6][1] = '5p read counts'
            h[6][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 5p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.'
            h[8][1] = '3p read counts'
            h[8][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 3p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.'
            h[9][1] = 'remaining reads'
            h[9][2] = 'this is the number of reads that did not map to any of the 5p or 3p sequences'
        else:
            h[5][2] = 'this is the sum of read counts for the mature and star miRNAs.'
            h[6][1] = 'mature read count(s)'
            h[6][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the mature miRNA, including 2 nts upstream and 5 nts downstream. If more than one mature sequence is given this will be a comma separated list. In parenthesis are normalized read counts shown.'
            h[8][1] = 'star read count'
            h[8][2] = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the star miRNA, including 2 nts upstream and 5 nts downstream. This field is empty unless a reference star miRNA was given as input to quantifier.py. If more than one mature sequence is given this will be a comman separated list'
            h[9][1] = 'remaining reads'
            h[9][2] = 'this is the number of reads that did not map to any of the mature and star sequences'

        h[7][1] = '-    '
        h[7][2] = '-    '
        h[10][1] = '-'   # 'miRBase mature id'
        h[10][2] = '-'   # 'Clicking this field will link to miRBase.'
        h[11][1] = '-'
        h[11][2] = '-'
        h[12][1] = 'UCSC browser'
        h[12][2] = 'if a species name was input to moRNA Finder, then clicking this field will initiate a UCSC blat search of the miRNA precursor sequence against the reference genome.'
        h[13][1] = 'NCBI blastn'
        h[13][
            2] = 'clicking this field will initiate a NCBI blastn search of the miRNA precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).'

        if options.get('-P'):
            h[14][1] = 'miRBase 5p sequence(s)'
            h[14][
                2] = 'this is/are the 5p miRNA sequence(s) input to quantifier.py.'
            h[15][1] = 'miRBase 3p sequence(s)'
            h[15][
                2] = 'this is/are the 3p miRNA sequence(s) input to quantifier.py.'
        else:
            h[14][1] = 'miRBase mature sequence(s)'
            h[14][
                2] = 'this is/are the mature miRNA sequence(s) input to quantifier.py.'
            h[15][1] = 'miRBase star sequence(s)'
            h[15][
                2] = 'this is/are the star miRNA sequence(s) input to quantifier.py. This field is empty unless a reference star miRNA was given as input to quantifier.py.'

        h[16][1] = 'miRBase precursor sequence'
        h[16][2] = 'this is the precursor miRNA sequence input to quantifier.py.'
        if options.get('-a') == '':
            h[17][1] = 'genomic position'

    if csv:
        f = 1
        CSV.write('{}\n'.format(hl))
        for k in hash_sort_key(h, lambda x: x[0]):
            if f:
                f = 0
                CSV.write(h[k][1])
            else:
                CSV.write('\t{}'.format(h[k][1]))
        CSV.write('\n')
        return

    HTML.write('''
        <br>
        <br>
        <br><h2>{}</h2><br>
        <font face="Times New Roman" size="2">
        <table border="1">
    '''.format(hl))

    for k in hash_sort_key(h, lambda x: x[0]):
        if k == 2 or k == 3 or k == 4 or k == 7 or k == 10 or k == 11:
            continue

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
        else:
            HTML.write('{}{}{}{}{}\n\n'.format(
                p11,
                h[k][1],
                p2,
                h[k][2],
                q
            ))


def CloseHTML():
    HTML.write('''
        </table>
        </body>
        </html>
    ''')
    HTML.close()


def PrintQuantifier():
    global known, sig, HTML, counter, ltime, mat_pre_arf, hairpin2mature, mature2hairpin, files_mirnaex, hp, hm, hs
    not_seen = {}
    signature = {}
    reads = None

    # create HTML for quantifier module
    fname = 'expression_analyses/expression_analyses_{}/expression_{}.html'.format(
        ltime, ltime)
    HTML = open_or_die2(fname, 'w+')
    CreateHTML(HTML)
    PrintHtmlTableHeader()

    # now read in the mature_ref_this_species mapped against precursors from the quantifier module
    # store ids in hashes

    (mature, path0, extension0) = fileparse(options.get('-k'), '\..*')
    # filename

    mature += '_mapped.arf'  # change .fa suffix to _mapped.arf

    exprs, exprs2 = {}, {}
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

        create_hash_key_chain(hairpin2mature, '', line[2])
        hairpin2mature[line[2]] += line[0] + ","

        exprs[line[0]] = line[1]

        create_hash_key_chain(exprs2, '', line[2], line[0])
        exprs2[line[2]][line[0]] = line[1]

    IN.close()

    exprs_sample = {}

    # read in sample stuff to expression_value hash
    for f in files_mirnaex:
        IN = open_or_die2(f, 'rb')
        samples = []
        while True:
            l = IN.readline()
            if not l:
                break

            line = re.split('\t', l)
            if re.search(r'precursor', l):
                samples = line

            for index in range(4, len(line)):
                create_hash_key_chain(exprs_sample, line[index], samples[
                                      index], line[2], line[0])
                exprs_sample[samples[index]][line[2]][line[0]] = line[index]

        IN.close()

    exdir = 'expression_analyses/expression_analyses_{}'.format(ltime)

    # read in mature,precursor and star seuqnces;
    fname = '{}/mature.converted'.format(exdir)
    IN = open_or_die(fname, 'rb', 'file mature.converted not found\n')

    seq, _id = None, None
    while True:
        line = IN.readline()
        if not line:
            break

        m = re.findall(r'>(\S+)', line)
        if m:
            _id = m[0]
            seq = IN.readline().strip().lower()
            seq = tr(seq, 't', 'u')

            if options.get('-P') == '':
                if re.search(r'-5p', _id):
                    hm[_id] = seq
                elif re.search(r'-3p', line):
                    hm[_id] = seq
                else:
                    hs[_id] = seq
                    hm[_id] = seq
            else:
                hm[_id] = seq

    width = 0
    for k in hm.keys():
        if len(k) > width:
            width = len(k)

    width *= 5

    IN.close()

    fname = '{}/star.converted'.format(exdir)
    if os.path.isfile(fname):
        IN = open_or_die2(fname, 'rb')
        while True:
            line = IN.readline()
            if not line:
                break

            m = re.findall(r'>(\S+)', line)
            if m:
                _id = m[0]
                seq = IN.readline().strip().lower()
                seq = tr(seq, 't', 'u')
                hs[_id] = seq

        IN.close()

    fname = '{}/precursor.converted'.format(exdir)
    IN = open_or_die2(fname, 'rb')
    while True:
        line = IN.readline()
        if not line:
            break

        m = re.findall(r'>(\S+)', line)
        if m:
            _id = m[0]
            seq = IN.readline().strip()
            hp[_id] = seq
    IN.close()

    fname = 'expression_analyses/expression_analyses_{}/mature2hairpin'.format(
        ltime)
    IN = open_or_die2(fname, 'rb')
    hairpin2mature2 = {}
    while True:
        l = IN.readline().strip()
        if not l:
            break

        line = re.split(r'\t', l)
        hairpin2mature2[line[0]] = line[1]
    IN.close()

    if not options.get('-i'):
        die('options -i file is not specified.')

    IN = open_or_die2(options.get('-i'), 'rb')
    while True:
        l = IN.readline().strip()
        if not l:
            break

        line = re.split('\t', l)
        id1h = line[0]  # mature ID
        id2h = line[5]  # precursor id

        id1h = re.sub('-5p', '', id1h)
        id1h = re.sub('-3p', '', id1h)

        # here is assumed that multiple precursor ids have 3 - in their id,
        # seems to be ok so far
        m = re.match('^(\w+\-\w+\-\w+)\-\d+$', id2h)
        if m:
            id2h = m.groups()[0]

        # stringent mapping let7a only allowed to map pre-let7a if k is given
        if options.get('-l') == '' and not re.search(id2h, id1h, re.IGNORECASE) and not re.search(id1h, id2h, re.IGNORECASE):
            continue

        create_hash_key_chain(mat_pre_arf, '', line[5], line[0], 'start')
        mat_pre_arf[line[5]][line[0]]['start'] = int(line[7])

        create_hash_key_chain(mat_pre_arf, '', line[5], line[0], 'end')
        mat_pre_arf[line[5]][line[0]]['end'] = int(line[8])

        create_hash_key_chain(mat_pre_arf, '', line[5], line[0], 'seq')
        mat_pre_arf[line[5]][line[0]]['seq'] = line[9]

        create_hash_key_chain(mat_pre_arf, None, line[5], line[0], 'star')
        mat_pre_arf[line[5]][line[0]][
            'star'] = 1 if re.search(r'\*', line[0]) else 0

    IN.close()

    # think about a stringent option that designates if to show only mappings where the ID of precursor and mature is the same !!!
    # read in star mapped to mature sequence if a star file was given
    if options.get('-j') and os.path.isfile(options.get('-j')):
        IN = open_or_die2(options.get('-j'), 'rb')
        while True:
            l = IN.readline().strip()
            if not l:
                break

            line = re.split(r"\t", l)
            id1h = line[0]  # this is the mature ID
            id2h = line[5]  # this is the precursor ID

            # remove multiple endings if ambigous just for matching with
            # precursor
            id1h = re.sub(r'-5p', '', l)
            id1h = re.sub(r'-3p', '', l)

            # stringent mapping let7a only allowed to map pre-let7a if k is
            # given
            if options.get('-l') == '' and re.search(id2h, id1h, re.IGNORECASE) and re.search(id1h, id2h, re.IGNORECASE):
                continue

            create_hash_key_chain(mat_pre_arf, '', line[5], line[0], 'start')
            mat_pre_arf[line[5]][line[0]]['start'] = line[7]

            create_hash_key_chain(mat_pre_arf, '', line[5], line[0], 'end')
            mat_pre_arf[line[5]][line[0]]['end'] = line[8]

            create_hash_key_chain(mat_pre_arf, '', line[5], line[0], 'seq')
            mat_pre_arf[line[5]][line[0]]['seq'] = line[9]

            create_hash_key_chain(mat_pre_arf, None, line[5], line[0], 'star')
            if re.search(r'\*', line[0]):
                mat_pre_arf[line[5]][line[0]]['star'] = 1

                create_hash_key_chain(mat_pre_arf, None, line[5], 'sb')
                if not mat_pre_arf[line[5]]['sb']:
                    mat_pre_arf[line[5]]['sb'] = line[7]
            else:
                mat_pre_arf[line[5]][line[0]]['star'] = 0

                create_hash_key_chain(mat_pre_arf, None, line[5], 'sb')
                if not mat_pre_arf[line[5]]['mb']:
                    mat_pre_arf[line[5]]['mb'] = line[7]

        IN.close()

    # open miRBase.mrd from quantifier module ## precursor ids as hash
    oid = None
    if not options.get('-q'):
        die('option -q is not specified.')

    IN = open_or_die2(options.get('-q'), 'rb')
    while True:
        l = IN.readline()
        if not l:
            break

        m = re.match(r'^\>(\S+)', l)
        if m:
            _id = m.groups()[0]
            oid = m.groups()[0]
            _id = tr(_id, '|', '_')

            create_hash_key_chain(hash_q, None, _id, 'oid')
            hash_q[_id]['oid'] = oid

            create_hash_key_chain(hash_q, None, _id, 'id')
            hash_q[_id]['id'] = _id

            counter = 0
        else:
            m = re.match(r'^remaining read count\s*(\d+)', l)
            if m:
                create_hash_key_chain(hash_q, None, _id, 'remaining_rc')
                hash_q[_id]['remaining_rc'] = m.groups()[0]
            else:
                m = re.match(r'^total read count\s*(\d*)', l)
                if m:
                    create_hash_key_chain(hash_q, None, _id, 'freq_total')
                    hash_q[_id]["freq_total"] = int(m.groups()[0])
                else:
                    # read in here everything that mapped to the precursor
                    m = re.match(r'^(\S+) read count\s*(\d+)', l)
                    if m:
                        m = m.groups()
                        create_hash_key_chain(
                            hash_q, None, _id, 'mapped', m[0])
                        hash_q[_id]['mapped'][m[0]] = m[1]
                    else:
                        m = re.match(r'^pri_seq\s+(\S+)', l)
                        if m:
                            create_hash_key_chain(hash_q, None, _id, 'pri_seq')
                            hash_q[_id]['pri_seq'] = m.groups()[0]
                            d = {}
                            if hash_value_true(hash_q, _id, 'obs'):
                                d = ssplit(hash_q[_id]['obs'])
                            else:
                                d = ssplit(hash_q[_id]['exp'])

                            s = ssplit(hash_q[_id]['pri_seq'])
                            mseq = ''
                            sseq = ''

                            # this can also be done by just using the information of the mature_signature_arf file
                            # now set labels for mature and star seq in
                            # explanation string
                            for i in range(0, len(m.groups()[0])):
                                if d[i] != 'f':
                                    create_hash_key_chain(
                                        hash_q, '', _id, 'ucsc_seq')
                                    hash_q[_id]['ucsc_seq'] += s[i]
                                if d[i] == 'M' or d[i] == '5' or d[i] == '3':
                                    mseq += s[i]
                                elif d[i] == 'S':
                                    sseq += s[i]

                            # if there is an observed star sequence use this
                            sseq_obs = ""
                            if hash_value_true(hash_q, _id, 'obs'):
                                d = ssplit(hash_q[_id]['obs'])
                                for i in range(0, len[m[0]]):
                                    if d[i] == 'S':
                                        sseq_obs += s[i]

                            hash_q[_id]["mat_seq"] = mseq
                            hash_q[_id]["star_seq"] = sseq
                            hash_q[_id]["star_seq_obs"] = sseq_obs
                        else:
                            m = re.findall(r'^exp\s+(\S+)', l)
                            if m:
                                create_hash_key_chain(hash_q, None, _id, 'exp')
                                hash_q[_id]['exp'] = m[0]
                            else:
                                m = re.findall(r'^obs\s+(\S+)', l)
                                if m:
                                    create_hash_key_chain(
                                        hash_q, None, _id, 'obs')
                                    hash_q[_id]['obs'] = m[0]
                                else:
                                    m = re.findall(r'^pri_struct\s+(\S+)', l)
                                    if m:
                                        hash_q[_id]["pri_struct"] = m[0]
                                        reads = 1
                                        counter = 0

                                        continue
                                    else:
                                        m = re.findall(
                                            r'^(\S\S\S)(\S+)\s+(\S+)\s+(\S+)$', l)
                                        if m and reads:
                                            m = m[0]
                                            counter += 1
                                            create_hash_key_chain(
                                                hash_q, None, _id, 'reads', m[0], counter, 'rid')
                                            hash_q[_id]["reads"][m[0]][counter][
                                                "rid"] = "{}{}".format(m[0], m[1])
                                            create_hash_key_chain(
                                                hash_q, None, _id, 'reads', m[0], counter, 'seq')
                                            hash_q[_id]["reads"][m[0]][
                                                counter]["seq"] = m[2]
                                            create_hash_key_chain(
                                                hash_q, None, _id, 'reads', m[0], counter, 'mm')
                                            hash_q[_id]["reads"][m[0]][
                                                counter]["mm"] = m[3]
    IN.close()
    print_stderr('parsing miRBase.mrd file finished\n')

    sf = None
    mf = None

    for _id in hash_sort_key(exprs2, lambda x: hash_q[x[0]]['freq_total'] * -1):
        mature = {}
        star = {}
        ind = None
        s_star = hash_q[_id]['star_seq']

        for k2 in exprs2[_id].keys():  # k1 = precursor, k2 = mirna mapped to it
            if options.get('-P') == '':
                if re.search(r'-3p', k2):
                    star[k2] = exprs2[_id][k2]
                elif re.search(r'-5p', k2):
                    mature[k2] = exprs2[_id][k2]
                else:
                    ind = hash_q[_id]["pri_seq"].find(hash_q[_id]["mat_seq"])

                    if ind > len(hash_q[_id]["pri_seq"]) / 2:
                        star[k2] = exprs2[_id][k2]
                    elif ind > 0:
                        mature[k2] = exprs2[_id][k2]
                    else:
                        print_stderr(
                            'Could not determine where {} sits in precursor\nPutting it to the 5p species hash\n'.format(k2))
                        mature[k2] = exprs2[_id][k2]
            else:
                if re.search(r'\*', k2):
                    star[k2] = exprs2[_id][k2]
                else:
                    mature[k2] = exprs2[_id][k2]

        hash_q[_id]["pdf"] = 1
        sig += 1

        # set blat link for pre-miRNA sequence
        if org == "":
            blat = '<td> - </td>'
        else:
            blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq={}\" target=\"_blank\">blat</a></td>".format(hash_q[
                                                                                                                                                        _id]['pri_seq'])

        if (hash_q[_id]['star_seq_obs']):
            s_star = hash_q[_id]['star_seq_obs']

        s_mat = hash_q[_id]['mat_seq']

        # here the belonging precursor is shown
        known = "<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms={}\" target=\"_blank\">{}</a></td>".format(
            _id, oid)

        try:
            sf = hash_q[_id]["freq_star"]
        except KeyError:
            sf = ''

        try:
            mf = hash_q[_id]["freq_mature"]
        except KeyError:
            mf = ''

        try:
            if mf != exprs[oid]:
                mf = exprs[oid]
                if sf:
                    sf = mf
        except KeyError:
            pass

        # do not print precursors with 0 reads mapped to it
        if not hash_value_true(hash_q, _id, 'freq_total'):
            continue

        if options.get('-d') == '' or not hash_value_true(hash_q, _id, 'freq_total'):
            HTML.write(
                '<tr><td nowrap=\"nowrap\"><div style=\"visibility: hidden\">{}.pdf </div>{}</a></td>'.format(_id, _id))
        else:
            HTML.write(
                '<tr><td nowrap=\"nowrap\"><a href=\"pdfs_{}/{}.pdf\">{}</a></td>'.format(ltime, _id, _id))

        HTML.write(
            '\n<td>{}</td>\n<td>\n\n\n'.format(hash_q[_id]["freq_total"]))

        mf = "<table>"
        s_mat = "<table>"

        mf += "<tr><td WIDTH={}>id</td>".format(width)
        for sample in sorted(exprs_sample.keys()):
            mf += "<td>{}</td>".format(sample)

        if len(mature.keys()) != 0:
            for mx in sorted(mature.keys()):
                mf += "\n<tr><td  nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.py?terms={}\" target=\"_blank\">{}</a></td>".format(
                    _id, mx)
                for sample in sorted(exprs_sample.keys()):
                    mf += "<td><nobr>{}</nobr></td>".format(
                        exprs_sample[sample][_id][mx].strip())

                mf += '</tr>\n'
                s_mat += '<tr><td>{}</td></tr>'.format(hm[mx])
        else:
            mf += "\n<tr><td  nowrap=\"nowrap\"><a>na</a></td>"
            for sample in sorted(exprs_sample.keys()):
                mf += "<td><nobr>0</nobr></td>"

            mf += "</tr>\n"
            s_mat += "<tr><td>-</td></tr>"

        mf += "</table>\n"
        s_mat += "</table>"

        HTML.write('\n{}\n</td>'.format(mf))

        s_star = "<table>"
        sf = "<table>"
        sf += "<tr><td WIDTH={}>id</td>".format(width)

        for sample in sorted(exprs_sample.keys()):
            sf += "<td>{}</td>".format(sample)

        if len(star.keys()) != 0:
            for sx in star.keys():
                sf += "\n<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms={}\" target=\"_blank\">{}</a></td>".format(
                    _id, sx)

                for sample in sorted(exprs_sample.keys()):
                    sf += "<td><nobr>{}</nobr></td>".format(
                        exprs_sample[sample][_id][sx])

                sf += "</tr>\n"

                s_star += "<tr><td>{}</td>\n</tr>".format(hs[sx])
        else:
            sf += "\n<tr><td  nowrap=\"nowrap\"><a>na</a></td>"
            for sample in sorted(exprs_sample.keys()):
                sf += "<td><nobr>0</nobr></td>"

            sf += "</tr>\n"
            s_star += "<tr><td>-</td></tr>"

        sf += "</table>\n"
        s_star += "</table>"

        HTML.write('\n<td>\n{}\n</td>'.format(sf))
        HTML.write('''

            <td nowrap="nowrap">{}</td>
            {}
            <td><a href={}{}&JOB_TITLE={}{} target="_blank">blast</a></td>
            <td>{}</td>\n<td>{}</td>\n<td>{}</td>
            </tr>
        '''.format(
            hash_q[_id]["remaining_rc"],
            blat,
            blast,
            hash_q[_id]['pri_seq'],
            hash_q[_id]["oid"],
            blast_query,
            s_mat,
            s_star,
            hash_q[_id]['pri_seq']
        ))


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
    global xc, yc
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


def CreatePDFQuantifier(_hash, filename):
    global rna_d, xc, yc, bpo2r, bpo1r, bpo1, bpo2, y, sid, hairpin2mature
    fname = '{}/pdfs_{}/{}.pdf'.format(cwd, ltime, filename)
    c = canvas.Canvas(fname, pagesize=A4)
    spacer = len(sid)
    # c.setFont('Times-Roman', 20)
    font = 'Times-Roman'
    madd = 60
    drawLabel(c, xposshift + 20, y + 300 + downy,
              font, 8, 'miRBase precursor', black)
    drawLabel(c, xposshift + 110, y + 300 + downy, font, 8, ': ' + sid, black)

    spaces = ' ' * (spacer - len(str(_hash[sid]['freq_total'])))
    drawLabel(c, xposshift + 20, y + 230 + madd +
              downy, font, 8, 'Total read count', black)
    drawLabel(c, xposshift + 110, y + 230 + madd + downy, font,
              8, ": " + str(_hash[sid]['freq_total']), black)
    jk = 10

    h2m = re.split(r',', hairpin2mature[sid])
    for h in h2m:
        if re.match(r'^\s*$', h):
            continue

        if options.get('-t'):
            if not re.search(options.get('-m'), h):
                continue

        spaces = ' ' * (spacer - len(h))
        drawLabel(c, xposshift + 20, y + 230 - jk + madd + downy,
                  font, 8, "{} read count".format(h), black)
        drawLabel(c, xposshift + 110, y + 230 - jk + madd + downy,
                  font, 8, ": " + _hash[sid]['mapped'][h], black)
        jk += 10

    spaces = ' ' * (spacer - len('remaining reads'))
    drawLabel(c, xposshift + 20, y + 230 - jk + madd +
              downy, font, 8, "remaining reads", black)
    drawLabel(c, xposshift + 110, y + 230 - jk + madd + downy,
              font, 8, ": " + _hash[sid]['remaining_rc'], black)
    jk += 10
    # c.setFont('Courier')

    # DrawStructure
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

        if re.search(r'/sequence', l):
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
    red_dist = abs(xc[mb + cshift] - xc[mb - 1 + cshift])

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
    scfactor = rx / (maxx - minx)

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
        Base(c, xc[i] + x, yc[i] + y, rna_d[i], assign_str[i - 1 + offset], 8)

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
        dx1 = (dx - scfactorl * dx) / 2
        dy1 = (dy - scfactorl * dy) / 2

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

    return c


def drawLabel(drawCanvas, x, y, fontName, fontSize, string, color):
    '''
    Draw a label
    '''
    c = drawCanvas
    c.setFont(fontName, fontSize)
    c.setFillColor(color)
    c.drawString(x, y, string)


def CreateHistogramQuantifier(drawCanvas):
    global hstruct, y, xposshift, lstruct_multi, totalreads
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

    p = c.beginPath()
    p.moveTo(xposshift + 20, y + 50)
    lastx = position_hash[0]
    lasty = (hstruct[0] / float(totalreads)) * 100 + y + 50
    c.drawPath(p)

    for i in range(0, len(sstruct)):
        c.setStrokeColor(assign_str[i])
        p = c.beginPath()
        p.moveTo(lastx, lasty)
        p.lineTo(position_hash[i + 1], (hstruct[i] /
                                        float(totalreads) * 100) + y + 50)
        c.drawPath(p)
        lastx = position_hash[i + 1]
        lasty = (hstruct[i] / float(totalreads)) * 100 + y + 50


def CreateAlignmentQuantifier(drawCanvas, _hash):
    global y, mat_pre_arf, rna, hash2sample, hash2order
    trb = 'Courier-Bold'
    c = drawCanvas
    y += 20
    for k1 in mat_pre_arf[sid].keys():
        if re.search(r'\*', k1):
            drawLabel(
                c, xposshift + (20 + ((mat_pre_arf[sid][k1]['start']) * multiplier)), y, trb, 6, k1, col_star_obs)
        else:
            drawLabel(
                c, xposshift + (20 + ((mat_pre_arf[sid][k1]['start']) * multiplier)), y, trb, 6, k1, col_mature)
        y -= 10

    for i in range(0, len(rna)):
        drawLabel(c, position_hash[i], y, trb, 6, rna[i], assign_str[i])

    drawLabel(c, xposshift + 25 + lstruct_multi, y, trb, 6, '-3\'', black)
    drawLabel(c, xposshift + 10, y, trb, 6, '5\'-', black)

    try:
        _obs = _hash[sid]['obs']
    except KeyError:
        _obs = None

    if _obs:
        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'obs', black)
    else:
        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'exp', black)

    if _obs:
        y -= 10
        for i in range(0, len(rna)):
            drawLabel(c, position_hash[i], y, trb,
                      6, rna[i], assign_str_exp[i])
        drawLabel(c, xposshift + 50 + lstruct_multi, y, trb, 6, 'exp', black)

    y -= 10

    # structx = ssplit(struct)
    # sadd = 0

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

                    y += 20
                    for k1 in mat_pre_arf[sid].keys():
                        if re.search(r'\*', k1):
                            drawLabel(
                                c, xposshift + (20 + ((mat_pre_arf[sid][k1]['start']) * multiplier)), y, trb, 6, k1, col_star_obs)
                        else:
                            drawLabel(
                                c, xposshift + (20 + ((mat_pre_arf[sid][k1]['start']) * multiplier)), y, trb, 6, k1, col_mature)
                        y -= 10
                    y -= 20
            y -= 10
    else:
        for k in hash_sort_key(hash2order, lambda x: x[1]):
            tag = hash2sample[k]
            drawLabel(position_hash[0], y, trb, 6, hash2seq[k], black)

            # matches and read numbers
            drawLabel(c, xposshift + 30 + lstruct_multi, y, trb,
                      6, int(hash2c[tag][hash2key[k]]), black)
            drawLabel(c, xposshift + 70 + lstruct_multi,
                      y, trb, 6, hash2mm[k], black)
            drawLabel(c, xposshift + 110 + lstruct_multi,
                      y, trb, 6, hash2sample[k], black)

            y -= 10
            if y < 100:
                # add a new page
                c.showPage()
                y = 800
                y += 20
                for k1 in mat_pre_arf[sid].keys():
                    if re.search(r'\*', k1):
                        drawLabel(
                            c, xposshift + (20 + ((mat_pre_arf[sid][k1]['start']) * multiplier)), y, trb, 6, k1, col_star_obs)
                    else:
                        drawLabel(
                            c, xposshift + (20 + ((mat_pre_arf[sid][k1]['start']) * multiplier)), y, trb, 6, k1, col_mature)
                    y -= 10

                for i in range(0, len(rna)):
                    drawLabel(c, position_hash[i], y,
                              trb, 6, rna[i], assign_str[i])

                y -= 10


def ClosePDF(drawCanvas):
    '''
    Save the canvas to file
    '''
    drawCanvas.save()


def CreateStructurePDFQuantifier(_hash):
    global hash2order, rna, hstruct, sstruct, y, yorig, sid, assign_str, assign_str_exp, lstruct_multi, totalreads, mb, lb, hash2c, hash2key, hash2order, hash2mm, hash2seq, hash2sample

    filename = None
    print_stderr('Creating PDF files\n')
    for k in hash_sort_key(_hash, lambda x: _hash[x[0]]['freq_total'] * -1):
        if not _hash[k]['pdf']:
            continue
        if not _hash[k]['freq_total']:
            continue

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
        try:
            sb = mat_pre_arf[sid]['sb']     # starting
        except KeyError:
            sb = 0
        lmature = 0  # string of mature
        try:
            mb = mat_pre_arf[sid]['mb']     # mature begin
        except KeyError:
            mb = 0
        sstruct = 0     # structure string
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

        desc2 = ssplit(_hash[sid]['exp'])
        for i in range(0, len(desc2)):
            if (desc2[i] == "f"):   # # assign_str now starts at 0 not at one
                assign_str_exp[i] = "black"
                assign_str[i] = "black"
            elif (desc2[i] == "l"):
                assign_str_exp[i] = col_loop
                assign_str[i] = col_loop
            elif (desc2[i] == "S"):
                assign_str_exp[i] = col_star_obs
                assign_str[i] = col_star_obs
            else:
                assign_str_exp[i] = col_mature
                assign_str[i] = col_mature

        pdir = '{}/pdfs_{}'.format(cwd, ltime)
        if not os.path.isdir(pdir):
            os.mkdir(pdir)

        _fname = '{}/pdfs_{}/{}.tmp'.format(cwd, ltime, filename)
        FOLD = open_or_die2(_fname, 'w+')
        sstruct = _hash[sid]["pri_struct"]
        lstruct_multi = ((len(sstruct) + 2) * multiplier)

        FOLD.write("5{}3\n".format(pri_seq))
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
                        if options.get('-W'):
                            dc /= weighted[_hash[sid]
                                           ["reads"][tag][read]["rid"]]

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

        drawCanvas = CreatePDFQuantifier(_hash, filename)

        # DrawStructure(filename) is combine into CreatePDFQuantifier

        if totalreads != '0':
            CreateHistogramQuantifier(drawCanvas)

        CreateAlignmentQuantifier(drawCanvas, _hash)

        y -= 20

        ClosePDF(drawCanvas)

        # ATTR filename_ss.ps does not exist
        # os.unlink("{}/pdfs_{}/{}_ss.ps".format(cwd, ltime, filename))
        try:
            os.unlink("{}/pdfs_{}/{}.tmp".format(cwd, ltime, filename))
            os.unlink("{}/pdfs_{}/rna.ps".format(cwd, ltime))
        except:
            pass
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
        sys.argv[1:], "ugv:f:ck:os:t:w:er:q:dx:y:ab:i:j:lm:M:PW:")
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

    if options.get('-W') and file_s(options.get('-W')):
        IN = open_or_die2(options.get('-W'), 'rb')
        while True:
            line = IN.readline()
            if not line:
                break

            m = re.findall(r'(\S+)\s+(\d+)', line)
            if m:
                m = m[0]
                weighted[m[0]] = m[1]

        IN.close()

    # everything else given to it corresponds to the samples
    files_mirnaex = re.split(",", options.get('-M'))
    for l in files_mirnaex:
        print_stderr(l, " file with miRNA expression values\n")

    if not options.get('-y'):
        die("no timestamp given with parameter y\n")

    ltime = options.get('-y')

    if options.get('-x') and not options.get('-q'):
        die("\nError:\n\toption -x can only be used together with option -q\n\n")

    if options.get('-w'):
        pdf_path = "http://localhost:8001/links/moR/{}".format(
            options.get('-w'))

    cwd = current_dir()
    if not options.get('-w'):
        pdf_path = 'file://{}'.format(cwd)

    # organism parameter
    org = organisms.get(options.get('-t'), '')

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

    PrintQuantifier()

    CloseHTML()

    os.system("cp expression_analyses/expression_analyses_{}/expression_{}.html expression_{}.html".format(
        ltime,
        ltime,
        ltime
    ))

    if options.get('-d') != '':
        mirbase = 1
        # Create pdf
        CreateStructurePDFQuantifier(hash_q)

    sys.exit(0)
