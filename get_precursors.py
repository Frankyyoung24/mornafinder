#!/usr/bin/env python

import getopt
import os
import re
import sys

from port import die, open_or_die2, tr

usage = '''
Usage: {} -r result.csv  [options]

[options]
    -s [int]\toutput only precursors with min-score >= [int]
    -t [int]\toutput only precursors with score     <  [int]
    -d      \toutput dna instead of rna
    -p      \tmake simple id with the name only
    -h      \tprint this help message
    -m      \tget_mature instead of precursor
    -k      \tget_star instead of precursor
    -T      \tTrackname for bedfiles
    -o      \toutdir

'''.format(sys.argv[0])

if __name__ == '__main__':
    opt, args = getopt.getopt(sys.argv[1:], "r:s:dpht:mkT:o:")
    options = dict(opt)

    if not options.get('-r') or not os.path.isfile(options.get('-r')) or options.get('-h') == '':
        die(usage)

    od = "result_mirnas"
    if options.get('-o'):
        od = options.get('-o')

    if not os.path.isdir(od):
        os.mkdir(od)

    if not re.search(r'\S+', options.get('-T', '')):
        options['-T'] = 'notTrackname'

    l1 = re.split('result', options.get('-r'))
    timestamp, dmp = re.split(".csv", l1[1])

    known, novel, _not, line, thres, score, strandcol, bedh1, pcoord = (
        None, None, None, None, None, None, None, None, '', )

    bedh = "browser position BPOS\nbrowser hide all\ntrack name=\"TNAME\" description=\"TDESC\" visibility=2\nitemRgb=\"On\";\n"
    first = 1

    line = []

    thres = -50

    if options.get('-s') is not None:
        thres = options.get('-s')

    score = thres
    _max = 'na'
    maxs = 999999999999999999999999999
    if options.get('-t'):
        _max = options.get('-t')
        maxs = _max

    IN = open_or_die2(options.get('-r'), 'rb')
    seqcol = 15
    if options.get('-m') == '':
        seqcol = 13

    if options.get('-k') == '':
        seqcol = 14

    names = ('0', '1', '2', '3', '4', '5', '6', '7', '8',
             '9', '10', '11', '12', 'mature', 'star', 'pres')

    while True:
        l = IN.readline()
        if not l:
            break

        if re.search('novel miRNAs predicted by moRNA Finder', l):
            novel = 1
            known = 0
            _not = 0
            line = IN.readline()
            first = 1
            OUT = open_or_die2("{}/novel_{}{}_score{}_to_{}.fa".format(od,
                                                                       names[seqcol], timestamp, score, _max), 'wb')
            BED = open_or_die2("{}/novel_{}{}_score{}_to_{}.bed".format(od,
                                                                        names[seqcol], timestamp, score, _max), 'wb')
            continue
        elif re.search('miRBase miRNAs detected', l):
            OUT.close()
            BED.close()
            novel = 0
            known = 1
            _not = 0
            line = IN.readline()
            first = 1
            OUT = open_or_die2("{}/known_{}{}_score{}_to_{}.fa".format(od,
                                                                       names[seqcol], timestamp, score, _max), 'wb')
            BED = open_or_die2("{}/known_{}{}_score{}_to_{}.bed".format(od,
                                                                        names[seqcol], timestamp, score, _max), 'wb')
            continue
        elif re.search('miRBase miRNAs not detected by moRNA Finder', l):
            novel = 0
            known = 0
            _not = 1
            first = 1
            OUT.close()
            BED.close()
            OUT = open_or_die2("{}/not_{}{}_score{}_to_{}.fa".format(od,
                                                                     names[seqcol], timestamp, score, _max), 'wb')
            BED = open_or_die2("{}/not_{}{}_score{}_to_{}.bed".format(od,
                                                                      names[seqcol], timestamp, score, _max), 'wb')
            continue

        if re.match(r'^\s*$', l):
            continue

        if novel or known:
            l = l.strip()
            line = re.split('\t', l)
            coord = 'na'
            if len(line) > 16 and line[16]:
                coord = line[16]

            if known:
                if float(line[1]) >= thres and float(line[1]) < maxs:
                    if options.get('-d') == '':
                        line[seqcol] = tr(line[seqcol], 'uU', 'tT')

                    if options.get('-p') == '':
                        m = re.search(r'\|([a-zA-Z0-9_-]*)$', line[0])
                        if m:
                            line[0] = m.groups()[0]

                        OUT.write(">{}\n{}\n".format(
                            line[0], line[seqcol].upper()))

                    else:
                        OUT.write(">{}_{}_x{}_coord:{}_score:{}\n{}\n".format(
                            line[0],
                            line[9],
                            line[5],
                            coord,
                            line[1],
                            line[seqcol].upper()
                        ))

                    coord = coord.strip()

                    strandcol = "0,0,255"

                    m = re.match('^(\S+):(\d+)\.\.(\d+):(\S)$', coord)
                    if m:
                        m = m.groups()
                        pcoord = "{}:{}-{}".format(m[0], m[1], m[2])
                        if m[3] == "+":
                            strandcol = "255,0,0"

                        if first:
                            first = 0
                            bedh1 = bedh
                            bedh1 = re.sub('TNAME', '{}.known_miRNAs'.format(
                                options['-T']), bedh1, count=1)
                            bedh1 = re.sub('TDESC', 'known miRNAs detected by moRNA Finder for {}'.format(
                                options['-T']), bedh1, count=1)
                            bedh1 = re.sub('BPOS', pcoord, bedh1, count=1)
                            BED.write(bedh1)

                        BED.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            m[0],
                            m[1],
                            m[2],
                            line[0],
                            line[1],
                            m[3],
                            m[1],
                            m[2],
                            strandcol
                        ))
            else:
                if float(line[1]) >= thres and float(line[1]) < maxs:
                    if options.get('-d') == '':
                        line[seqcol] = tr(line[seqcol], 'uU', 'tT')

                    if options.get('-p') == '':
                        m = re.search(r'\|([a-zA-Z0-9_-]*)$', line[0])
                        if m:
                            line[0] = m.groups()[0]

                        OUT.write(">{}\n{}\n".format(
                            line[0], line[seqcol].upper()))

                    else:
                        OUT.write(">{}_x{}_coord:{}_score:{}\n{}\n".format(
                            line[0],
                            line[5],
                            coord,
                            line[1],
                            line[seqcol].upper()
                        ))

                    coord = coord.strip()
                    strandcol = "0,0,255"
                    m = re.match(r'^(\S+):(\d+)\.\.(\d+):(\S)$', coord)
                    if m:
                        m = m.groups()
                        if m[3] == "+":
                            strandcol = "255,0,0"

                        if first:
                            first = 0
                            bedh1 = bedh
                            bedh1 = re.sub('TNAME', '{}.novel_miRNAs'.format(
                                options.get('-T')), bedh1, count=1)
                            bedh1 = re.sub('TDESC', 'novel miRNAs detected by moRNA Finder for {}'.format(
                                options.get('-T')), bedh1, count=1)
                            bedh1 = re.sub('BPOS', pcoord, bedh1, count=1)
                            BED.write(bedh1)

                        BED.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            m[0],
                            m[1],
                            m[2],
                            line[0],
                            line[1],
                            m[3],
                            m[1],
                            m[2],
                            strandcol
                        ))
                continue

        if _not:
            l = l.strip()
            line = re.split('\t', l)
            coord = 'na'
            if len(line) > 16 and line[16]:
                coord = line[16]

            if options.get('-d') == '':
                try:
                    line[seqcol] = tr(line[seqcol], 'uU', 'tT')
                except IndexError:
                    pass

            if options.get('-p') == '':
                pass
            else:
                OUT.write(">{}_x{}\n{}\n".format(
                    line[0], line[4], line[seqcol]))

            coord = coord.strip()
            strandcol = "0,0,255"
            m = re.match(r'^(\S+):(\d+)\.\.(\d+):(\S)$', coord)
            if m:
                if m[3] == "+":
                    strandcol = "255,0,0"

                if first:
                    first = 0
                    bedh1 = bedh
                    bedh1 = re.sub('TNAME', '{}.not_detected_miRNAs'.format(
                        options.get('-T')), bedh1, count=1)
                    bedh1 = re.sub('TDESC', 'miRNAs not detected by moRNA Finder for {}'.format(
                        options.get('-T'), bedh1, count=1))
                    bedh1 = re.sub('BPOS', pcoord, bedh1, count=1)
                    BED.write(bedh1)

                BED.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    m[0],
                    m[1],
                    m[2],
                    line[0],
                    line[1],
                    m[3],
                    m[1],
                    m[2],
                    strandcol
                ))

            continue

    OUT.close()
