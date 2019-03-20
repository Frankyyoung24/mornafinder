#!/usr/bin/env python
from __future__ import print_function

import os
import re
import sys
import time
import urllib2

ltime = None
home = os.environ['HOME']


def print_stderr(*messages):
    print(*messages, file=sys.stderr, end='')


def die(*messages):
    print_stderr(*messages)
    sys.exit(-1)


def open_or_die2(f, mode):
    try:
        f = open(f, mode)
        return f
    except IOError:
        if mode.find('r') > 0:
            die('can not open file : %@\n'.format(f))
        elif mode.find('w'):
            die('can not create file : %@\n'.format(f))
        elif mode.find('a'):
            die('can not append to file : %@\n'.format(f))
        else:
            die('can not access file : %@\n'.format(f))
    return None


def rem_mor(f):
    filename = "{}/{}".format(home, f)
    if os.path.isfile(filename):
        os.system('cp ~/{} ~/{}_{}'.format(f, f, ltime))
        print_stderr("~/{} backup is ~/{}_{}\n".format(f, f, ltime))
        os.system('mv ~/{} ~/{}.bak'.format(f, f))
        OUT = open_or_die2(filename, 'wb')
        IN = open_or_die2("{}/{}.bak".format(home, f), 'rb')
        line = []
        tmp = None
        while True:
            l = IN.readline()
            if not l:
                break

            tmp = ""
            if re.search(r'moRNAFinder', l):
                line = re.split(':', l)
                for l in line:
                    if not re.search('moRNAFinder', l):
                        tmp += l
            else:
                OUT.write(l)
    else:
        print("file {}/{}.bak not found\n".format(home, f))


def checkBIN(a, b):
    os.system("{}> tmp 2>tmp2".format(a))

    IN = open_or_die2('tmp', 'rb')
    found = 1
    while True:
        l = IN.readline()
        if not l:
            break

        if re.search(b, l):
            found = 0

    IN.close()
    if found:
        IN = open_or_die2('tmp2', 'rb')
        while True:
            l = IN.readline()
            if not l:
                break

            if re.search(b, l):
                found = 0

    IN.close()
    return found


# this routine checks if files exists on remote servers
def check(url):
    try:
        r = urllib2.urlopen(url)
        if r.code == 200:
            return 1
        else:
            return 0
    except urllib2.HTTPError:
        return 0


if __name__ == "__main__":

    cwd = os.path.realpath(os.path.curdir)

    print_stderr('''\n
##################################################################
#
# This is the moRNA Finder installer.
# It will work under a bash,csh and ksh shell.
# It will try to download all necessary files and install them.
# Please restart your shell to make changes take effect
#
##################################################################

''')

    ltime = int(time.time())

    if len(sys.argv) > 1:
        new = sys.argv[1]
    else:
        new = ''

    if not re.search('no', new):
        print_stderr(
            'making backup of .bashrc,.bash_profile and .cshrc and removing all entries of moRNA Finder in those files\n')
        rem_mor(".bashrc")
        rem_mor(".bash_profile")
        rem_mor(".cshrc")

    gcc = os.popen('gcc --version 2 > &1').read()
    if re.search(r'(GCC)', gcc, re.IGNORECASE):
        die('\nError:\n\tno gcc compiler installed. Please install a gcc compiler\n')
    else:
        m = re.match(r'^gcc\s*\S*\s*(\d+\S+)\s*', gcc)
        if m:
            print_stderr('gcc version: {}'.format(m.groups()[0]))

    wget = os.popen('wget').read()
    curl = os.popen('curl 2>&1').read()

    dtool = ''
    dopt = ''

    if re.search(r'URL', wget, re.IGNORECASE):
        dtool = "wget"
    elif curl:
        dtool = "curl"
        dopt = " -O"
    else:
        die("No commandline download tool found on your system. Please install wget or curl on your machine\n")

    cdir = os.popen('pwd 2>&1').read().strip()
    err = None
    dfile = ''

    # only attach to config file if not yet existing
    _in = os.popen('grep  "{}/moRNAFinder:*" ~/.bashrc'.format(cdir)).read()
    if not _in:
        os.system("echo 'export PATH=$PATH:{}' >> ~/.bashrc".format(cdir))

    _in = os.popen('grep "{}/:*" ~/.bash_profile'.format(cdir)).read()
    if not _in:
        os.system("echo 'export PATH=$PATH:{}' >> ~/.bash_profile".format(cdir))

    in2 = None
    if os.path.isfile("~/.cshrc"):
        in2 = os.popen('grep  "{}:*" ~/.cshrc'.format(cdir)).read()

    if not os.path.isdir("essentials"):
        os.system('mkdir essentials')

    os.chdir("essentials")

    a = os.popen('uname -a').read()
    bowtie = None
    bowtie_version = "1.2.0"
    bowtie_version_update = "1.2"
    ret = checkBIN("bowtie", "Usage")
    if ret == 0:
        print_stderr("bowtie already installed, nothing to do ...\n")
    else:
        if not os.path.isdir("bowtie-{}".format(bowtie_version)):
            print_stderr(
                "Downloading bowtie {} binaries\n\n".format(bowtie_version))
            if re.search('Darwin', a, re.IGNORECASE):   # download Mac version
                bowtie = "bowtie-{}-macos-x86_64.zip".format(
                    bowtie_version_update)
            elif re.search('x86_64', a, re.IGNORECASE):
                bowtie = "bowtie-{}-linux-x86_64.zip".format(
                    bowtie_version_update)
            else:
                bowtie = "bowtie-{}-source.zip".format(bowtie_version_update)

            if not os.path.isfile(bowtie):
                if check("http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/{}/{}".format(bowtie_version, bowtie)):
                    err = os.system("{} http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/{}/{} {}".format(
                        dtool, bowtie_version, bowtie, dopt))

                    if err:
                        die("\nError:\n\t{} could not be downloaded\n\n\n".format(bowtie))

                elif check("http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/old/{}/{}".format(bowtie_version, bowtie)):
                    err = os.system("{} http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/old/{}/{} {}".format(
                        dtool, bowtie_version, bowtie, dopt))
                    if err:
                        die("\nError:\n\t{} could not be downloaded\n\n\n".format(bowtie))
                else:
                    die("\nError:\n\t{} not found on server http://netcologne.dl.sourceforge.net/project/bowtie-bio/bowtie/ \n\n\n".format(bowtie))

            if not os.path.isfile(bowtie):
                die("{} download failed \nPlease try to download bowtie manually from here http://bowtie-bio.sourceforge.net/index.shtml\n\n\n".format(bowtie))

            print_stderr("Installing bowtie binaries\n\n")
            err = os.system("unzip {} > /dev/null  2>&1".format(bowtie))

            if err:
                die("unzip {} was not successful. Please remove essentials dir and try again.\n".format(
                    bowtie))

        _in = os.popen(
            'grep  "{}/essentials/bowtie-{}:*" ~/.bashrc'.format(cdir, bowtie_version)).read()
        if not _in:
            os.system(
                "echo 'export PATH=$PATH:{}/essentials/bowtie-{}' >> ~/.bashrc".format(cdir, bowtie_version))

        _in = os.popen(
            'grep  "{}/essentials/bowtie-{}:*" ~/.bash_profile'.format(cdir, bowtie_version)).read()
        if not _in:
            os.system(
                "echo 'export PATH=$PATH:{}/essentials/bowtie-{}' >> ~/.bash_profile".format(cdir, bowtie_version))

        _in2 = os.popen(
            'grep  "{}/essentials/bowtie-{}:*" ~/.cshrc'.format(cdir, bowtie_version)).read()

    ret = checkBIN("RNAfold -h 2", "usage")
    if ret == 0:
        print_stderr("RNAfold already installed, nothing to do ...\n")
    else:
        if not os.path.isdir("ViennaRNA-1.8.4"):
            dfile = "ViennaRNA-1.8.4.tar.gz"
            if not os.path.isfile(dfile):
                print_stderr("Downloading Vienna package now\n\n")
                if check("http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-1.8.4.tar.gz"):
                    err = os.system(
                        "{} http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-1.8.4.tar.gz {}".format(dtool, dopt))
                    if err:
                        die("Download of Vienna package not successful\n\n")
                else:
                    die("Vienna package not found at  http://www.tbi.univie.ac.at/RNA/packages/source/ViennaRNA-1.8.4.tar.gz\nPlease try to download the Vienna package from here http://www.tbi.univie.ac.at/RNA/RNAfold.html\n")

            if not os.path.isfile("ViennaRNA-1.8.4.tar.gz"):
                die("Vienna package download failed\n")

            print_stderr("Installing Vienna package now \n\n")
            os.system('tar xvvzf ViennaRNA-1.8.4.tar.gz > /dev/null 2>&1')
            os.chdir("ViennaRNA-1.8.4")
            os.system(
                './configure --prefix={}/essentials/ViennaRNA-1.8.4/install_dir > /dev/null 2>&1'.format(cdir))
            os.system('make > /dev/null 2>&1')
            os.system('make install > /dev/null 2>&1')

            os.chdir("..")

    _in = os.popen(
        'grep "{}/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.bashrc'.format(cdir)).read()
    if not _in:
        print_stderr("Vienna package path has been added to $PATH variable\n")
        os.system(
            "echo 'export PATH=$PATH:{}/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.bashrc".format(cdir))

    _in = os.popen(
        'grep "{}/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.bash_profile'.format(cdir)).read()
    if not _in:
        print_stderr("Vienna package path has been added to $PATH variable\n")
        os.system(
            "echo 'export PATH=$PATH:{}/essentials/ViennaRNA-1.8.4/install_dir/bin' >> ~/.bash_profile".format(cdir))

    in2 = os.popen(
        'grep  "{}/essentials/ViennaRNA-1.8.4/install_dir/bin:*" ~/.cshrc'.format(cdir))
    if not in2:
        pass

    ret = checkBIN("randfold", "let7")
    if ret == 0:
        print_stderr("\nrandfold already installed, nothing to do ...\n")
    else:
        dfile = "squid.tar.gz"
        if not os.path.isfile(dfile):
            print_stderr("Downloading SQUID library now\n")
            os.system(
                "{} http://eddylab.org/software/squid/squid.tar.gz {}".format(dtool, dopt))

        if not os.path.isfile("squid.tar.gz"):
            die("squid could not be downloaded\n Please try to download the library from here http://selab.janelia.org/software.html")

        if not os.path.isdir("squid"):
            print_stderr("Extracting squid and configuring it now\n")
            os.system('tar xxvzf squid.tar.gz > /dev/null 2>&1')
            os.chdir("{}/essentials/squid-1.9g".format(cdir))
            os.system('./configure > /dev/null 2>&1')
            os.system('make > /dev/null 2>&1')
            os.chdir("..")

        dfile = "randfold-2.0.tar.gz"
        if not os.path.isfile(dfile):
            print_stderr("Downloading randfold now\n")
            os.system(
                '{} http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/downloads/randfold-2.0.tar.gz {}'.format(dtool, dopt))

        if not os.path.isfile("randfold-2.0.tar.gz"):
            die("randfold could not be downloaded\nPlease try to download randfold from here http://bioinformatics.psb.ugent.be/software/details/Randfold\n")

        if not os.path.isdir("randfold-2.0"):
            print_stderr("\nInstalling randfold now\n")
            os.system('tar xvvzf randfold-2.0.tar.gz > /dev/null 2>&1')
            os.chdir("{}/essentials/randfold-2.0".format(cdir))
            IN = open_or_die2("Makefile", 'rb')
            OUT = open_or_die2("Makefile_new", 'wb')

            while True:
                l = IN.readline()
                if not l:
                    break

                if re.search(r'INCLUDE=-I\.\s*', l, re.IGNORECASE):
                    OUT.write(
                        "INCLUDE=-I. -I{}/essentials/squid -L{}/essentials/squid/\n".format(cdir, cdir))
                else:
                    OUT.write('\n')

            IN.close()
            OUT.close()

            os.system('mv Makefile Makefile.orig')
            os.system('mv Makefile_new Makefile')
            os.system('make > /dev/null 2>&1')
            os.chdir("..")

        _in = os.popen(
            'grep  "{}/essentials/randfold-2.0:*" ~/.bashrc'.format(cdir))
        if not _in:
            print_stderr("Randfold path has been added to $PATH variable\n")
            os.system(
                "echo 'export PATH=$PATH:{}/essentials/randfold-2.0' >> ~/.bashrc".format(cdir))

        _in = os.popen(
            'grep  "{}/essentials/randfold-2.0:*" ~/.bash_profile'.format(cdir))
        if not _in:
            print_stderr("Randfold path has been added to $PATH variable\n")
            os.system(
                "echo 'export PATH=$PATH:{}/essentials/randfold-2.0' >> ~/.bash_profile".format(cdir))

        # $in2 = `grep  "$dir/essentials/randfold-2.0:*" ~/.cshrc`;
        # if($in2){
        #   `echo 'setenv PATH $PATH:$dir/essentials/randfold-2.0' >> ~/.cshrc`;

    # install pip & python requirement
    _in = os.popen('grep  "{}/.local/bin" ~/.bashrc'.format(home)).read()
    if not _in:
        os.system(
            "echo 'export PATH=$PATH:{}/.local/bin' >> ~/.bashrc".format(home))

    _in = os.popen('grep  "{}/.local/bin" ~/.bash_profile'.format(home)).read()
    if not _in:
        os.system(
            "echo 'export PATH=$PATH:{}/.local/bin' >> ~/.bash_profile".format(home))

    ret = checkBIN("pip", "Usage")
    if ret == 0:
        print_stderr("\npip have been installed.\n")
    else:
        os.chdir(cwd)
        print_stderr('\nInstalling python pip & requirement packages...\n')
        r = os.system('python get-pip.py --user')
        if r:
            die('Intall python pip failed. Please install the pip tool manually. Download address: https://bootstrap.pypa.io/get-pip.py\n')
        else:
            _in = os.popen(
                'grep  "{}/.local/bin" ~/.bashrc'.format(home)).read()
            if not _in:
                os.system(
                    "echo 'export PATH=$PATH:{}/.local/bin' >> ~/.bashrc".format(home))

            _in = os.popen(
                'grep  "{}/.local/bin" ~/.bash_profile'.format(home)).read()
            if not _in:
                os.system(
                    "echo 'export PATH=$PATH:{}/.local/bin' >> ~/.bash_profile".format(home))

    os.chdir(cwd)
    if not os.path.isfile('./requirements.txt'):
        die('requirements.txt file is not exist.\n')

    print_stderr('Installing python requirement packages...\n')
    r = os.system(
        '{}/.local/bin/pip install  --user -r requirements.txt > /dev/null 2>&1'.format(home))
    if r:
        die('moRNA Finder related package install failed. Please install the packages manually with command:\n\tpip install --user -r requirements.txt \n')

    print_stderr("Installation successful\nPlease start a new shell\n")
