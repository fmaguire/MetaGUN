#!/usr/bin/env python

__author__ = "LIU YongChu"
__update__ = "2010-09-05"
__update__ = "2012-06-11"

import pre
import sys, os, traceback
from string import *
from subprocess import *
from OftenUsedFun import *

global is_win
if 'win' in sys.platform:
    is_win = 1
else:
    is_win = 0
    
def main():
    argc = len(sys.argv)
    argv = sys.argv
    if argc == 1:
        print "Usage:  domain-search.py  input-sequences  search-results-filename"
        return

    # load default settings
    pre.default_settings()
    pre.seqs_file = argv[1]
    pre.hit_file = argv[2]
    
    if len(pre.seqs_file) == 0:
        print 'The sequence file is not defined.'
        return
    if not os.path.exists(pre.seqs_file):
        print 'The sequence file %s is not correctly defined.' % pre.seqs_file
        return

    print_section_title('Novel predict')

    print '  Blast searches against known database:', pre.blast_db
    alc_file = pre.seqs_file + '.alc'
    csqa_file = pre.seqs_file + '.csqa'
    pre.project_name = pre.seqs_file
    pre.hit_file = argv[2]
    
    cmd = '%s -p -n %s -a -s %s' % (pre.get_data, pre.project_name, pre.seqs_file)
    call(cmd, shell = True, stdout = PIPE)
    cmd = '%s -n %s -s %s -c %s -q 1' % (pre.get_data, pre.project_name, pre.seqs_file, alc_file)
    call(cmd, shell = True, stdout = PIPE)
    cmd = '%s -i %s -d %s -m 8 -v 1 -b 1 -e 1 -o %s' % (pre.rpsblast, csqa_file, pre.blast_db, pre.hit_file)
    call(cmd, shell = True, stdout = PIPE)
    rm_files([alc_file, csqa_file])
    print '  Blast search results save in:', pre.hit_file
    return           

if __name__=="__main__":
    main()
