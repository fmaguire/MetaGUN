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
        print "Usage:  bin-seq.py  input-sequences  bin-results-filename"
        return

    # load default settings
    pre.default_settings()
    
    pre.seqs_file = argv[1]
    pre.bin_file = argv[2]
    bindistri_file = pre.bin_file + 'distri'
    
    if len(pre.seqs_file) == 0:
        print 'The sequence file is not defined.'
        return
    if not os.path.exists(pre.seqs_file):
        print 'The sequence file %s is not correctly defined.' % pre.seqs_file
        return
    
    print '  Input sequence file:', pre.seqs_file

    print_section_title('Sequence Binning')
    print '  Sequence binning model:', pre.binmodel
    print '  Taxonomy map file:', pre.taxonomy
    cmd = '%s -m %s -s %s -t %s -b %s -d %s' % \
          (pre.bin_seqs, pre.binmodel, pre.seqs_file, pre.taxonomy, pre.bin_file, bindistri_file)
    call(cmd, shell = True, stdout = PIPE)
    print '  Binning result file:', pre.bin_file
    return

if __name__=="__main__":
    main()
