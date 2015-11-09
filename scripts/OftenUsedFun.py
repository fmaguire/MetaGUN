#!/usr/bin/env python

__author__ = "LiuYC"
__version__= "1.0"
__create__ = "20120905"
__update__ = "20120905"

import os, sys
from string import *
from subprocess import *

#############################################
#  delete intermediate files  
#############################################
def rm_files(files):
    for f in files:
        if os.path.exists(f):
            os.remove(f)
###############################################################################################
#  check if file exists
###############################################################################################
def check_file_existence(files):
     for this_file in files:
          assert os.path.exists(this_file), "%s not found, please check it." % this_file

###############################################################################################
#  print title for section description  
###############################################################################################
def print_section_title(descript):
    width = 80
    padding = int((width - len(descript))/2)
    print '\n', ' '.rjust(padding, '-'), descript, ' '.ljust(width-padding-len(descript), '-')

###############################################################################################
#  parse command line
###############################################################################################
def command_line_parser(argv, short_opt, long_opt, has_argv):
    option = [False, '']
    argc = len(argv)
    for i in range(0, argc):
        if argv[i] == short_opt or argv[i] == long_opt:
            if has_argv:
                if i+1 >= argc:
                    print 'Option %s is not specified.' % argv[i]
                    option[0] = True
                    return option
                else:
                    option[0] = True
                    option[1] = argv[i+1]
                    return option
            else:
                option[0] = True
                return option
    return option

if __name__ == "__main__":
    option = command_line_parser(sys.argv, '-a', '--argc', True)
    print option
