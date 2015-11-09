#!/usr/bin/env python

import sys, os
import pre
from subprocess import *
from string import *

def program_prepair():
	# define program used in this script
	programs = pre.default_settings()
	global svmscale, svmtrain
	global subset_py, grid_py, is_win32
	is_win32 = (sys.platform == 'win32')

	svmscale = pre.svm_scale
	svmtrain = pre.svm_train
	subset_py = pre.subset_py
	grid_py = pre.grid_py

def get_para(argv = sys.argv):
        """Train svm model. version 3.0 update: 2012-06-11. It can be used under both windows and unix.
Usage: train-cds-model.py  [options]  featurefile

    -s        subset number (default, 5000)
    -m        cachesize (default 500 (M))
    -b        1/0, 1 for probability model (default, 1)
    -l        x scaling lower limit (default, -1)
    -u        s scaling upper limit (default, +1)
    -log2c    start,stop,step (default, 0,7,1)
    -log2g    start,stop,step (default, 0,-7,-1)
    """
	if len(argv) < 2:
                print get_para.__doc__
		sys.exit(1)

        global subset_number, cachesize, prob, lower, upper
        global log2c, log2g
        global feature
        
        subset_number = 5000
        cachesize = 500
        prob = 1
        lower = -1
        upper = 1
        log2c = '0,7,1'
        log2g = '0,-7,-1'
        
	argc = len(argv)
	i = 1
	while i < argc:
                if argv[i][0] != '-':
                        break
                if argv[i] == '-s':
                        i = i + 1
                        subset_number = int(argv[i])
                elif argv[i] == '-m':
                        i = i + 1
                        cachesize = int(argv[i])
                elif argv[i] == '-b':
                        i = i + 1
                        prob = int(argv[i])
                elif argv[i] == '-l':
                        i = i + 1
                        lower = int(argv[i])
                elif argv[i] == '-u':
                        i = i + 1
                        upper = int(argv[i])
                elif argv[i] == '-log2c':
                        i = i + 1
                        log2c = argv[i]
                elif argv[i] == '-log2g':
                        i = i + 1
                        log2g = argv[i]
                else:
                        i = i + 1
                i = i + 1

        if i >= argc:
                print get_para.__doc__
                exit(1)
        
	feature = argv[i]

	print feature, subset_number, cachesize

def scale(lower, upper):
	print '*\nScale featrue'
	print '\tlower = %i, upper = %i' % (lower, upper)
	global scaleF, rangeF
	scaleF = os.path.split(feature)[1] + '.scale'
	rangeF = os.path.split(feature)[1] + '.range'
	cmd = '%s -l %i -u %i -s %s "%s" > %s' % (svmscale, lower, upper, rangeF, feature, scaleF)
	call(cmd, shell = True)

def grid(log2c, log2g):
	print '*\nSubset featrue for grid search'
	global subF, crossF, outF
	subF = scaleF + '.sub'
	crossF = subF + '.cross-acc'
	outF = subF + '.out'
	cmd = '%s -s 1 %s %s %s' % (subset_py, scaleF, subset_number, subF)
	call(cmd, shell = True)
	
	print '*\nGrid search for the best c and gamma parameter'
	print '\tlog2c = %s, log2g = %s' % (log2c, log2g)
	global c, g, rate
	cmd = '%s -log2c %s -log2g %s -xtest %s %s' % (grid_py, log2c, log2g, subF, subF)
	f = os.popen(cmd, 'r')
	line = ''
	while True:
		last_line = line
		line = f.readline()
		if not line: break
	c, g, rate = map(float, last_line.split())
	print '*\nBest c=%s, g=%s subset Harmonic-mean=%s' % (c,g,rate)
	
def train():
	print '*\nTraining probability CDS model'
	modelF = os.path.split(feature)[1] + '.model'
	cmd = '%s -m %i -b %s -c %s -g %s "%s" "%s"' \
	% (svmtrain, cachesize, prob, c, g, scaleF, modelF)
	print cmd
	call(cmd, shell = True)
	print '*\nCDS Model was Trained'
	
def delF():
	os.remove(subF)
	os.remove(scaleF)
	
def main():	
	argv = sys.argv
	argc = len(argv)
	if argc < 2:
                print get_para.__doc__
                return
	
	# prepairing programs to be used
	program_prepair()

	# get parameters from command line	
	get_para()

	# model trainning
	scale(lower, upper)
	grid(log2c, log2g)
	train()
	
	# delete intermediate files
	delF()

if __name__ == "__main__":
        main()
