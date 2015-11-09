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

############################################################
#  detailed helps
############################################################
def detail_help():
    """MetaGUN  Latest update: 2012/09/06.
metagun.py  [OPTIONS]
Options:
   -h, --help                     print this detailed help message and exit;
   -n, --name  STRING             project name, suffix will be added suffix after that, default, 'sample';
   -s, --seqs  FILE               input sequence file, multiple FASTA format;
   -c, --cds-model  PATH          directory path of cds-model, where supervised CDS models stored
   -t, --tis-model  PATH          directory path of tis-model, where supervised TIS models stored
   -m, --bin-model  FILE          binning model file
   -x, --tax-map  FILE            taxonomic map file for binning
   -b, --bin-result  FILE         file stores the output binning results
   -p, --svm-porb  INT            0/1, set 1 to switch on the probability svm model, default: 0
   -B, --blast-hits  FILE         rpsblast hit file, if given, the blast Pfam procedure will be skipped
       --svm-cv  FLOAT            cut probability value in universal prediction to determine if an ORF
                                  is coding or not, default, 0.5
       --svm-cv2 FLOAT            cut probability value in novel prediction to determine if an ORF
                                  is coding or not, default, 0.5
       --svm-ssize  INT           sub-sampling size of the training sets for training the SVM paramters,
                                  default 5000
       --blast-ev  FLOAT          cut e-value of the combined blast hits, default, 1e-5
       --seeds-ev  FLOAT          cut e-vaule for choosing seed blast hits, default, 1e-40
   -U, --run-uni-pred  INT        0/1, set 1 to run universal prediction module, default, 1
   -N, --run-novel-pred  INT      0/1, set 1 to run novel prediction model, default, 1
   -M, --run-metatisa  INT        0/1, set 1 to run MetaTISA for TIS refinement, default, 1.
   -u, --unix  INT                0/1, set 1 for the unix platform and 0 for the windows platform, default, 1.
"""

############################################################
#  parse command line
############################################################
def parse_command(argv):
    """MetaGUN  Latest update: 2012/09/06.
metagun.py  [-n project_name]  [-s input_sequence_file]  [-h]
Options:
    -n STRING    project name, suffix will be added suffix after that, default, 'sample';
    -s FILE      input sequence file, multiple FASTA format;
    -b FILE      specify the fragment classifiction result file;
    -B FILE      specify the domain search result file;

    Generally, run:
        metagun.py -n example -s example.seq
    If fragment classification results and domain search results are obatined separately, run:
        metagun.py -n example -s example.seq -b example.bin -B example.blastcdd
"""
    opt = command_line_parser(argv, '-n', '--name', True)
    if opt[0]:
        pre.project_name = opt[1]
    opt = command_line_parser(argv, '-s', '--seqs', True)
    if opt[0]:
        pre.seqs_file = opt[1]
    opt = command_line_parser(argv, '-c', '--cds-model', True)
    if opt[0]:
        pre.cdsmodel = opt[1]
    opt = command_line_parser(argv, '-t', '--tis-model', True)
    if opt[0]:
        pre.tismodel = opt[1]
    opt = command_line_parser(argv, '-m', '--bin-model', True)
    if opt[0]:
        pre.binmodel = opt[1]
    opt = command_line_parser(argv, '-x', '--tax-map', True)
    if opt[0]:
        pre.taxonomy = opt[1]
    opt = command_line_parser(argv, '-b', '--bin-result', True)
    if opt[0]:
        pre.exists_bin_file = True
        pre.bin_file = opt[1]
    opt = command_line_parser(argv, '-p', '--svm-prob', True)
    if opt[0]:
        pre.prob_status = int(opt[1])
    opt = command_line_parser(argv, '-B', '--blast-hits', True)
    if opt[0]:
        pre.exists_hit_file = True
        pre.hit_file = opt[1]
    opt = command_line_parser(argv, '', '--svm-cut-value', True)
    if opt[0]:
        pre.svm_cut_value = float(opt[1])
    opt = command_line_parser(argv, '', '--svm-cut-value2', True)
    if opt[0]:
        pre.svm_cut_value2 = float(opt[1])
    opt = command_line_parser(argv, '', '--svm-sub-size', True)
    if opt[0]:
        pre.svm_sub_size = int(opt[1])
    opt = command_line_parser(argv, '', '--blast-e-value', True)
    if opt[0]:
        pre.blast_ev = float(opt[1])
    opt = command_line_parser(argv, '', '--seeds-e-value', True)
    if opt[0]:
        pre.seeds_ev = float(opt[1])
    opt = command_line_parser(argv, '-U', '--run-uni-pred', True)
    if opt[0]:
        if int(opt[1]) == 0:
            pre.run_uni_pred = False
    opt = command_line_parser(argv, '-N', '--run-novel-pred', True)
    if opt[0]:
        if int(opt[1]) == 0:
            pre.run_novel_pred = False
    opt = command_line_parser(argv, '-M', '--run-metatisa', True)
    if opt[0]:
        if int(opt[1]) == 0:
            pre.run_metatisa = False

############################################################
#  predict coding sequences by SVM model
############################################################
def svm_predict(range_file, model_file, fea_file, orf_file, cut_value, prob_status, cds_pred_file):
    """SVM predict using previously trained model."""
    fea_scale_file = fea_file + '.scale'
    svm_pred_file = fea_file + '.predict'
    
    # scale
    cmd = '%s -r %s %s > %s' % (pre.svm_scale, range_file, fea_file, fea_scale_file)
    call(cmd, shell = True, stdout = PIPE)
    
    # predict
    cmd = '%s -b %i %s %s %s' % (pre.svm_predict, prob_status, fea_scale_file, model_file, svm_pred_file)
    call(cmd, shell = True, stdout = PIPE)
    
    # get predicted CDS
    cmd = '%s -b %i -o %s -c %f %s %s' % (pre.extr_pred_cds, prob_status, cds_pred_file, cut_value, orf_file, svm_pred_file)
    call(cmd, shell = True, stdout = PIPE)
    
    rm_files([fea_scale_file, svm_pred_file])


############################################################
#  predict coding sequences by universal SVM model
############################################################
def uni_pred(uni_pred_file):
    """Predict genes by the supervised universal CDS model."""
    print_section_title('Universal predict')
    
    # get all bins
    allbins = []
    fr = open(pre.bin_file, 'r')
    for line in fr.readlines():
        binname = line.split()[-1]
        if binname not in allbins:
            allbins.append(line.split()[-1])
    fr.close()
    allbins.sort()
     
    # start to predict
    is_unix = 0
    if is_win == 0:
        is_unix = 1
    print '  Get all candidate ORFs ...'
    cmd = '%s -p -n %s -s %s -t %s -b %s -u' % \
          (pre.get_data, pre.project_name, pre.seqs_file, pre.tismodel, pre.bin_file)
    if is_unix == 0:
        cmd = '%s -p -n %s -s %s -t %s -b %s' % \
              (pre.get_data, pre.project_name, pre.seqs_file, pre.tismodel, pre.bin_file)
    call(cmd, shell = True, stdout = PIPE)

    print '  Predict using universal model ...'
    print '  SVM probability mode:', pre.prob_status
    print '  SVM scroe cut-vaule:', pre.svm_cut_value
    cds_model_path = pre.cdsmodel
    name = pre.project_name
    cut_value = pre.svm_cut_value
    uni_w = open(uni_pred_file, 'w')
    uni_w.close()
    count = 1
    types = ['00', '01', '10', '11']
    for g in allbins:
        print '\t', g
        for t in types:
            fea_file = name + '.' + g + '.fea_' + t
            orf_file = name + '.' + g + '.orf_' + t
            if not os.path.exists(fea_file):
                continue
            model_file = cds_model_path + os.sep + g + '_' + t + '.model'
            range_file = cds_model_path + os.sep + g + '_' + t + '.range'
            cds_pred_file = fea_file + '.cds'
            svm_predict(range_file, model_file, fea_file, orf_file, cut_value, pre.prob_status, cds_pred_file)

            cmd = '%s -a %s -b %s --a+b -o %s' % (pre.metalocs_operate, cds_pred_file, uni_pred_file, uni_pred_file)
            call(cmd, shell = True, stdout = PIPE)
            rm_files([fea_file, orf_file, cds_pred_file])
        count = count + 1
    return

############################################################
#  count the number of lines for input file
############################################################
def line_counter(input_file):
    count = -1
    for count, line in enumerate(open(input_file, 'rU')):
        pass
    return count+1

def trainingset_counter(input_file):
    trainingset = {'+1': 0, '-1': 0}
    if not os.path.exists(input_file):
        return trainingset
    fin = open(input_file, 'r')
    for line in fin.readlines():
        trainingset[line[0:2]] = trainingset[line[0:2]] + 1
    return trainingset
    
############################################################
#  predict coding sequences by novel SVM model
############################################################
def novel_pred(novel_pred_file):
    """Predict genes with sample-specific svm model. If too few high similarity hits to known database, this procedure will be skipped."""
    print_section_title('Novel predict')

    # Get significantly conserved hits against known protein database
    if pre.exists_hit_file == False:
        print '  Blast searches against known database:', pre.blast_db
        alc_file = pre.project_name + '.alc'
        csqa_file = pre.project_name + '.csqa'
        pre.hit_file = csqa_file + '.blast'
          
        cmd = '%s -p -n %s -a -s %s' % (pre.get_data, pre.project_name, pre.seqs_file)
        call(cmd, shell = True, stdout = PIPE)
        cmd = '%s -n %s -s %s -c %s -q 1' % (pre.get_data, pre.project_name, pre.seqs_file, alc_file)
        call(cmd, shell = True, stdout = PIPE)
        cmd = '%s -i %s -d %s -m 8 -v 1 -b 1 -e 1 -o %s' % (pre.rpsblast, csqa_file, pre.blast_db, pre.hit_file)
        call(cmd, shell = True, stdout = PIPE)
        print cmd
        rm_files([alc_file, csqa_file])
    else:
        print '  Using given blast result file:', pre.hit_file

    # Get training sets for training SVM model based on significant hits
    is_unix = 0
    if is_win == 0:
        is_unix = 1
    print '  Seeds e-value cut:', pre.seeds_ev
    bh_file = pre.hit_file + '.bh'
    cmd = '%s -f %s -o %s -e %s' % (pre.parse_blast, pre.hit_file, bh_file, pre.seeds_ev)
    call(cmd, shell = True, stdout = PIPE)
    cmd = '%s -n %s -s %s -c %s -h -t %s -b %s -u' % \
          (pre.get_data, pre.project_name, pre.seqs_file, bh_file, pre.tismodel, pre.bin_file)
    if is_unix == 0:
        cmd = '%s -n %s -s %s -c %s -h -t %s -b %s' % \
          (pre.get_data, pre.project_name, pre.seqs_file, bh_file, pre.tismodel, pre.bin_file)
    call(cmd, shell=True, stdout=PIPE)

    #  Get predicting data
    cmd = '%s -p -n %s -s %s -t %s -b %s -u %i' % (pre.get_data, pre.project_name, pre.seqs_file, pre.tismodel, pre.bin_file, is_unix)
    call(cmd, shell = True, stdout = PIPE)

    #  Get binning groups
    fin_bin = open(pre.bin_file, 'r')
    bin_groups = []
    for line in fin_bin.readlines():
        group = line.split()[-1]
        if group not in bin_groups:
            bin_groups.append(group)

    threshold4training = 50
    print '  SVM subset size:', pre.svm_sub_size
    print '  Threshold of positive instance size for training:', threshold4training
    nn_1_file = novel_pred_file+'.1'
    fout_nn_1 = open(nn_1_file, 'w')
    for group in bin_groups:
        print '\t', group 
        enough4training = {'00' : False, '01' : False, '10' : False, '11' : False}
        orf_types = ['00', '01', '10', '11']
        # check if the positive instance is enough for training svm model
        for tt in orf_types:
            if trainingset_counter(pre.project_name + '.' + group + '.trainfea_' + tt)['+1'] > threshold4training:
                enough4training[tt] = True
#        print ' ', group, enough4training
        
        # train svm model based on genes with significant hits and use them for predicting
        this_nn_1_file = nn_1_file[0:nn_1_file.rfind('.')] + '.' + group + '.novel-pred'
        fout_this_novel = open(this_nn_1_file, 'w')
        fout_this_novel.close()
        
        for tt in orf_types:
            trainfea_file = pre.project_name + '.' + group + '.trainfea_' + tt
            fea_file = pre.project_name + '.' + group + '.fea_' + tt
            orf_file = pre.project_name + '.' + group + '.orf_' + tt
            range_file = trainfea_file + '.range'
            model_file = trainfea_file + '.model'          
            cds_pred_file = fea_file + '.cds'
            
            if enough4training[tt] == False:
                rm_files([trainfea_file, fea_file, orf_file])
                continue
            
            # train
            cmd = '%s -m 1000 -s %i -b %i -l 0 -u 1 %s' % (pre.train_cds_model_py, pre.svm_sub_size, pre.prob_status, trainfea_file)
            call(cmd, shell=True, stdout=PIPE)
            rm_files([trainfea_file, trainfea_file+'.scale.sub.out', trainfea_file+'.scale.sub.cross-acc'])

            # predict
            svm_predict(range_file, model_file, fea_file, orf_file, pre.svm_cut_value2, pre.prob_status, cds_pred_file)
            cmd = '%s -a %s -b %s --a+b -o %s' % (pre.metalocs_operate, cds_pred_file, this_nn_1_file, this_nn_1_file)
            call(cmd, shell = True, stdout = PIPE)
            rm_files([range_file, model_file, fea_file, orf_file, cds_pred_file])
        cmd = '%s -a %s -b %s --a+b -o %s' % (pre.metalocs_operate, this_nn_1_file, nn_1_file, nn_1_file)
        call(cmd, shell=True, stdout=PIPE)
        rm_files([this_nn_1_file])
    fout_nn_1.close()

    print '  Blast e-valule cut:', pre.blast_ev
    bl_file = pre.hit_file + '.bl'
    cmd = '%s -f %s -o %s -e %s' % (pre.parse_blast, pre.hit_file, bl_file, pre.blast_ev)
    call(cmd, shell=True, stdout=PIPE)

    cmd = '%s -a %s -b %s --a+b -o %s' % (pre.metalocs_operate, nn_1_file, bl_file, novel_pred_file)
    call(cmd, shell=True, stdout=PIPE)
#    print '  Final predictions store in file:', novel_pred_file
    rm_files([bl_file, bh_file, nn_1_file])
    return
    
def main():
    argc = len(sys.argv)
    argv = sys.argv
    if argc == 1:
        print parse_command.__doc__
        return
    if command_line_parser(argv, '-h', '--help', False)[0]:
        print detail_help.__doc__
        return

    # load default settings
    pre.default_settings()
    
    # parse settings from command line
    parse_command(argv)
    
    if len(pre.seqs_file) == 0:
        print 'The sequence file is not defined.'
        return
    if not os.path.exists(pre.seqs_file):
        print 'The sequence file %s is not correctly defined.' % pre.seqs_file
        return

    print_section_title('Start MetaGUN Prediction')
    print '  Input sequence file:', pre.seqs_file
    print '  Project name:', pre.project_name

    print_section_title('Sequence Binning')                         
    if pre.exists_bin_file == False:
        pre.bin_file = pre.project_name + '.bin'
        bindistri_file = pre.project_name + '.bindistri'
        print '  Sequence binning model:', pre.binmodel
        print '  Taxonomy map file:', pre.taxonomy
        cmd = '%s -m %s -s %s -t %s -b %s -d %s' % \
              (pre.bin_seqs, pre.binmodel, pre.seqs_file, pre.taxonomy, pre.bin_file, bindistri_file)
        call(cmd, shell = True, stdout = PIPE)
        print '  Binning result file:', pre.bin_file
    else:
        print '  Use given binning result:', pre.bin_file
        
    uni_pred_file = pre.project_name + '.uni-pred'
    novel_pred_file = pre.project_name + '.novel-pred'
    
    if pre.run_uni_pred:
        uni_pred(uni_pred_file)
    if pre.run_novel_pred:
        novel_pred(novel_pred_file)

    metagun_file = pre.project_name + '.metagun'
    if pre.run_uni_pred and pre.run_novel_pred:
        cmd = '%s -a %s -b %s --a+b -o %s' % (pre.metalocs_operate, novel_pred_file, uni_pred_file, metagun_file)
        call(cmd, shell=True, stdout=PIPE)
        rm_files([uni_pred_file, novel_pred_file])
    elif pre.run_uni_pred and not pre.run_novel_pred:
        os.rename(uni_pred_file, metagun_file)
    elif not pre.run_uni_pred and pre.run_novel_pred:
        os.rename(novel_pred_file, metagun_file)
    if not pre.run_metatisa:
        print_section_title('MetaGUN Prediction Over')
        print '  Final predictions output', metagun_file
        return
    
    print_section_title('TIS Refinement')
    metatisa_file = metagun_file + '.metatisa'
    cmd = '%s -s %s -l %s -b %s -p %s -f MED -o %s -c %s -w %i' % \
          (pre.metatisa, pre.seqs_file, metagun_file, pre.bin_file, pre.tismodel, metatisa_file, pre.metatisa_settings, is_win)
    call(cmd, shell=True, stdout=PIPE)
    if os.path.exists(metagun_file):
        os.remove(metagun_file)
    os.rename(metatisa_file, metagun_file)
    print_section_title('MetaGUN Prediction Over')
    print '  Final predictions output', metagun_file
    return           

if __name__=="__main__":
    main()
