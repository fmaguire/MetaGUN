#!/usr/bin/env python

__author__ = "Liu Yong-Chu"
__version__ = "Revision: 3.0"
__create__ = "2012-05-23"
__update__ = "2012-06-11"

import sys, os
from OftenUsedFun import *

def default_settings():
     """Default settings for MetaGUN."""

     script_file_path = __file__[0:__file__.rfind(os.sep)]
     binary_file_path = os.sep.join([script_file_path, '..', 'bin'])
     is_win = False
     if 'win' in sys.platform:
          is_win = True

     ########           Default programs
     global svm_train, svm_scale, svm_predict
     global bin_seqs, get_data, extr_pred_cds, rpsblast, parse_blast, metalocs_operate, metatisa
     global train_cds_model_py, subset_py, grid_py
     
     svm_train = os.sep.join([binary_file_path, 'svm-train'])
     svm_scale = os.sep.join([binary_file_path, 'svm-scale'])
     svm_predict = os.sep.join([binary_file_path, 'svm-predict'])
     
     bin_seqs = os.sep.join([binary_file_path, 'bin-seqs'])
     get_data = os.sep.join([binary_file_path, 'get-data'])
     extr_pred_cds = os.sep.join([binary_file_path, 'extr-pred-cds'])
     rpsblast = os.sep.join([binary_file_path, 'rpsblast'])
     parse_blast = os.sep.join([binary_file_path, 'parse-blast'])
     metalocs_operate = os.sep.join([binary_file_path, 'metalocs-operate'])
     metatisa = os.sep.join([binary_file_path, 'metatisa'])
     
     train_cds_model_py = os.sep.join([script_file_path, 'train-cds-model.py'])
     subset_py = os.sep.join([script_file_path, 'subset.py'])
     grid_py = os.sep.join([script_file_path, 'grid-xtest.py'])

     if is_win:
          check_file_existence([svm_train+'.exe', svm_scale+'.exe', svm_predict+'.exe', bin_seqs+'.exe', get_data+'.exe'])
          check_file_existence([extr_pred_cds+'.exe', rpsblast+'.exe', parse_blast+'.exe', metalocs_operate+'.exe'])
     else:
          check_file_existence([svm_train, svm_scale, svm_predict, bin_seqs, get_data])
          check_file_existence([extr_pred_cds, rpsblast, parse_blast, metalocs_operate])
     check_file_existence([train_cds_model_py, subset_py, grid_py])

     ########           Default parameters
     global project_name, taxonomy, binmodel, cdsmodel, tismodel, blast_db, metatisa_settings
     global min_orf_len, max_orf_len, svm_cut_value, svm_cut_value2, svm_sub_size
     global exists_bin_file, exists_hit_file, ORFsets_status, prob_status, is_help
     global bin_file, seqs_file, hit_file, blast_ev, seeds_ev
     global run_uni_pred, run_novel_pred, run_metatisa

     dat_file_path = os.sep.join([script_file_path, '..', 'dat'])
     project_name = 'sample'
#     taxonomy = os.sep.join([dat_file_path, 'binmodel', 'test.bin-map'])
#     binmodel = os.sep.join([dat_file_path, 'binmodel', 'test.binmodel'])
     taxonomy = os.sep.join([dat_file_path, 'binmodel', '261-genomes.bin-map'])
     binmodel = os.sep.join([dat_file_path, 'binmodel', '261-genomes.k8.binmodel'])
     cdsmodel = os.sep.join([dat_file_path, 'cdsmodel'])
     tismodel = os.sep.join([dat_file_path, 'tismodel'])
     blast_db = os.sep.join([dat_file_path, 'Cdd', 'Cdd'])
     metatisa_settings = os.sep.join([dat_file_path, 'tismodel', 'metatisa-settings.txt'])
     
     min_ORF_len = 60
     max_ORF_len = 1500
     svm_cut_value = 0.5
     svm_cut_value2 = 0.5
     svm_sub_size = 10000

     exists_bin_file = False
     exists_hit_file = False
     ORFsets_status = 0
     prob_status = 1
     is_help = False
     
     run_uni_pred = True
     run_novel_pred = True
     run_metatisa = True
     
     bin_file = ''
     seqs_file = ''
     hit_file = ''
     blast_ev = 1e-10
     seeds_ev = 1e-40

if __name__ == "__main__":
     default_settings()
     print "\n----  All variables  ----"
     print "\n".join(["\t%s %s" % (k.ljust(25), v) for k, v in vars().items()])
else:
     default_settings()
