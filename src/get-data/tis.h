//  author: liuyc
// program: for getting TIS score
// history: created -- 2010/4/7
#include "seqop.h"
#include "matrix.h"
#include <fstream>
#include <iostream>
#include <set>
using namespace std;

#ifndef TIS_H
#define TIS_H

class TIS
{
public:
	TIS()
	{
		upbps = downbps = max_order = this_order = headnum_cutoff = -1;
		pT = 0, pFU = 0, pFD = 0, cutoff = 0;
	}

	// settings
	double pT, pFU, pFD, cutoff;
	int upbps, downbps, max_order, this_order;
	int headnum_cutoff;
	int minimum_upstream;
	string prior_para_path;
	set<int> Starts, Stops;

	void set_codons();
	void get_para_path(const string& path) { prior_para_path = path; }

	// parameters to be calculate
	Matrix_T<double> trueTIS, upfalse, downfalse;
	int max_cand_num;
	double prob_upfalse, prob_downfalse, prob_true;

	// operates
	void initiate_settings();
	Matrix_T<double> get_matrix(const string& matrix_file);
	void get_prior_paras(const string& genus);
	void read_prior_tismodel(const string& tismodel_file);
	void scoring_XTG(const char* ssseq, double& t, double& d, double& u, int start);
	void scoring_XTG(const char* ssseq, double& t, double& d, double& u);
};

#endif
