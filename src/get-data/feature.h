#ifndef FEATURE_H
#define FEATURE_H

#include "seqop.h"
#include <fstream>
#include <iomanip>

namespace Feature
{
	vector<double> getfea(const string& seq, int maxlen, double upfalse_ts, double true_ts, double downfalse_ts, bool is_leftmost);
	void foutfea(ofstream& fout, const vector<double>& fea, bool is_coding);
};

#endif
