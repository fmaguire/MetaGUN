
#include "feature.h"

vector<double> Feature::getfea(const string& seq, int maxlen, double upfalse_ts, double true_ts, double downfalse_ts, bool is_leftmost)
{
	vector<double> fea = SeqOp::codonEDP(seq.data(), 0, (int)seq.size());
	int seqlen = ((int)seq.size() > maxlen ? maxlen : (int)seq.size());
	fea.push_back(seqlen);

	if(is_leftmost) 
	{
		fea.push_back(upfalse_ts);
		fea.push_back(true_ts);
		fea.push_back(downfalse_ts);
	}

	return fea;
}

void Feature::foutfea(std::ofstream &fout, const std::vector<double> &fea, bool is_coding)
{
	fout << std::setprecision(4);
	fout << (is_coding ? "+1 " : "-1 ");
	for(int k = 0; k < (int)fea.size(); ++ k)
		if(fea[k] != 0)
			fout << k+1 << ':' << fea[k] << ' ';
	fout << endl;
}
