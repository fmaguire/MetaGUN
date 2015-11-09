
#include "tis.h"

void TIS::initiate_settings()
{
	max_order = 2;
	this_order = 2;
	upbps = 50;
	downbps = 15;
	headnum_cutoff = 200;
	minimum_upstream = 10;
}

void TIS::set_codons()
{
	Starts.erase(Starts.begin(), Starts.end());
	Stops.erase(Stops.begin(), Stops.end());

	// ATG, CTG, GTG, TTG
	Starts.insert(14), /*Starts.insert(30), */Starts.insert(46), Starts.insert(62);

	// TAA, TAG, TGA
	Stops.insert(48), Stops.insert(50), Stops.insert(56);
}

Matrix_T<double> TIS::get_matrix(const std::string &matrix_file)
{
	ifstream fin(matrix_file.data());
	if(!fin) { cout << "can't open " << matrix_file << endl; exit(-1); }

	int lin = 0, col = 0;
	string line;

	fin >> lin >> col;
	getline(fin, line);

	Matrix_T<double> ma(lin, col, 0.0);
	for(int c = 0; c < col; ++ c)
	{
		for(int l = 0; l < lin; ++ l)
			fin >> ma(l, c);
		getline(fin, line);
	}
	fin.close();
	return ma;
}

void TIS::get_prior_paras(const std::string &genus)
{
	// for windows
	string path = prior_para_path + '\\' + genus + '\\';
	// for unix
//	string path = prior_para_path + '/' + genus + '/';
	string prior_prob_file = path + genus + ".prioriProb";
	string prior_true_file = path + genus + ".trueMa";
	string prior_upfalse_file = path + genus + ".ufMa";
	string prior_downfalse_file = path + genus + ".dfMa";

	ifstream fin_prior_prob(prior_prob_file.data());
	if(!fin_prior_prob) { cout << "can't open " << prior_prob_file << endl; exit(-1); }
	fin_prior_prob >> prob_true >> prob_upfalse >> prob_downfalse;
	fin_prior_prob.close();

	trueTIS = get_matrix(prior_true_file);
	upfalse = get_matrix(prior_upfalse_file);
	downfalse = get_matrix(prior_downfalse_file);
}

void TIS::read_prior_tismodel(const std::string &tismodel_file)
{
	std::ifstream fin(tismodel_file.data());
	if(!fin) { cout << "can't open " << tismodel_file << endl; exit(-1); }
	std::string line;
	while(!fin.eof()) {
		line.clear();
		getline(fin, line);
		if(line.empty())
			continue;
		if(line == "-------Prior probabilites------") {
			fin >> prob_true >> prob_upfalse >> prob_downfalse >> cutoff;
			getline(fin, line);
		} 
		else if(line == "------True TIS Markov model------") {
			int lin=0, col=0;
			fin >> lin >> col;
			getline(fin, line);
			Matrix_T<double> mm(lin, col, 0.0);
			for(int c=0; c<col; ++c) {
				for(int l=0; l<lin; ++l)
					fin >> mm(l, c);
				getline(fin, line);
			}
			trueTIS = mm;
		}
		else if(line == "------Upfalse TIS Markov model------") {
			int lin=0, col=0;
			fin >> lin >> col;
			getline(fin, line);
			Matrix_T<double> mm(lin, col, 0.0);
			for(int c=0; c<col; ++c) {
				for(int l=0; l<lin; ++l)
					fin >> mm(l, c);
				getline(fin, line);
			}
			upfalse = mm;
		}
		else if(line == "------Downfalse TIS Markov model------") {
			int lin=0, col=0;
			fin >> lin >> col;
			getline(fin, line);
			Matrix_T<double> mm(lin, col, 0.0);
			for(int c=0; c<col; ++c) {
				for(int l=0; l<lin; ++l)
					fin >> mm(l, c);
				getline(fin, line);
			}
			downfalse = mm;
		}
	}
	fin.close();
}

void TIS::scoring_XTG(const char* ssseq, double& t, double& d, double& u, int start)
{
	// set to zero
	int beg = 0;
	t = 0, d = 0, u = 0;

	int index = 0;
	int f = 0;
	// find the index of first (ORDER-1) nucleotide(s)
	for( ; f < this_order ; ++f ){
		index += (int)pow(4.0,this_order-f)*(ssseq[beg+f]-'0');
	}
	// get the score of the frist (this_order-1) nucleotides(s)
	for( f = 0; f < 4; ++f ){
		t += exp(trueTIS(index+f, start)); 
		d += exp(downfalse(index+f, start));
		u += exp(upfalse(index+f, start));
	}
	// multipe the prior probabilities, as plus for log
	t = log(t) + log(prob_true),  d = log(d) + log(prob_downfalse), u = log(u) + log(prob_upfalse);

	// formulating the this_order-th Markov model
	int i = 0;
	for( ; i+start < trueTIS.getColum(); ++i){
		int index = 0;
		int f = 0;
		for( ; f < this_order ; ++f ){
			index += (int)pow(4.0,this_order-f)*(ssseq[beg+i+f]-'0');
		}
		double pt = 0, pt2 = 0, pd = 0, pu = 0;
		for( f = 0; f < 4; ++f ){
			pt += exp(trueTIS(index+f, i+start)); 
			pd += exp(downfalse(index+f, i+start));
			pu += exp(upfalse(index+f, i+start));
		}
		t += trueTIS(index+(ssseq[beg+i+this_order]-'0'), i+start) - log(pt);
		d += downfalse(index+(ssseq[beg+i+this_order]-'0'), i+start) - log(pd);
		u += upfalse(index+(ssseq[beg+i+this_order]-'0'), i+start) - log(pu);
	}
	t = exp(t), d = exp(d); u = exp(u);
	double sum = t+d+u;
	t /= sum, d /= sum, u /= sum;
}

void TIS::scoring_XTG(const char* ssseq, double& t, double& d, double& u)
{
	// set to zero
	int beg = 0;
	t = 0, d = 0, u = 0;

	int index = 0;
	int f = 0;
	// find the index of first (ORDER-1) nucleotide(s)
	for( ; f < this_order ; ++f ){
		index += (int)pow(4.0,this_order-f)*(ssseq[beg+f]-'0');
	}
	
	// get the score of the frist (this_order-1) nucleotides(s)
	for( f = 0; f < 4; ++f ){
		t += exp(trueTIS(index+f,0)); 
		d += exp(downfalse(index+f,0));
		u += exp(upfalse(index+f,0));
	}

	// multipe the prior probabilities, as plus for log
	t = log(t) + log(prob_true),  d = log(d) + log(prob_downfalse), u = log(u) + log(prob_upfalse);

	// formulating the this_order-th Markov model
	int i = 0;
	for( ; i < trueTIS.getColum(); ++i ){
		int index = 0;
		int f = 0;
		for( ; f < this_order ; ++f ){
			index += (int)pow(4.0,this_order-f)*(ssseq[beg+i+f]-'0');
		}
		double pt = 0, pt2 = 0, pd = 0, pu = 0;
		for( f = 0; f < 4; ++f ){
			pt += exp(trueTIS(index+f,i)); 
			pd += exp(downfalse(index+f,i));
			pu += exp(upfalse(index+f,i));
		}
		t += trueTIS(index+(ssseq[beg+i+this_order]-'0'),i) - log(pt);
		d += downfalse(index+(ssseq[beg+i+this_order]-'0'),i) - log(pd);
		u += upfalse(index+(ssseq[beg+i+this_order]-'0'),i) - log(pu);
	}
	t = exp(t), d = exp(d); u = exp(u);
	double sum = t+d+u;
	t /= sum, d /= sum, u /= sum;
}

