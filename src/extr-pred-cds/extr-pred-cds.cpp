// created:  2009-5-27
// revised:  2009-5-27
// extract predicted cds and get the corresponding svm label and probabilities

#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <iomanip>
using namespace std;

class CDS 
{
public:
	CDS() { start = stop = phase = complt = label = 0; 
	prob = 0; strand = '+'; }
	int start, stop, phase, complt;
	char strand;
	int label;
	double prob;
};

std::map<string, vector<CDS> > get_cdss_b1(const string& pred_file, const string& ORF_file)
{
	// get label
	int true_num = 0, false_num = 0;
	std::ifstream fin;
	fin.open(pred_file.data());
	if(!fin){cout<<"can't open "<<pred_file<<endl; exit(0);}

	std::vector<std::pair<int, double> > label_probs;

	string line;
	getline(fin, line);
	istringstream is(line);
	string label;
	int a = 0, b = 0;
	is >> label >> a >> b;

	int col = 1;
	if(b == 1) col = 2;
	while(!fin.eof())
	{
		int label = 0;
		double prob = 0;
		fin >> label;
		for(int i=0; i<col; ++i)
			fin >> prob;
		getline(fin, line);
		if(label == 1){
			++ true_num;
		} else
			++false_num;
		label_probs.push_back(std::make_pair(label, prob));
	}
	fin.close();

	// map to the ORFs
	std::ifstream fin1;
	fin1.open(ORF_file.data());
	if(!fin1.good()) { cout << "open " << ORF_file << " failed." << endl; exit(0); }

	std::map<string, vector<CDS> > CDSs;
	int index = 0;
	int trueN = 0;
	while(!fin1.eof()) {
		string line;
		getline(fin1, line);
		if(line.empty()) continue;

		for(int k = 0; k < (int)line.size(); ++k)
			if(line[k] == '|')
				line[k] = ' ';
		string pid;
		CDS cl;
		istringstream is(line);
		is >> pid >> cl.start >> cl.stop >> cl.strand >> cl.phase >> cl.complt;
		cl.label = label_probs[index].first;
		cl.prob = label_probs[index].second;

		std::map<string, vector<CDS> >::iterator iter = CDSs.find(pid);
		if(iter == CDSs.end())
			CDSs[pid].push_back(cl);
		else
		{
			bool exist = false;
			for(int a = 0; a < (int)iter->second.size(); ++ a)
			{
				CDS& tmpc = iter->second[a];
				if( (tmpc.strand == '+' && cl.strand == '+' && tmpc.stop == cl.stop)
					|| (tmpc.strand == '-' && cl.strand == '-' && tmpc.start == cl.start) )
				{
					exist = true;
					if(cl.prob > tmpc.prob)
						tmpc = cl;
					break;
				}
			}
			if(!exist)
				iter->second.push_back(cl);
		}
		++ index;
	}
	cout << std::left << setw(8) << index << " ORFs in " << CDSs.size() << " sequence reads" << endl;
	return CDSs;
}

std::map<string, vector<CDS> > get_cdss_b0(const string& pred_file, const string& ORF_file)
{
	// get label
	int true_num = 0, false_num = 0;
	std::ifstream fin(pred_file.data());
	if(!fin){cout<<"can't open "<<pred_file<<endl; exit(0);}

	std::vector<int> labels;
	while(!fin.eof())
	{
		int label = 0;
		fin >> label;
		if(label != 1 && label != -1)
			continue;
		labels.push_back(label);
		if(label == 1)
			++ true_num;
		else 
			++ false_num;
	}
	fin.close();

	// map to the ORFs
	std::ifstream fin1;
	fin1.open(ORF_file.data());
	if(!fin1.good()) { cout << "open " << ORF_file << " failed." << endl; exit(0); }

	std::map<string, vector<CDS> > CDSs;
	int index = 0;
	int trueN = 0;
	while(!fin1.eof()) {
		string line;
		getline(fin1, line);
		if(line.empty()) continue;

		for(int k = 0; k < (int)line.size(); ++k)
			if(line[k] == '|')
				line[k] = ' ';
		string pid;
		CDS cl;
		istringstream is(line);
		is >> pid >> cl.start >> cl.stop >> cl.strand >> cl.phase >> cl.complt;
		cl.label = labels[index];

		std::map<string, vector<CDS> >::iterator iter = CDSs.find(pid);
		if(iter == CDSs.end())
			CDSs[pid].push_back(cl);
		else
		{
			bool exist = false;
			for(int a = 0; a < (int) iter->second.size(); ++ a)
			{
				CDS& tmpc = iter->second[a];
				if( (tmpc.strand == '+' && cl.strand == '+' && tmpc.stop == cl.stop)
					|| (tmpc.strand == '-' && cl.strand == '-' && tmpc.start == cl.start) )
				{
					exist =true;
					if(tmpc.stop - tmpc.start + 1 < cl.stop - cl.start + 1)
						tmpc = cl;
					break;
				}
			}
			if(!exist)
				iter->second.push_back(cl);
		}
		++ index;
	}
	cout << std::left << setw(8) << index << " ORFs in " << CDSs.size() << " sequence reads" << endl;
	cout<< setw(8) << true_num << " predicted to be CDS" << endl
		<< setw(8) << false_num << " predicted to be Non-Coding" << endl;
	return CDSs;
}

void adjust_lable(map<string, vector<CDS> >& CDSs, const double cut_value = 0.5)
{
	int true_num = 0, false_num = 0;
	map<string, vector<CDS> >::iterator it = CDSs.begin();
	for(; it != CDSs.end(); ++it)
	{
		for(unsigned int k = 0; k < it->second.size(); ++k) {
			if(it->second[k].prob >= cut_value) {
				it->second[k].label = 1;
				++ true_num;
			} else {
				it->second[k].label = -1;
				++ false_num;
			}
		}
	}
	cout << std::left;
	cout<< setw(8) << true_num << " predicted to be CDS" << endl
		<< setw(8) << false_num << " predicted to be Non-Coding" << endl;
}

void output(map<string, vector<CDS> >& true_cdss, const string& file)
{
	ofstream fout(file.data());
	if(!fout) { cout << "can't write " << file << endl; exit(0); }

	fout << std::setprecision(4);
	map<string, vector<CDS> >::iterator iter = true_cdss.begin();
	for(; iter != true_cdss.end(); ++iter) {
		fout << '>' << iter->first << endl;
		for(int k = 0; k < (int) iter->second.size(); ++k) {
			if(iter->second[k].label == -1)
				continue;
			fout << setw(8) << iter->second[k].start << setw(8) << iter->second[k].stop
				<< setw(4) << iter->second[k].strand << setw(3) << iter->second[k].phase 
				<< setw(3) << iter->second[k].complt / 10 << iter->second[k].complt % 10 
				<< setw(15) << iter->second[k].prob << endl;
		}
	}
}

void help()
{
	cout << "extr-pred-cds    2012-11-14    arguments:" << endl
		 << "extr-pred-cds  [options]  allORF-file  svm-predict-file" << endl
		 << "Options:" << endl
		 << "    -c FLOAT     cut value of the identification, Default, 0.5. ORFs with a probability\n"
		 << "                 higher than that are regarded as CDSs, while others are non-coding.\n"
		 << "    -b INT       0/1, 1 for the probability svm model, 0 for ordinary svm-model, default, 1.\n"
		 << "    -o FILE      out file name, Default, CDSs.txt\n";
	exit(0);
}

int main(int argc, char* argv[])
{
	double cut_value = 0.5;
	string out_file = "CDSs.txt";
	string ORF_file, pred_file;
	int b = 1;

	if(argc < 3 ) help();
	int i = 1;
	for(; i < argc - 2; ++i)
	{
		if(argv[i][0] != '-') continue;
		if(++i >= argc - 2) {
			cout << "option " << argv[i-1] << " is not specified" << endl;
			return 0;
		}
		switch(argv[i-1][1])
		{
		case 'c': cut_value = atof(argv[i]); 
			break;
		case 'o': out_file = string(argv[i]); 
			break;
		case 'b': b = atoi(argv[i]); 
			break;
		}
	}
	ORF_file = argv[i];
	pred_file= argv[i+1];

	if(b == 1)
	{
		cout<< "cut_value: " << cut_value << endl;
		map<string, vector<CDS> > CDSs = get_cdss_b1(pred_file, ORF_file);
		adjust_lable(CDSs, cut_value);
		output(CDSs, out_file);
		return 0;
	}
	else if(b == 0)
	{
		map<string, vector<CDS> > CDSs = get_cdss_b0(pred_file, ORF_file);
		output(CDSs, out_file);
		return 0;
	}
	else
	{
		cout << "No " << b << " option for -b" << endl; 
		return -1;
	}
}
