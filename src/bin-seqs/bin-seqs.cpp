
#include <iostream>
#include <string>
#include "bin.h"
#include "OftenUsedFun.h"
using namespace std;

void exit_with_help()
{
	cout << "bin-seqs    2012-11-14    arguments:\n"
		 << "    -m  FILE      sequence binning model.\n"
		 << "    -t  FILE      sequence binning taxonomic relationships.\n"
		 << "    -s  FILE      input sequences, FASTA format.\n"
		 << "    -b  STRING    sequence binning results. default, bins.txt.\n"
		 << "    -d  STRING    distribution of binning results. default, bindistri.txt.\n"
		 << endl;
	exit(-1);
}

int main(int argc, char* argv[])
{
	// default settings
	string binmodel_file, binmap_file, seqs_file;
	int kmer_len = 0;
	string bin_file = "bins.txt", bindistri_file = "bindistri.txt";
	bool is_binary = true;

	// parse command line
	if(argc == 1) exit_with_help();
	pair<bool, string> opt = command_line_parser(argc, argv, "-m", "", true);
	if(opt.first) {
		binmodel_file = opt.second;
	} else {
		cout << "sequence binning model is not defined." << endl; 
		return -1;
	}

	opt = command_line_parser(argc, argv, "-t", "", true);
	if(opt.first) {
		binmap_file = opt.second;
	} else {
		cout << "sequence binning taxonomic relationships is not defined." << endl;
		return -1;
	}

	opt = command_line_parser(argc, argv, "-s", "", true);
	if(opt.first) {
		seqs_file = opt.second;
	} else {
		cout << "sequence file is not defined." << endl;
		return -1;
	}

	opt = command_line_parser(argc, argv, "-b", "", true);
	if(opt.first)
		bin_file = opt.second;
	opt = command_line_parser(argc, argv, "-d", "", true);
	if(opt.first)
		bindistri_file = opt.second;

	map<string, string> seqs = read_fasta(seqs_file);
	map<string, vector<float> > binmodel = get_binmodel(binmodel_file, kmer_len, is_binary);
	map<string, pair<string, string> > binmap = get_binmap(binmap_file);

	map<string, string> reads_bing;
	reads_bing = bin_seqs(binmodel, seqs, kmer_len);
	output_bin(reads_bing, binmap, bin_file, bindistri_file);
	return 0;
}

