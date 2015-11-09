
#include "seqop.h"
#include "meta.h"
#include "OftenUsedFun.h"
#include <iostream>

void exit_with_help()
{
	cout<< "get-data   2012-11-14   arguments:\n"
		<< "common options:\n"
		<< "    -p, --pred             to switch on predict mode, default, OFF;\n"
		<< "    -n, --name  STRING     project name, default: 'sample';\n"
		<< "    -s, --seqs  FILE       fasta format sequence file;\n"
		<< "    -t, --tism  PATH       directory path of TIS models;\n"
		<< "    -m, --minL  INT        minimum length of ORFs. default: 60;\n"
		<< "    -M, --maxL  INT        maximun length of ORFs. default: 1500;\n"
		<< "options when pred mode is OFF:\n"
		<< "    -c, --cds  FILE        coding ORF file, for the training mode.\n"
		<< "    -q, --csq  INT         0: off; 1: output amino acid sequences according to the given locations of -ccds option;\n"
		<< "                           2: nucleotides sequences, default: 0.\n"
		<< "options when pred mode is ON:\n"
		<< "    -b, --bin  FILE        sequence binning file.\n"
		<< "    -h, --shd              to switch on labeling coding and non-coding ORFs in the shadow rule. defaule: OFF.\n"
		<< "    -1, --one              to ouput all features into one file. default: OFF\n"
		<< "        --set              to switch on the ORFsets mode. default: OFF\n"
		<< "        --sts              to output statistics of the input fragments. default: OFF.\n"
		<< "    -a, --alc              to output the locations of all the longest ORFs. default: OFF.\n"
		<< "    -u, --unix             to run the program under unix platform, defalut, OFF.\n"
		<< endl;
	exit(-1);
}

int main(int argc, char* argv[])
{
	// default settings
	string name = "sample";
	string bin_file, cds_file, tis_path, seqs_file;

	bool is_pred = false, is_bin = false, is_outaa = false, is_outnu = false, is_unix = false;
	bool is_outsts = false, is_outalc = false, is_ORFsets = false, is_shadow = false, is_one = false;
	int min_ORF_len = 60, max_ORF_len = 1200, is_csq = 0;

	// parse command line
	if(argc == 1)
		exit_with_help();
	pair<bool, string> opt = command_line_parser(argc, argv, "-p", "--pred", false);
	is_pred = opt.first;
	opt = command_line_parser(argc, argv, "-n", "--name", true);
	if(opt.first) name = opt.second;
	opt = command_line_parser(argc, argv, "-s", "--seqs", true);
	if(opt.first) seqs_file = opt.second;
	opt = command_line_parser(argc, argv, "-t", "--tism", true);
	if(opt.first) tis_path = opt.second;
	opt = command_line_parser(argc, argv, "-m", "--minL", true);
	if(opt.first) min_ORF_len = atoi(opt.second.data());
	opt = command_line_parser(argc, argv, "-M", "--maxL", true);
	if(opt.first) max_ORF_len = atoi(opt.second.data());
	opt = command_line_parser(argc, argv, "-c", "--cds", true);
	if(opt.first) cds_file = opt.second;
	opt = command_line_parser(argc, argv, "-q", "--csq", true);
	if(opt.first) is_csq = atoi(opt.second.data());
	opt = command_line_parser(argc, argv, "-b", "--bin", true);
	if(opt.first) {
		bin_file = opt.second;
		is_bin = true;
	}
	is_shadow = command_line_parser(argc, argv, "-h", "--shd", false).first;
	is_ORFsets = command_line_parser(argc, argv, "", "--set", false).first;
	is_outsts = command_line_parser(argc, argv, "", "--sts", false).first;
	is_outalc = command_line_parser(argc, argv, "-a", "--alc", false).first;
	is_unix = command_line_parser(argc, argv, "-u", "--unix", false).first;

	if(seqs_file.empty())
	{
		cout << "sequences file is not defined." << endl;
		return -1;
	}

	MetaSeqs meta;
	meta.min_ORF_len = min_ORF_len;
	meta.max_ORF_len = max_ORF_len;
	meta.is_unix = is_unix;
	meta.get_metaseqs(seqs_file);
	if(is_outsts)
	{
		meta.get_seqs_info(name);
		return 0;
	}
	meta.find_all_LORFs();
	if(is_outalc)
	{
		meta.fout_ORF_locs(string(name + ".alc"));
		return 0;
	}
	if(is_outaa)
	{
		meta.fout_ORF_seqs(string(name + ".alsa"), true);
		return 0;
	}
	if(is_outnu)
	{
		meta.fout_ORF_seqs(string(name + ".alsn"), false);
		return 0;
	}
	if(is_csq != 0)
	{
		if(cds_file.empty())
		{
			cout << "the coding sequences file is not specified." << endl;
			return -1;
		}
		cout << "minimal ORF length: " << meta.min_ORF_len << endl;
		meta.label_ORFs(cds_file);
		if(is_csq == 1)
			meta.fout_cod_seqs(string(name + ".csqa"), true);
		else
			meta.fout_cod_seqs(string(name + ".csqn"), false);
		return 0;
	} 

	if(bin_file.empty())
	{
		cout << "the binning file is not specified." << endl;
		return -1;
	}	
	if(tis_path.empty())
	{
		cout << "the tis path is not specified. " << endl;
		return -1;
	}
	meta.get_bin_info(bin_file);

	if(is_pred)
	{
		if(is_ORFsets)
		{
			meta.find_all_ORFsets();
			meta.tis_scoring_ORFsets(tis_path);
		}
		else
			meta.tis_scoring(tis_path);
		/*meta.genus_binreads.size() == 1 ? is_bin = false : is_bin = true;*/
		if(is_one)
			is_bin = false;
		if(is_bin)
			meta.get_allfea(name, 2, is_ORFsets);
		else
			meta.get_allfea(name, is_ORFsets);
	}
	else
	{
		if(cds_file.empty())
		{
			cout << "the coding sequence file is not specified. " << endl;
			return -1;
		}
		meta.tis_scoring(tis_path);
		if(is_shadow) {
			meta.label_ORFs_shadow_rule(cds_file);
		}
		else
			meta.label_ORFs(cds_file);
		cout << "is_bin: " << is_bin << endl;
		if(is_bin)
			meta.get_trainfea(name, 2);
		else
			meta.get_trainfea(name);
	}
	return 0;
}
