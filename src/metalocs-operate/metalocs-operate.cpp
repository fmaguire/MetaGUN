
// Operates of two metalocs sets,
// version 2.0:
// date: 2010/06/21
// author: Liu Yong-Chu

#include "OftenUsedFun.h"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

class MetaLoc
{
public:
	MetaLoc() {start = stop = phase = complt = -1; strand = '+'; prob = 0.0; descript = ""; }
	bool ispos() { if(strand == '+') return true; else return false; }
	int start, stop, phase, complt;
	double prob;
	char strand;
	string descript;
};

map<string, vector<MetaLoc> > get_metalocs(const string& file)
{
	std::ifstream fin(file.data());
	if(!fin) { cout << "can't open "<<file<<endl; exit(0); }

	map<string, vector<MetaLoc> > metalocs;
	int locs_num = 0;
	string line;
	getline(fin, line);
	while(!fin.eof()) 
	{
		if(line.find('>') != std::string::npos)
		{
			string pid = line.substr(line.find('>')+1);
			getline(fin, line);
			while(line.find('>') == std::string::npos && !fin.eof())
			{
				vector<string> vv = split_string(line);
				MetaLoc c;
				c.start = atoi(vv[0].data()), c.stop = atoi(vv[1].data()), c.strand = vv[2][0];
				metalocs[pid].push_back(c);
				++ locs_num;
				getline(fin, line);
			}
		}
		else
			getline(fin, line);
	}
	fin.close();
	return metalocs;
}

map<string, vector<MetaLoc> > inter_metaloc(map<string, vector<MetaLoc> >& metaloc_a, 
											map<string, vector<MetaLoc> >& metaloc_b)
{
	map<string, vector<MetaLoc> > inter;
	map<string, vector<MetaLoc> >::iterator itera = metaloc_a.begin();
	for(; itera != metaloc_a.end(); ++ itera)
	{
		map<string, vector<MetaLoc> >::iterator iterb = metaloc_b.find(itera->first);
		if(iterb == metaloc_b.end())
			continue;
		for(int a = 0; a < (int)itera->second.size(); ++ a)
		{
			MetaLoc& ma = itera->second[a];
			bool bothexist = false;
			for(int b = 0; b < (int)iterb->second.size(); ++ b)
			{
				MetaLoc& mb = iterb->second[b];
				if((ma.ispos() && mb.ispos() && ma.stop == mb.stop)
					|| (!ma.ispos() && !mb.ispos() && ma.start == mb.start))
				{
					bothexist = true;
					break;
				}
			}
			if(bothexist)
				inter[itera->first].push_back(ma);
		}
	}
	return inter;
}

map<string, vector<MetaLoc> > union_metaloc(map<string, vector<MetaLoc> >& metaloc_a, 
											map<string, vector<MetaLoc> >& metaloc_b)
{
	map<string, vector<MetaLoc> > union_metalocs = metaloc_b;
	map<string, vector<MetaLoc> >::iterator itera = metaloc_a.begin();
	for(; itera != metaloc_a.end(); ++ itera)
	{
		map<string, vector<MetaLoc> >::iterator iterb = metaloc_b.find(itera->first);
		if(iterb == metaloc_b.end())
		{
			union_metalocs.insert(std::make_pair(itera->first, itera->second));
			continue;
		}
		for(int a = 0; a < (int)itera->second.size(); ++ a)
		{
			MetaLoc& ma = itera->second[a];
			bool only_a_exist = true;
			for(int b = 0; b < (int)iterb->second.size(); ++ b)
			{
				MetaLoc& mb = iterb->second[b];
				if((ma.ispos() && mb.ispos() && ma.stop == mb.stop)
					|| (!ma.ispos() && !mb.ispos() && ma.start == mb.start))
				{
					only_a_exist = false;
					break;
				}
			}
			if(only_a_exist)
				union_metalocs[itera->first].push_back(ma);
		}
	}
	return union_metalocs;
}

map<string, vector<MetaLoc> > differ_metaloc(map<string, vector<MetaLoc> >& metaloc_a, 
											 map<string, vector<MetaLoc> >& metaloc_b)
{
	map<string, vector<MetaLoc> > a_minus_b;
	map<string, vector<MetaLoc> >::iterator itera = metaloc_a.begin();
	for(; itera != metaloc_a.end(); ++ itera)
	{
		map<string, vector<MetaLoc> >::iterator iterb = metaloc_b.find(itera->first);
		if(iterb == metaloc_b.end())
		{
			a_minus_b.insert(std::make_pair(itera->first, itera->second));
			continue;
		}

		for(int a = 0; a < (int)itera->second.size(); ++ a)
		{
			MetaLoc& ma = itera->second[a];
			bool only_a_exist = true;
			for(int b = 0; b < (int)iterb->second.size(); ++ b)
			{
				MetaLoc& mb = iterb->second[b];
				if((ma.ispos() && mb.ispos() && ma.stop == mb.stop)
					|| (!ma.ispos() && !mb.ispos() && ma.start == mb.start))
				{
					only_a_exist = false;
					break;
				}
			}
			if(only_a_exist)
				a_minus_b[itera->first].push_back(ma);
		}
	}
	return a_minus_b;
}

void comparison(map<string, vector<MetaLoc> >& metaloc_a, map<string, vector<MetaLoc> >& metaloc_b)
{
	double total_cds_a = 0, total_cds_b = 0, exist_ab = 0;
	map<string, vector<MetaLoc> >::iterator itera = metaloc_a.begin();
	for(; itera != metaloc_a.end(); ++ itera)
	{
		total_cds_a += (double) itera->second.size();
		map<string, vector<MetaLoc> >::iterator iterb = metaloc_b.find(itera->first);
		if(iterb == metaloc_b.end())
			continue;
		for(int a = 0; a < (int)itera->second.size(); ++ a)
		{
			MetaLoc& ma = itera->second[a];
			bool bothexist = false;
			for(int b = 0; b < (int)iterb->second.size(); ++ b)
			{
				MetaLoc& mb = iterb->second[b];
				if((ma.ispos() && mb.ispos() && ma.stop == mb.stop)
					|| (!ma.ispos() && !mb.ispos() && ma.start == mb.start))
				{
					bothexist = true;
					break;
				}
			}
			if(bothexist)
				++ exist_ab;
		}
	}

	map<string, vector<MetaLoc> >::iterator iterb = metaloc_b.begin();
	for(; iterb != metaloc_b.end(); ++ iterb)
		total_cds_b += (double) iterb->second.size();

	double only_exist_a = total_cds_a - exist_ab;
	double only_exist_b = total_cds_b - exist_ab;
	double sn = exist_ab / total_cds_a * 100;
	double sp = exist_ab / total_cds_b * 100; 
	double ha = 2 * sn * sp / (sn + sp);
//	cout << "#pred" << '\t' << "#benc" << '\t' << "#true" << '\t' << "#false" 
//		<< '\t' << "#miss" << '\t' << "%sn" << '\t' << "%sp" << '\t' << "%ha" << endl;
	cout << total_cds_b << '\t' << total_cds_a << '\t' << exist_ab << '\t' << only_exist_b 
		<< '\t' << only_exist_a << std::setprecision(4) << '\t' << sn << '\t' << sp << '\t' << ha << endl;
}

void get_bin_info(const string& binfile, map<string, pair<string, string> >& readbin,
				  map<string, vector<string> >& genus_binreads, map<string, vector<string> >& group_binreads)
{
	ifstream fin(binfile.data());
	if(!fin) { cout << "can't open " << binfile << endl; exit(-1); }

	readbin.erase(readbin.begin(), readbin.end());
	while(!fin.eof())
	{
		string id, genus, group, line;
		fin >> id >> genus >> group;
		getline(fin, line);
		if(!id.empty()) 
		{
			id = id.substr(1);
			readbin[ id ] = std::make_pair(genus, group);
		}
	}
	fin.close();

	genus_binreads.erase(genus_binreads.begin(), genus_binreads.end());
	group_binreads.erase(group_binreads.begin(), group_binreads.end());
	map<string, pair<string, string> >::iterator iterread = readbin.begin();
	for(; iterread != readbin.end(); ++ iterread)
	{
		genus_binreads[iterread->second.first].push_back(iterread->first);
		group_binreads[iterread->second.second].push_back(iterread->first);
	}
	cout << '\t' << readbin.size() << " reads are classified into " << genus_binreads.size() 
		<< " genus and " << group_binreads.size() << " groups." << endl;
}

void comparison(map<string, vector<MetaLoc> >& metaloc_a, map<string, vector<MetaLoc> >& metaloc_b,
				const string& bin_file)
{
	map<string, pair<string, string> > readbin;
	map<string, vector<string> > genus_binreads, group_binreads;
	get_bin_info(bin_file, readbin, genus_binreads, group_binreads);

	map<string, vector<string> >::iterator iter_g = group_binreads.begin();
	for(; iter_g != group_binreads.end(); ++ iter_g)
	{
		double total_cds_a = 0, total_cds_b = 0, exist_ab = 0;
		for(int gi = 0; gi < (int)iter_g->second.size(); ++ gi)
		{
			map<string, vector<MetaLoc> >::iterator itera = metaloc_a.find(iter_g->second[gi]);
			map<string, vector<MetaLoc> >::iterator iterb = metaloc_b.find(iter_g->second[gi]);

			if(itera == metaloc_a.end() && iterb == metaloc_b.end())
				continue;
			if(itera == metaloc_a.end() && iterb != metaloc_b.end())
			{
				total_cds_b += (int)iterb->second.size();
				continue;
			}
			if(itera != metaloc_a.end() && iterb == metaloc_b.end())
			{
				total_cds_a += (int)itera->second.size();
				continue;
			}

			total_cds_a += (int)itera->second.size();
			total_cds_b += (int)iterb->second.size();

			for(int a = 0; a < (int)itera->second.size(); ++ a)
			{
				MetaLoc& ma = itera->second[a];
				bool bothexist = false;
				for(int b = 0; b < (int)iterb->second.size(); ++ b)
				{
					MetaLoc& mb = iterb->second[b];
					if( (ma.ispos() && mb.ispos() && ma.stop == mb.stop) || (!ma.ispos() && !mb.ispos() && ma.start == mb.start) )
					{
						bothexist = true;
						break;
					}
				}
				if(bothexist)
					++ exist_ab;
			}
		}
		double only_exist_a = total_cds_a - exist_ab;
		double only_exist_b = total_cds_b - exist_ab;
		double sn = exist_ab / total_cds_a * 100;
		double sp = exist_ab / total_cds_b * 100; 
		double ha = 2 * sn * sp / (sn + sp);
		cout << iter_g->first << endl;
//		cout << '\t' << "#pred" << '\t' << "#benc" << '\t' << "#true" << '\t' << "#false" 
//			<< '\t' << "#miss" << '\t' << "%sn" << '\t' << "%sp" << '\t' << "%ha" << endl;
		cout << '\t' << total_cds_b << '\t' << total_cds_a << '\t' << exist_ab << '\t' << only_exist_b 
			<< '\t' << only_exist_a << std::setprecision(4) << '\t' << sn << '\t' << sp << '\t' << ha << endl;
		cout << std::setprecision(8);
	}
}

void output_metalocs(map<string, vector<MetaLoc> >& metaloc, const string& outfile)
{
	ofstream fout(outfile.data());
	if(!fout) { cout << "can't open " << outfile << endl; exit(-1); }

	int read_num = 0, cds_num = 0;
	map<string, vector<MetaLoc> >::iterator iter = metaloc.begin();
	for(; iter != metaloc.end(); ++ iter)
	{
		++ read_num;
		fout << '>' << iter->first << endl;
		for(int i = 0; i < (int)iter->second.size(); ++ i)
		{
			++ cds_num;
			MetaLoc& m = iter->second[i];
			fout << std::setw(8) << m.start << std::setw(8) << m.stop << std::setw(4) << m.strand << endl;
//				<< std::setw(4) << m.phase << std::setw(3) << m.complt / 10 << m.complt % 10 << endl;
		}
	}
	fout.close();
	cout << '\t' << cds_num << " cds in " << read_num << " reads outputs in " << outfile << endl;
}

void exit_with_help()
{
	cout << "metalocs-operates [options]" << endl
		 << "Version 3.0    2012/08/09    arguments:\n"
		 << "This program implements operates for metalocs according to the STOP coordinates, as follows:\n"
		 << "   -a FILE       input metalocs set A, always be regarded as the benchmark set.\n"
		 << "   -b FILE       input metalocs set B.\n"
		 << "   -c            accuracy comparison option.\n"
		 << "   -g FILE       assess accuracies on each group(usually genome), this file format is the same as binning file.\n"
		 << "   -o FILE       output file name.\n"
		 << "   --ab          intersection of set A and B, AB.\n"
		 << "   --a+b         union of set A and B, A+B.\n"
		 << "   --a/b         set A minus set B, A-B.\n"
		 << "   --b/a         set B minus set A, B-A.\n"
		 << "  comparison output format:\n"
		 << "  #pred\t#bench\t#true\t#false\t#miss\t%sn\t%sp\t%ha\n";
	exit(-1);
}

int main(int argc, char* argv[])
{
	// default settings
	string metaloc_file_a, metaloc_file_b, output_file;
	bool is_comparison = false, is_intersection= false, is_union = false, is_a_minus_b = false, is_b_minus_a = false;

	bool is_bin = false;
	string bin_file;
	
	// parse command line
	if(argc == 1) {
		exit_with_help();
		return -1;
	}

	pair<bool, string> opt;
	opt = command_line_parser(argc, argv, "-a", "", true);
	if(!opt.first) { cout << "MetaLocs file A is not defined, please check it." << endl; return -1; }
	metaloc_file_a = opt.second;
	opt = command_line_parser(argc, argv, "-b", "", true);
	if(!opt.first) { cout << "MetaLocs file B is not defined, please check it." << endl; return -1; }
	metaloc_file_b = opt.second;
	opt = command_line_parser(argc, argv, "-o", "", true);
	if(opt.first) output_file = opt.second;

	opt = command_line_parser(argc, argv, "-c", "", false);
	if(opt.first) is_comparison = true;
	opt = command_line_parser(argc, argv, "", "--ab", false);
	if(opt.first) is_intersection = true;
	opt = command_line_parser(argc, argv, "", "--a+b", false);
	if(opt.first) is_union = true;
	opt = command_line_parser(argc, argv, "", "--a/b", false);
	if(opt.first) is_a_minus_b = true;
	opt = command_line_parser(argc, argv, "", "--b/a", false);
	if(opt.first) is_b_minus_a = true;

	opt = command_line_parser(argc, argv, "-g", "", true);
	if(opt.first) {
		is_bin = true;
		bin_file = opt.second;
	}

	map<string, vector<MetaLoc> > metalocs_a = get_metalocs(metaloc_file_a);
	map<string, vector<MetaLoc> > metalocs_b = get_metalocs(metaloc_file_b);

	// Operate
	if(is_comparison) {
		if(!is_bin)
			comparison(metalocs_a, metalocs_b);
		else
			comparison(metalocs_a, metalocs_b, bin_file);
		return 0;
	} else if(is_union) {
		cout << "A plus B:" << endl;
		map<string, vector<MetaLoc> > metalocs_apb = union_metaloc(metalocs_a, metalocs_b);
		output_metalocs(metalocs_apb, output_file);
		return 0;
	} else if(is_intersection) {
		cout << "A inter B:" << endl;
		map<string, vector<MetaLoc> > metalocs_ab = inter_metaloc(metalocs_a, metalocs_b);
		output_metalocs(metalocs_ab, output_file);
		return 0;
	} else if(is_a_minus_b) {
		cout << "A minus B: " << endl;
		map<string, vector<MetaLoc> > metalocs_a_b = differ_metaloc(metalocs_a, metalocs_b);
		output_metalocs(metalocs_a_b, output_file);
		return 0;
	} else if(is_b_minus_a) {
		cout << "B minus A: " << endl;
		map<string, vector<MetaLoc> > metalocs_b_a = differ_metaloc(metalocs_b, metalocs_a);
		output_metalocs(metalocs_b_a, output_file);
		return 0;
	}
}
