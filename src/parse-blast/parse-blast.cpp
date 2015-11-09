
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

class Blast
{
public:
	string query, hit;
	double identity, e_value, bit_score;
	int match_bps, mismatch_bps, gap, query_start, query_stop, hit_start, hit_stop;
};

vector<string> parse_string(const string& s, const string& key = " ")
{
	vector<string> after_parse;
	int lpos = 0;
	int rpos = (int)s.find(key);
	while(rpos != std::string::npos)
	{
		string w = s.substr(lpos, rpos - lpos);
		lpos = rpos + (int)key.size();
		rpos = (int)s.find(key, lpos);
		after_parse.push_back(w);
	}
	string w = s.substr(lpos, rpos - lpos);
	after_parse.push_back(w);
	return after_parse;
}

map<string, vector<Blast> > parse_blast_output(const string& blast_file)
{
	ifstream fin(blast_file.data());
	if(!fin) { cout << "can't open " << blast_file << endl; exit(-1); }

	int hit_num = 0;
	map<string, vector<Blast> > blast_hit_map;
	while(!fin.eof())
	{
		Blast b;
		fin >> b.query >> b.hit >> b.identity >> b.match_bps >> b.mismatch_bps >> b.gap
			>> b.query_start >> b.query_stop >> b.hit_start >> b.hit_stop >> b.e_value >> b.bit_score;
		string line;
		getline(fin, line);
		if(b.query.empty())
			continue;
		blast_hit_map[b.query].push_back(b);
		++ hit_num;
	}
	fin.close();
	cout << hit_num << " hits for " << blast_hit_map.size() << " query genes." << endl; 
	return blast_hit_map;
}

map<string, vector<string> > parse_blast(map<string, vector<Blast> >& blast_hit_map, double e_value, double bit_score, double identity)
{
	map<string, vector<string> > satisfied_queries;
	int satisfied_num = 0;

	map<string, vector<Blast> >::iterator iterquery = blast_hit_map.begin();
	for(; iterquery != blast_hit_map.end(); ++ iterquery)
	{
		for(int ii = 0; ii < (int)iterquery->second.size(); ++ ii)
		{
			Blast& b = iterquery->second[ii];
			if(b.e_value <= e_value && b.bit_score >= bit_score && b.identity >= identity)
			{
				satisfied_queries[parse_string(iterquery->first, "|")[0]].push_back(iterquery->first);
				++ satisfied_num;
				break;
			}
		}
	}
	cout << satisfied_num << " queries satisfied the given threshold." << endl;
	return satisfied_queries;
}

void output(map<string, vector<string> >& satisfied_queries, const string& file)
{
	ofstream fout(file.data());
	if(!fout) { cout << "can't open " << file << endl; exit(-1); }

	map<string, vector<string> >::iterator iterr = satisfied_queries.begin();
	for(; iterr != satisfied_queries.end(); ++ iterr)
	{
		fout << '>' << iterr->first << endl;
		vector<string>::iterator iterg = iterr->second.begin();
		for(; iterg != iterr->second.end(); ++ iterg)
		{
			vector<string> vv = parse_string(*iterg, "|");
			fout << setw(8) << vv[1] << setw(8) << vv[2] << setw(4) << vv[3] << endl;
		}
	}
	fout.close();
}

void help()
{
	cout<< "parse-blast  version 1.0" << endl << endl
		<< "    -f    blast result file." << endl
		<< "    -e    threshold of e value, default, 1." << endl
		<< "    -b    threshold of bit score, default, 0." << endl
		<< "    -i    threshold of identity, default, 0." << endl
		<< "    -o    output file, default, out.txt." << endl;
}

int main(int argc, char* argv[])
{
	if(argc == 1) 
	{
		help();
		return -1;
	}

	double e_value = 1, bit_score = 0, identity = 0;
	string blast_file, out_file = "out.txt";

	int ii = 1;
	for(; ii < argc; ++ ii)
	{
		if(argv[ii][0] != '-') break;
		string option(argv[ii]);
		if(option == "-f")
			blast_file = string(argv[ii+1]);
		else if(option == "-e")
			e_value = atof(argv[ii+1]);
		else if(option == "-b")
			bit_score = atof(argv[ii+1]);
		else if(option == "-i")
			identity = atof(argv[ii+1]);
		else if(option == "-o")
			out_file = string(argv[ii+1]);
		else 
		{
			cout << "No " << option << " option." << endl;
			return -1;
		}
		++ ii;
	}

	if(blast_file.empty())
	{
		cout << "The BLAST file is not defined." << endl;
		return -1;
	}

	map<string, vector<Blast> > blast_results = parse_blast_output(blast_file);
	map<string, vector<string> > satisfied_queries = parse_blast(blast_results, e_value, bit_score, identity);
	output(satisfied_queries, out_file);
}
