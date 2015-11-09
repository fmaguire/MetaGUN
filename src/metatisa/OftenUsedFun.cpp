
#include "OftenUsedFun.h"

vector<string> parse_string(const string& s, const string& key)
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

// NOTE THAT std::string::npos = -1
vector<string> split_string(const string& s)
{
	vector<string> after_split;
	int lpos = 0, rpos = 0;
	int rpos1 = (int)s.find(' ');
	int rpos2 = (int)s.find('\t');
	if(rpos1 == -1 && rpos2 == -1) rpos = -1;
	else if(rpos1 == -1 && rpos2 != -1) rpos = rpos2;
	else if(rpos1 != -1 && rpos2 == -1) rpos = rpos1;
	else rpos = std::min(rpos1, rpos2);

	while(rpos != std::string::npos)
	{
		if(lpos != rpos)
			after_split.push_back(s.substr(lpos, rpos - lpos));
		lpos = rpos + 1;
		rpos1 = (int)s.find(' ', lpos);
		rpos2 = (int)s.find('\t', lpos);
		if(rpos1 == -1 && rpos2 == -1) rpos = -1;
		else if(rpos1 == -1 && rpos2 != -1) rpos = rpos2;
		else if(rpos1 != -1 && rpos2 == -1) rpos = rpos1;
		else rpos = std::min(rpos1, rpos2);
	}
	after_split.push_back(s.substr(lpos, rpos - lpos));
	return after_split;
}

string num2str(const int n)
{
	if(n==0) return "0";	
	int nn = (n>0 ? n : -n);
	string s;
	while(nn) {
		char nc = (nn%10) + 48;
		s = nc + s;
		nn /= 10;
	}
	if(n<0) s = "-" + s;
	return s;
}

string num2str(const int n, const int len)
{
	string s = num2str(n);
	if(len < (int)s.size()) {
		cout << "given length is not enough." << endl;
		return s;
	}

	string ss(len, '0');
	int i = 0;
	while(++i <= (int)s.size())
		ss[(int)ss.size() - i] = s[(int)s.size() -i];
	return ss;
}

vector<string> readlines(const string& file)
{
	ifstream fin(file.data());
	if(!fin) { cout << "can't open " << file << endl; exit(-1); }

	vector<string> lines;
	string line;
	while(!fin.eof())
	{
		line.clear();
		getline(fin, line);
		if(line.empty())
			continue;
		lines.push_back(line);
	}
	fin.close();
	return lines;
}

map<string, string> read_fasta(const string& file)
{
	cout << "Get input sequences ... " << endl;
	ifstream fin(file.data());
	if(!fin) { cout << "can't open " << file << endl; exit(-1); }

	map<string, string> seqs;
	string line;
	getline(fin, line);
	while(!fin.eof())
	{
		if(line[0] == '>')
		{
			string title = line.substr(1);
			string seq;
			getline(fin, line);
			while(line[0] != '>' && !fin.eof())
			{
				seq += line;
				getline(fin, line);
			}
			seqs[title] = seq;
		}
		else
			getline(fin, line);
	}
	fin.close();
	cout << '\t' << seqs.size() << " input fasta sequences." << endl;
	return seqs;
}

pair<bool, string> command_line_parser(int argc, char* argv[], const string& short_opt, 
									   const string& long_opt, const bool has_argv)
{
	for(int i=0; i<argc; ++i) {
		if(string(argv[i])==short_opt || string(argv[i])==long_opt) {
			if(has_argv) {
				if(i+1 >= argc) {
					cout << "Option " << argv[i] << " is not specified." << endl;
					return std::make_pair(true, "");
				} else
					return std::make_pair(true, string(argv[i+1]));
			} else
				return std::make_pair(true, "");
		}
	}
	return std::make_pair(false, "");
}
