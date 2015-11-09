
// Define often used functions
// Last update: 2012-11-10
#ifndef OFTENUSEDFUN_H
#define OFTENUSEDFUN_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
using namespace std;

vector<string> parse_string(const string& s, const string& key = " ");
vector<string> split_string(const string& s);

string num2str(const int n);
string num2str(const int n, const int len);

vector<string> readlines(const string& file);
map<string, string> read_fasta(const string& file);

pair<bool, string> command_line_parser(int argc, char* argv[], const string& short_opt, 
									   const string& long_opt, const bool has_argv);

#endif 
