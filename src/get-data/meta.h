
#ifndef META_H
#define META_H

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include "seqop.h"
#include "feature.h"
#include "matrix.h"
#include "tis.h"

class Codons
{
public:
	set<int> starts, stops;
	void set_ss_codons()
	{
		starts.erase(starts.begin(), starts.end());
		stops.erase(stops.begin(), stops.end());
		// ATG(14), CTG(30), GTG(46), TTG(62)
		starts.insert(14), starts.insert(46), starts.insert(62);
		// TGA(48), TAG(50), TAA(56)
		stops.insert(48), stops.insert(50), stops.insert(56);
	}
};

class Locs
{
public:
	Locs() { begin = end = -1; strand = '+'; }
	int begin, end;
	char strand;
	int ispos() { if(strand == '+') return true; else return false; }
	int start() { if(strand == '+') return begin; else return end; }
	int stop() { if(strand == '+') return end; else return begin; }
	int size() { return end - begin + 1; }
	void adjust_start(const int new_pos) { if(strand == '+') begin = new_pos; else end = new_pos; }
};

class MetaORF : public Locs
{
public:
	string coding_type, cut_type;
	bool iscoding() { if(coding_type == "C") return true; else return false; }
	bool is_leftmost, is_stopin;
	int leftmost_xtg, phase;

	// TIS
	bool has_tis_cand;
	double downfalse_ts, upfalse_ts, true_ts;
	int true_tis_site;

	MetaORF()
	{
		//Locs();
		leftmost_xtg = phase = -1;
		is_leftmost = is_stopin = false;
		cut_type = "--";
		has_tis_cand = false;
		downfalse_ts = upfalse_ts = true_ts = 1e-10;
		true_tis_site = -1;
	}
};

class Read
{
public:
	Read() {}
	string seq, rid, descript;
	map<string, vector<MetaORF> > ORFs;
	map<string, vector<MetaORF> > ORFsets;
};

// to find all longest ORFs on positive strand
vector<MetaORF> find_pos_LORFs(const string& seq, const int add_size, const int step,
							   bool is_pos_strand, const int min_ORF_len = 60);

// to find all the longest ORFs on both positive and negative strands
vector<MetaORF> find_LORFs(const string& seq, const int add_size, const int step, const int min_ORF_len = 60);

// get meta ORFs from file
map<string, vector<MetaORF> > get_metaORF(const string& locfile);

class XTGScore
{
public:
	double true_score, upfalse_score, downfalse_score;
	XTGScore() { true_score = upfalse_score = downfalse_score = 0; }
};

// scoring a XTG, which is a TIS candidate, using a TriTISA TIS model
XTGScore scoring_xtg(const string& seq, int xtg_site, const TIS& tis);

class MetaSeqs
{
public:
	map<string, Read> metaseqs;
	map<string, pair<string, string> > readbin;
	map<string, vector<string> > genus_binreads;
	map<string, vector<string> > group_binreads;

	int all_ORF_num, cod_ORF_num, noncod_ORF_num, reloc_ORF_num;
	map<string, int> cut_type_num;

	int max_ORF_len, min_ORF_len;
	bool is_unix;

	void initiate_settings()
	{
		all_ORF_num = cod_ORF_num = noncod_ORF_num = reloc_ORF_num = 0;
		cut_type_num["00"] = cut_type_num["01"] = cut_type_num["10"] = cut_type_num["11"] = 0;
		min_ORF_len = 60, max_ORF_len = 1200;
		is_unix = false;
	}

	MetaSeqs() { initiate_settings(); }

	void get_metaseqs(const string& seqs_file);
	void find_all_LORFs();
	void find_all_ORFsets();
	void label_ORFs(const string& coding_cds_file);
	void label_ORFs_shadow_rule(const string& coding_cds_file);

	void get_bin_info(const string& bin_file);
	void tis_scoring(const string& prior_para_path);
	void tis_scoring_ORFsets(const string& prior_para_path);

	void get_allfea(const string& name, bool is_ORFsets);
	void get_allfea(const string& name, int level, bool is_ORFsets);
	void get_trainfea(const string& name);
	void get_trainfea(const string& name, int level);
	void get_trainfea_group(const std::string &name);

	void fout_ORF_seqs(const string& outfile, bool is_aa = false);
	void fout_cod_seqs(const string& outfile, bool is_amino = true);
	void fout_ORF_locs(const string& outfile, bool is_coding_only = false);
	void fout_shd_ORF_locs(const string& outfile_basename);
	int seqs_num() { return (int)metaseqs.size(); }
	void get_seqs_info(const string& name);
};

#endif
