
#ifndef BIN_H
#define BIN_H

#include "seqop.h"
#include <map>
#include <fstream>
#include <algorithm>
#include <iomanip>

// learning the frequency of the k-mer oligonucleotide
std::map<int, float> get_kmer_freqs(const string& genomeSeq, const int kmer_len)
{
	std::map<int, float> freqs;
	float kmer_sum = 0;

	// positive strand
	int i = 0;
	for(i = 0; i < (int)genomeSeq.size() - kmer_len + 1; ++i) 
	{
		++ freqs[ SeqOp::motif2Digital(genomeSeq.substr(i, kmer_len)) ];
		++ kmer_sum;
	}
	
	// normalize
	map<int, float>::iterator iterf = freqs.begin();
	for(; iterf != freqs.end(); ++ iterf)
		iterf->second /= kmer_sum;

	return freqs;
}

// train bin-model 
map<string, map<int, float> > train_binmodel(const string& fnalist_file, const int kmer_len)
{
	// read in the genome list
	ifstream fin(fnalist_file.data());
	if(!fin) { cout << "can't open " << fnalist_file << endl; exit(0); }

	vector<string> fnafiles;
	while(!fin.eof())
	{
		string line;
		getline(fin, line);
		if(!line.empty())
			fnafiles.push_back(line);
	}
	fin.close();

	std::sort(fnafiles.begin(), fnafiles.end());
	cout << "Modeling " << kmer_len << "-mer frequency from " << fnafiles.size() << " genomes" << endl;

	map<string, map<int, float> > binmodel;
	for(int i = 0; i < (int)fnafiles.size(); ++ i)
	{
		string fnapath = fnafiles[i];
		// for windows
		string refseqno = (fnapath.rfind('\\') == std::string::npos ? fnapath : fnapath.substr(fnapath.rfind('\\')+1));
		// for unix
		//string refseqno = (fnapath.rfind('/') == std::string::npos ? fnapath : fnapath.substr(fnapath.rfind('/')+1));
		refseqno = refseqno.substr(0, refseqno.find('.'));

		cout << fnapath << endl;
		map<int, float> freqs = get_kmer_freqs(SeqOp::readfna(fnapath), kmer_len);
		binmodel[ refseqno ] = freqs;
	}
	return binmodel;
}

// output binmodel, 
void output_binmodel(map<string, map<int, float> >& binmodel, const int kmer_len,
					 const string& file, const bool isbinary = true)
{
	if(isbinary)
	{
		ofstream fout(file.data(), ios::binary);
		if(!fout) { cout << "can't open " << file << endl; exit(0); }

		int genome_num = (int)binmodel.size();
		fout.write((char*)(&genome_num), sizeof(genome_num));
		fout.write((char*)(&kmer_len), sizeof(kmer_len));

		map<string, map<int, float> >::iterator iterg = binmodel.begin();
		for(; iterg != binmodel.end(); ++ iterg)
		{
			string g = iterg->first;
			char space = ' ';
			for(int ii = 0; ii < 20; ++ii)
				ii < (int)g.size() ? fout.write((char*)(&(g[ii])), sizeof(g[ii])) : fout.write((char*)(&space), sizeof(space));
			int kmer_num = (int)iterg->second.size();
			fout.write((char*)(&kmer_num), sizeof(kmer_num));
			map<int, float>::iterator iterf = iterg->second.begin();
			for(; iterf != iterg->second.end(); ++ iterf)
			{
				fout.write((char*)(&iterf->first), sizeof(iterf->first));
				fout.write((char*)(&iterf->second), sizeof(iterf->second));
			}
		}
		fout.close();
		return;
	}

	ofstream fout(file.data());
	if(!fout) { cout << "can't open " << file << endl; exit(0); }
	int genome_num = (int)binmodel.size();
	fout << '\t' << genome_num << '\t' << kmer_len << endl;

	map<string, map<int, float> >::iterator iterg = binmodel.begin();
	for(; iterg != binmodel.end(); ++ iterg)
	{
		fout << '>' << iterg->first << '\t' << iterg->second.size();
		int count = 0;
		for(map<int, float>::iterator iterf = iterg->second.begin(); iterf != iterg->second.end(); ++ iterf)
		{
			if(count % 100 == 0)
				fout << "\n\t";
			fout << iterf->first << ':' << iterf->second << ' ';
			++ count;
		}
		fout << endl;
	}
	fout.close();
}

// read binmodel from file
map<string, vector<float> > get_binmodel(const string& binmodel_file, int& kmer_len, const bool isbinary = true)
{
	cout << "get bin-model from " << binmodel_file << " ..." << endl;
	map<string, vector<float> > binmodel;
	if(isbinary)
	{
		ifstream fin(binmodel_file.data(), ios::binary);
		int genome_num = 0;
		fin.read((char*)(&genome_num), sizeof(genome_num));
		fin.read((char*)(&kmer_len), sizeof(kmer_len));
		for(int i = 0; i < genome_num; ++ i)
		{
			string g;
			for(int ii = 0; ii < 20; ++ ii)
			{
				char ch = ' ';
				fin.read((char*)(&ch), sizeof(ch));
				if(ch != ' ')
					g += ch;
			}
			int kmer_num = 0;
			fin.read((char*)(&kmer_num), sizeof(kmer_num));

			int kmer_size = (int)pow(4.0, kmer_len);
			float pseudocount = (float)(1/kmer_size)*(1e-3f);
			vector<float> freqs(kmer_size, pseudocount);
			for(int k = 0; k < kmer_num; ++ k)
			{
				int index = 0;
				float f = 0.0f;
				fin.read((char*)(&index), sizeof(index));
				fin.read((char*)(&f), sizeof(f));
				freqs[index] = f;
			}
			binmodel[ g ] = freqs;
		}
	}

	else
	{
		ifstream fin(binmodel_file.data());
		if(!fin) { cout << "can't open " << binmodel_file << endl; exit(0); }

		int genome_num = 0;
		fin >> genome_num >> kmer_len;
		int kmer_size = (int)pow(4.0, kmer_len);

		string line;
		getline(fin, line);
		while(!fin.eof())
		{
			if(line[0] == '>') 
			{
				int lpos = (int) line.find('>');
				int rpos = (int) line.find('\t');
				string genome = line.substr(lpos + 1, rpos - lpos - 1);
				int count = atoi(line.substr(rpos + 1).data());
				string word;
				float pseudocount = (float)(1/kmer_size)*(1e-3f);
				std::vector<float> freqs(kmer_size, pseudocount);
				for(int c = 0; c < count; ++c) 
				{
					fin >> word;
					lpos = (int) word.find(':');
					int key = atoi(word.substr(0, lpos).data());
					float value = (float)atof(word.substr(lpos+1).data());
					freqs[ key ] = value;
				}
				binmodel[ genome ] = freqs;
				while(line.find('>') != std::string::npos && !fin.eof())
					getline(fin, line);
			}
			else
				getline(fin, line);
		}
		fin.close();
	}
	return binmodel;
}

map<string, pair<string, string> > get_binmap(const string& binmap_file)
{
	cout << "get bin-map from " << binmap_file << " ..." << endl;
	ifstream fin(binmap_file.data());
	if(!fin) { cout << "can't open " << binmap_file << endl; exit(0); }

	map<string, pair<string, string> > binmap;
	while(!fin.eof())
	{
		string refseqno, genus, group;
		string line;
		fin >> refseqno >> genus >> group;
		getline(fin, line);
		if(!refseqno.empty())
			binmap.insert(std::make_pair(refseqno, std::make_pair(genus, group)));
	}
	return binmap;
}

double lnprob4seq(vector<float>& freqs, const string& seq, const int kmer_len)
{
	double lnprob = 0;
	for(int i = 0; i < (int)seq.size() - kmer_len + 1; ++ i)
	{
		int index = SeqOp::motif2Digital(seq.substr(i, kmer_len));
		lnprob += (double)log( freqs[index] );
	}
	return lnprob;
}

map<string, string > bin_seqs(map<string, vector<float> >& binmodel, map<string, string>& seqs, const int kmer_len)
{
	cout << "Sequences binning ... " << endl << '\t';

	map<string, string> read_bing;
	int count = 0;
	map<string, string>::iterator iterread = seqs.begin();
	for(; iterread != seqs.end(); ++ iterread)
	{
		if((count++) % 1000 == 0)
			cout << "..*";
		string& seq = iterread->second;
		double maxlnprob = -1e100;
		string bing;
		map<string, vector<float> >::iterator iterf = binmodel.begin();
		for(; iterf != binmodel.end(); ++ iterf)
		{
			double lnprob = lnprob4seq(iterf->second, seq, kmer_len);
			if(lnprob > maxlnprob)
			{
				maxlnprob = lnprob;
				bing = iterf->first;
			}
		}
		read_bing[iterread->first] = bing;
	}
	cout << endl;

	return read_bing;
}

map<string, double> get_binning_profile(map<string, string>& read_bing, map<string, vector<float> >& binmodel)
{
	map<string, double> bin_profile;

	map<string, vector<float> >::iterator iterbm = binmodel.begin();
	for(; iterbm != binmodel.end(); ++iterbm)
		bin_profile[iterbm->first] = 1e-100;
	double total = 1e-100;

	map<string, string>::iterator iterr = read_bing.begin();
	for(; iterr != read_bing.end(); ++ iterr)
	{
		++ total;
		++ bin_profile[iterr->second];
	}
	map<string, double>::iterator iterb = bin_profile.begin();
	for(; iterb != bin_profile.end(); ++ iterb)
		iterb->second /= total;
	return bin_profile;
}

double binning_profile_diff(map<string, double>& bp1, map<string, double>& bp2)
{
	double diff = 0;
	map<string, double>::iterator iterbp1 = bp1.begin(), iterbp2 = bp2.begin();
	for(; iterbp1 != bp1.end(), iterbp2 != bp2.end(); ++ iterbp1, ++ iterbp2)
		diff += fabs(iterbp1->second - iterbp2->second);
	return diff;
}

map<string, string> bin_seqs_prior(map<string, vector<float> >& binmodel, map<string, string>& seqs, const int kmer_len,
									 map<string, string>& read_bing)
{
	// get new binning profiles
	map<string, double> bin_profile = get_binning_profile(read_bing, binmodel);
	map<string, string> this_read_bing;

	// binning based on current binning profile
	int count = 0;
	map<string, string>::iterator iterread = seqs.begin();
	for(; iterread != seqs.end(); ++ iterread)
	{
		if((count++) % 1000 == 0)
			cout << "..*";
		string& seq = iterread->second;
		double maxlnprob = -1e100;
		string bing;
		map<string, vector<float> >::iterator iterf = binmodel.begin();
		for(; iterf != binmodel.end(); ++ iterf)
		{
			double lnprob = lnprob4seq(iterf->second, seq, kmer_len);

			map<string, double>::iterator iterbp = bin_profile.find(iterf->first);
			double ln_prior_gp = 0;
			if(iterbp == bin_profile.end())
				ln_prior_gp = -50;
			else 
				ln_prior_gp = log(iterbp->second);
			lnprob += ln_prior_gp;

			if(lnprob > maxlnprob)
			{
				maxlnprob = lnprob;
				bing = iterf->first;
			}
		}
		this_read_bing[iterread->first] = bing;
	}
	cout << endl;

	map<string, double> this_bin_profile = get_binning_profile(this_read_bing, binmodel);
	double profile_diff = binning_profile_diff(this_bin_profile, bin_profile);
	cout << "profile diff: " << profile_diff << endl;

	return this_read_bing;
}

// output bin results
void output_bin(map<string, string>& read_bing, map<string, pair<string, string> >& binmap,
				const string& bin_file, const string& bindistri_file)
{
	ofstream foutb(bin_file.data());
	if(!foutb) { cout << "can't open " << bin_file << endl; exit(-1); }

	map<string, double> bindistri_genome;
	map<string, double> bindistri_genus;
	map<string, double> bindistri_group;
	double read_total = 1e-100;

	map<string, string>::iterator iterread = read_bing.begin();
	for(; iterread != read_bing.end(); ++ iterread)
	{
		string title = iterread->first;
		string bin_genome = iterread->second;
		map<string, pair<string, string> >::iterator iterbm = binmap.find(bin_genome);
		if(iterbm == binmap.end())
		{
			cout << "can't find bin-map item for " << bin_genome << endl;
			continue;
		}
		string bin_genus = iterbm->second.first, bin_group = iterbm->second.second;
		foutb << '>' << title << "\t" << bin_genus << "\t" << bin_group << endl;
		++ read_total;
		++ bindistri_genome[bin_genome];
		++ bindistri_genus[bin_genus];
		++ bindistri_group[bin_group];
	}
	foutb.close();

	ofstream foutd(bindistri_file.data());
	if(!foutd) { cout << "can't open " << bindistri_file << endl; exit(-1); }

	map<string, double>::iterator iterd = bindistri_group.begin();
	foutd << std::setprecision(4) << ">group_distribution" << endl;
	for(iterd = bindistri_group.begin(); iterd != bindistri_group.end(); ++ iterd)
		foutd << setw(30) << iterd->first << setw(15) << iterd->second / read_total * 100 << setw(15) << iterd->second << endl;

	foutd << endl << ">genus_distribution" << endl;
	for(iterd = bindistri_genus.begin(); iterd != bindistri_genus.end(); ++ iterd)
		foutd << setw(30) << iterd->first << setw(15) << iterd->second / read_total * 100 << setw(15) << iterd->second << endl;

	foutd << endl << ">genome_distribution" << endl;
	for(iterd = bindistri_genome.begin(); iterd != bindistri_genome.end(); ++ iterd)
		foutd << setw(30) << iterd->first << setw(15) << iterd->second / read_total * 100 << setw(15) << iterd->second << endl;
	foutd.close();
}
#endif

