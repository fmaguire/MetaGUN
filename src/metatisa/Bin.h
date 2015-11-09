#ifndef BIN_H
#define BIN_H

#include <time.h>
#include "Trans.h"
#include "Operate.h"
#include "TypeDefBase.h"

using namespace std;

static const int unitC = 100;
static int motif_len = 8;
static bool binaryFormat = true;

string distri_file = "bin_distri.txt";

// learning the frequency of the m-mer oligonucleotide for a gnome given the length m
std::map<int, float> motif_frequency(const string& genomeSeq, const int len)
{
	int genome_size= (int)genomeSeq.size();
	std::map<int, float> freqs;
	int i=0;
	for(; i<genome_size-len+1; ++i)
		++ freqs[ motif2Digital(genomeSeq.substr(i,len)) ];
	return freqs;
}

void outputFreqs2file(std::map<int, float>& freqs, std::ofstream& fout)
{
	std::map<int, float>::iterator iter = freqs.begin();
	int count = 0;
	for(; iter != freqs.end(); ++ iter) {
		if(iter->second == 0)
			continue;
		if(binaryFormat)
		{
			fout.write((char*)(&(iter->first)), sizeof(iter->first));
			fout.write((char*)(&(iter->second)), sizeof(iter->second));
		}
		else 
		{
			if( count % unitC == 0 )
				fout << "\n\t";
			fout << iter->first << ':' << iter->second << ' ';
			++ count;
		}
	}
	if(!binaryFormat)	fout << endl;
}

void get_Freqs_table(const string& file,
					 std::vector<string>& genomes,
					 std::vector<std::vector<float> >& table)
{
	genomes.erase(genomes.begin(), genomes.end());
	table.erase(table.begin(), table.end());

	std::ifstream fin;
	if(binaryFormat)
		fin.open(file.data(), ios::binary);
	else
		fin.open(file.data());
	if(!fin) { cout << "can't open " << file << endl; exit(0); }

	int genomeC = 0;
	if(binaryFormat) {
		fin.read((char*)(&genomeC), sizeof(genomeC));
		fin.read((char*)(&motif_len), sizeof(motif_len));
	}
	else 
		fin >> genomeC >> motif_len;
	int motifSize = (int)pow(4.0, motif_len);

	if(binaryFormat)
	{
		for(int k = 0; k < genomeC; ++k) {
			string genome;
			int freqSize = 0;
			float sum = 0.0f;
			char g[20];
			fin.read((char*)(&g), sizeof(g));
			fin.read((char*)(&freqSize), sizeof(freqSize));
			genome = string(g);
			genomes.push_back(genome);
			vector<float> freqs(motifSize, 1e-3f);
			int f = 0;
			for(f = 0; f < freqSize; ++f) {
				int index = 0;
				float p = 0.0f;
				fin.read((char*)(&index), sizeof(index));
				fin.read((char*)(&p), sizeof(p));
				freqs[index] = p;
				sum += p;
			}
			for(f = 0; f < (int)freqs.size(); ++f)
				freqs[f] /= sum;
			table.push_back(freqs);
		}
		fin.close();
		cout << "bin-model of " << genomes.size() << " genomes read in" << endl;
		return;
	}

	string line;
	getline(fin, line);
	while(!fin.eof()) {
		if(line.find('>') != std::string::npos) {
			int lpos = (int) line.find('>');
			int rpos = (int) line.find('\t');
			genomes.push_back( line.substr(lpos + 1, rpos - lpos) );
			int count = atoi(line.substr(rpos + 1).data());
			string word;
			std::vector<float> freqs(motifSize, 1e-5f);
			float sum = 0;
			for(int c = 0; c < count; ++c) {
				fin >> word;
				lpos = (int) word.find(':');
				int key = atoi(word.substr(0, lpos).data());
				float value = (float)atof(word.substr(lpos+1).data());
				sum += value;
				freqs[ key ] = value;
			}
			for(int k = 0; k < (int)freqs.size(); ++ k)
				freqs[k] /= sum;
			table.push_back( freqs );
			while(line.find('>') != std::string::npos && !fin.eof())
				getline(fin, line);
		}
		else
			getline(fin, line);
	}
	fin.close();
	cout << "bin-model of " << genomes.size() << " genomes read in" << endl;
}

// caculate the probabilities of a given sequence to belong to a series of genomes
std::vector<float> given_seq_prob(const string& seq, 
								   std::vector<std::vector<float> >& table,
								   const int motifLen)
{
	int no_genome = (int)table.size();
	std::vector<float> probs(no_genome,1.0);	

	int i=0;
	for(i=0; i<(int)seq.size()-motifLen+1; ++i) {
		int key = motif2Digital(seq.substr(i, motifLen));
		int j = 0;
		// 为了防止溢出，每次归一化
		float sum = 0;
		for(j=0; j<no_genome; ++j) {
			probs[j] *= table[j][key];
			sum += probs[j];
		}
		for(j=0; j<no_genome; ++j)
			probs[j] /= sum;
	}
	return probs;
}

void get_group_map(const std::string& file, std::vector<string>& genomes, std::vector<string>& groups)
{
	std::ifstream fin(file.data());
	if(!fin) {cout<<"can't read "<<file<<endl; exit(0);}
	genomes.erase(genomes.begin(), genomes.end());
	groups.erase(groups.begin(), groups.end());
	while(!fin.eof()) {
		string g, k;
		fin >> g >> k;
		if( !g.empty() && !k.empty())
			genomes.push_back(g), groups.push_back(k);
	}
	fin.close();
	cout << "genomes : " << genomes.size() << endl;

	// sort
	std::vector<string> g_sort = genomes;
	std::sort(g_sort.begin(), g_sort.end());
	std::vector<string> k_sort(groups.size());
	for(int i=0; i<(int)g_sort.size(); ++i) {
		int j = 0;
		for(; j<(int)genomes.size(); ++j)
			if(g_sort[i] == genomes[j])
				break;
		k_sort[i] = groups[j];
	}
	genomes = g_sort;
	groups = k_sort;
}

std::multimap<string, string> bin_fasta(std::vector<std::vector<float> >& table,
										map<string, pair<string, string> >& id_seq,
										const std::vector<string>& genomes,
										const std::vector<string>& groups)
{
	std::set<string> group_types;
	for(int k=0; k<(int)groups.size(); ++k)
		group_types.insert(groups[k]);
	cout << group_types.size() << " group types" << endl;
	
	int count = 0;
	std::multimap<string, string> GenusIdMap;
	std::map<string, float> distri;
	time_t t = time(NULL);

	map<string, pair<string, string> >::iterator iter = id_seq.begin();
	for(; iter != id_seq.end(); ++ iter)
	{
		if(++count % 10000==0)
			cout << count << '(' << time(NULL) -t << ')' << " --> ";
		std::vector<float> probs = given_seq_prob(iter->second.first, table, motif_len);
		int max_i = 0;
		float max_v = 0;
		for(int k=0; k<(int)probs.size(); ++k) {
			if(std::max(max_v, probs[k]) != max_v) {
				max_i = k, max_v = probs[k]; }
		}
		string bin_k = groups[max_i];
		GenusIdMap.insert( std::make_pair(bin_k, iter->first) );

		std::set<string>::iterator iter = group_types.begin();
		for(; iter != group_types.end(); ++iter) {
			if(bin_k == *iter) {
				++ distri[ bin_k ];
			}
		}
	}

	cout << count << " shotgun sequence" << endl;
	
	std::ofstream fout;
	fout.open(distri_file.data());
	if(!fout) { cout << "can't open bin_distribution.txt" << endl; exit(0); }

	std::map<string, float>::iterator it = distri.begin();
	for(; it != distri.end(); ++it)
		fout << std::setw(30) << it->first 
		<< std::setw(15) << it->second << endl;

	return GenusIdMap;
}

#endif
