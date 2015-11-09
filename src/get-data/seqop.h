
// author: Liu Yong-Chu
// version 2.0
// history: created on 2009/11/05
//			revised on 2010/04/03:	split into 2 files, seqop.h and seqop.cpp

#ifndef SEQOP_H
#define SEQOP_H

#include <string>
#include <vector>
#include <assert.h>
#include <math.h>
using namespace std;

// define the sequence operate

namespace SeqOp
{
	char Letter2Digital(const char res);		// convert letter format of DNA sequence to digital format
	std::string Letter2Digital(const string& seq);

	char Digital2Letter (const char res);	// convert digital format of DNA sequence to letter format
	string Digital2Letter(const string& seq);

	char revs(const char res);				// transfer to opposite strand
	string revs(const string& s);

	int motif2Digital(const string& seq);	// transfermation between motif and digital
	string digital2Motif(const int index, const int len);

	string drawseq(const string& genome,int start,unsigned length);	// extract sub-sequence 
	string drawORF(const string& genome, int left, int right, char strand);	// extract an ORF

	string toAminoSeq(const char* seq, int startPos, int endPosition);	// convert sequence to amino acid sequence

	vector<double> aminoEDP(const char* seq, int startPos, int endPosition);		// EDP for amino acid (20)
	vector<double> codonEDP(const char* seq, int startPos, int endPosition);		// EDP for codon (61)

	double gccontent(const string& seq);		// gc content
	double theta(const string& seq);			// theta = (gc1+gc3)/(gc1+gc2+gc3)
};

#endif
