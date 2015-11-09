
// Methods for operating genome sequence

#ifndef OPERATE_H
#define OPERATE_H

#include <string>
#include "Trans.h"
using namespace std;

string ReadFNA(const string& fna);

string drawseq(const string& genome,int start,unsigned length);

string drawORF(const string& genome, int left, int right, char strand);

// Entropy Density Profile of Amino (20)
vector<double> AAEDP(const char* seq, int startPos, int endPosition);

// Entropy Density Profile of Condon (61)
vector<double> CodonEDP(const char* seq, int startPos, int endPosition);

// DNA to Amino Acid
//string toAminoSeq(const char* seq, int startPos, int endPosition);

#endif
