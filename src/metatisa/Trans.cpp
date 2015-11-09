
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Methods for format transformation
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include "Trans.h"

// convert letter format of DNA sequence to digital format
char Letter2Digital(const char res)
{
	switch(res)
	{
	case 'A' : case 'a' : case '0' : return '0';
	case 'C' : case 'c' : case '1' : return '1';
	case 'G' : case 'g' : case '2' : return '2';
	case 'T' : case 't' : case '3' : return '3';
	case 'N' : case 'n' : return '2';
	case 'X' : case 'x' : return '0';
	case 'H' : case 'h' : return '3';
	case 'M' : case 'm' : return '1';
	case 'K' : case 'k' : return '2';
	case 'D' : case 'd' : return '0';
	case 'R' : case 'r' : return '2';
	case 'Y' : case 'y' : return '3';
	case 'S' : case 's' : return '1';
	case 'W' : case 'w' : return '0';
	case 'B' : case 'b' : return '1';
	case 'V' : case 'v' : return '2';
	default  : 
		cout << "!UNEXPECTED nucleotide acid!" << endl;
		return '$';
	}
}

string Letter2Digital(const string& seq)
{
	string revSeq = seq;
	size_t i=0;
	for( i=0; i<seq.size(); ++i )
		revSeq[i] = Letter2Digital( seq[i] );
	return revSeq;
}

// convert digital format of DNA sequence to letter format
char Digital2Letter (const char res)
{
	switch(res)
	{
	case '0' : case 'A' : case 'a' : return 'A';
	case '1' : case 'C' : case 'c' : return 'C';
	case '2' : case 'G' : case 'g' : return 'G';
	case '3' : case 'T' : case 't' : return 'T';
	default  : 
		cout << "!UNEXPECTED nucleotide acid!" << endl;
		return '$';
	}
}

string Digital2Letter(const string& seq)
{
	string revSeq = seq;
	size_t i=0; 
	for( i=0; i<revSeq.size(); ++i )
		revSeq[i] = Digital2Letter( seq[i] );
	return revSeq;
}

// transfer to opposite strand
char revs(const char res)
{
	switch(res )
	{
	case 'A' : case 'a' : return 'T';
	case 'C' : case 'c' : return 'G';
	case 'G' : case 'g' : return 'C';
	case 'T' : case 't' : return 'A';
	case '0' : return '3';
	case '1' : return '2';
	case '2' : return '1';
	case '3' : return '0';
	default :
		cout << "!UNEXPECTED nucleotide acid!" << endl;
		return '$';
	}
}

string revs(const string& s)
{
	size_t i=0;
	string s1;
	s1.resize(s.size());
	for (i=0; i<s.size(); ++i)
	{
		s1[i]=revs(s[s.size()-1-i]);
	}
	return s1;
}

// convert string data_type of the motifs to int data_type
int motif2Digital(const string& seq)
{
	string s=Letter2Digital(seq);
	int index=0;
	for(int i=0;i<(int)seq.size();++i)
		index += (s[i]-48)*(int)pow(4.0,i);
	return index;
}

// convert digital data_type to string data_type motif, given the index and the length of motif
// mainly for output
string digital2Motif(const int index, const int len)
{
	string motif(len, '0');
	int quot = index;
	int order = 0;
	while(quot!=0) {
		motif[order] += (quot%4);
		quot /= 4;
		order++;
	}
	return Digital2Letter(motif);
}
