
#include "seqop.h"

char SeqOp::Letter2Digital(const char res)
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
		assert("!UNEXPECTED nucleotide acid!");
		return '$';
	}
}

string SeqOp::Letter2Digital(const std::string &seq)
{
	string revSeq = seq;
	size_t i=0;
	for( i=0; i<seq.size(); ++i )
		revSeq[i] = SeqOp::Letter2Digital(seq[i]);
	return revSeq;
}

char SeqOp::Digital2Letter (const char res)
{
	switch(res)
	{
	case '0' : case 'A' : case 'a' : return 'A';
	case '1' : case 'C' : case 'c' : return 'C';
	case '2' : case 'G' : case 'g' : return 'G';
	case '3' : case 'T' : case 't' : return 'T';
	default  : 
		assert("!UNEXPECTED nucleotide acid!");
		return '$';
	}
}

string SeqOp::Digital2Letter(const string& seq)
{
	string revSeq = seq;
	size_t i=0; 
	for( i=0; i<revSeq.size(); ++i )
		revSeq[i] = SeqOp::Digital2Letter(seq[i]);
	return revSeq;
}

char SeqOp::revs(const char res)
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
		assert("!UNEXPECTED nucleotide acid!");
		return '$';
	}
}

string SeqOp::revs(const string& s)
{
	size_t i=0;
	string s1;
	s1.resize(s.size());
	for (i=0; i<s.size(); ++i)
	{
		s1[i] = SeqOp::revs(s[s.size()-1-i]);
	}
	return s1;
}

int SeqOp::motif2Digital(const string& seq)
{
	string s=Letter2Digital(seq);
	int index=0;
	for(int i=0;i<(int)seq.size();++i)
		index += (s[i]-48)*(int)pow(4.0,i);
	return index;
}

string SeqOp::digital2Motif(const int index, const int len)
{
	string motif(len, '0');
	int quot = index;
	int order = 0;
	while(quot!=0) {
		motif[order] += (quot%4);
		quot /= 4;
		order++;
	}
	return SeqOp::Digital2Letter(motif);
}

string SeqOp::drawseq(const string& genome,int start,unsigned length) 
{ 
	if (start<0) start+=(int)genome.size();
	if (start>(int)genome.size()) start-=(int)genome.size();
	if (start+length<=genome.size()) return genome.substr(start,length);
	else return genome.substr(start,genome.size()-start)+genome.substr(0,start+length-genome.size());
}

string SeqOp::drawORF(const string& genome, int left, int right, char strand)
{
	bool pos = (strand == '+' ? true : false);
	if(pos) return SeqOp::drawseq(genome, left, right-left+1);
	else return SeqOp::revs(drawseq(genome, left, right-left+1));
}

string SeqOp::toAminoSeq(const char* seq, int startPos, int endPosition)
{
	double sum = 0;
	string tmp;
	int i = startPos;
	for(; i <= endPosition - 3; i+=3 )
	{
		if( (seq[i]>'9' || seq[i]<'0') || (seq[i+1]>'9' || seq[i+1]<'0') || (seq[i+2]>'9' || seq[i+2]<'0') )
			continue;
		++sum;
		switch( ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48 )
		{
		case 36 : case 37 : case 38 : case 39 :
			tmp += 'A'; break;
		case 57 : case 59 :
			tmp += 'C'; break;
		case 33 : case 35 :
			tmp += 'D'; break;
		case 32 : case 34 :
			tmp += 'E'; break;
		case 61 : case 63 :
			tmp += 'F'; break;
		case 40 : case 41 : case 42 : case 43 :
			tmp += 'G'; break;
		case 17 : case 19 :
			tmp += 'H'; break;
		case 12 : case 13 : case 15 : 
			tmp += 'I'; break;
		case 0 : case 2 : 
			tmp += 'K'; break;
		case 28 : case 29 : case 30 : case 31 : case 60 : case 62 :
			tmp += 'L'; break;
		case 14 : 
			tmp += 'M'; break;
		case 1 : case 3 : 
			tmp += 'N'; break;
		case 20 : case 21 : case 22 : case 23 :
			tmp += 'P'; break;
		case 16 : case 18 : 
			tmp += 'Q'; break;
		case 8 : case 10 : case 24 : case 25 : case 26 : case 27 :
			tmp += 'R'; break;
		case 9 : case 11 : case 52 : case 53 : case 54 : case 55 :  
			tmp += 'S'; break;
		case 4 : case 5 : case 6 : case 7 : 
			tmp += 'T'; break;
		case 44 : case 45 : case 46 : case 47 :
			tmp += 'V'; break;
		case 58 :
			tmp += 'W'; break;
		case 49 : case 51 :
			tmp += 'Y'; break;
		case 48 : case 50 : case 56 :
		{
			--sum;
			break;
		}
		default :;
		}
	}
	return tmp;
}

vector<double> SeqOp::aminoEDP(const char* seq, int startPos, int endPosition)
{
	int size = 20;
	std::vector< double > EDP( size, 0 );
	double sum = 0;
	int i = startPos;
	for(; i <= endPosition - 3; i+=3 )
	{
		if( (seq[i]>'9' || seq[i]<'0') || (seq[i+1]>'9' || seq[i+1]<'0') || (seq[i+2]>'9' || seq[i+2]<'0') )
			continue;
		++sum;
		switch( ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48 )
		{
		case 36 : case 37 : case 38 : case 39 :
			++EDP[0]; break;	//tmp += 'A'; break;
		case 57 : case 59 :
			++EDP[1]; break;	//tmp += 'C'; break;
		case 33 : case 35 :
			++EDP[2]; break;	//tmp += 'D'; break;
		case 32 : case 34 :
			++EDP[3]; break;	//tmp += 'E'; break;
		case 61 : case 63 :
			++EDP[4]; break;	//tmp += 'F'; break;
		case 40 : case 41 : case 42 : case 43 :
			++EDP[5]; break;	//tmp += 'G'; break;
		case 17 : case 19 :
			++EDP[6]; break;	//tmp += 'H'; break;
		case 12 : case 13 : case 15 : 
			++EDP[7]; break;	//tmp += 'I'; break;
		case 0 : case 2 : 
			++EDP[8]; break;	//tmp += 'K'; break;
		case 28 : case 29 : case 30 : case 31 : case 60 : case 62 :
			++EDP[9]; break;	//tmp += 'L'; break;
		case 14 : 
			++EDP[10]; break;	//tmp += 'M'; break;
		case 1 : case 3 : 
			++EDP[11]; break;	//tmp += 'N'; break;
		case 20 : case 21 : case 22 : case 23 :
			++EDP[12]; break;	//tmp += 'P'; break;
		case 16 : case 18 : 
			++EDP[13]; break;	//tmp += 'Q'; break;
		case 8 : case 10 : case 24 : case 25 : case 26 : case 27 :
			++EDP[14]; break;	//tmp += 'R'; break;
		case 9 : case 11 : case 52 : case 53 : case 54 : case 55 :  
			++EDP[15]; break;	//tmp += 'S'; break;
		case 4 : case 5 : case 6 : case 7 : 
			++EDP[16]; break;	//tmp += 'T'; break;
		case 44 : case 45 : case 46 : case 47 :
			++EDP[17]; break;	//tmp += 'V'; break;
		case 58 :
			++EDP[18]; break;	//tmp += 'W'; break;
		case 49 : case 51 :
			++EDP[19]; break;	//tmp += 'Y'; break;
		case 48 : case 50 : case 56 :
		{
			--sum;
			break;
		}
		default :;
		}
	}

	int j = 0;
	for(; j < size; ++j )
	{
		if( EDP[j] != 0 )
		{
			double pi = EDP[j] / sum;
			EDP[j] = pi * log( pi );
		}
	}
	sum = 0;
	for( j = 0; j < size; ++j )
		sum += EDP[j];

	if( sum != 0 )
	{
		for( j = 0; j < size; ++j)		
			EDP[j] /= sum;
	}
	return EDP;
}

vector<double> SeqOp::codonEDP(const char* seq, int startPos, int endPosition)
{
	int size = 64;
	std::vector< double > EDP( size, 0 );
	double sum = 0;
	int i = startPos;
	for(; i <= endPosition - 3; i+=3 )	// we mustn't count for the stop codon
	{
		if( (seq[i]>'9' || seq[i]<'0') || (seq[i+1]>'9' || seq[i+1]<'0') || (seq[i+2]>'9' || seq[i+2]<'0') )
			continue;
		++sum;
		int codon = ((seq[i] - 48)<<4) + ((seq[i + 1] - 48)<<2) + seq[i + 2] - 48;
		EDP[ codon ]++;
	}
	EDP.erase( EDP.begin()+56 );	// delet the element of "TGA"
	EDP.erase( EDP.begin()+50 );	// delet the element of "TAG"
	EDP.erase( EDP.begin()+48 );	// delet the element of "TAA"

	size = (int)EDP.size();
	int j = 0;
	for(; j < size; ++j )
	{
		if( EDP[j] != 0 )
		{
			double pi = EDP[j] / sum;
			EDP[j] = pi * log( pi );
		}
	}
	sum = 0;
	for( j = 0; j < size; ++j )
		sum += EDP[j];

	if( sum != 0 )
	{
		for( j = 0; j < size; ++j)		
			EDP[j] /= sum;
	}
	return EDP;
}

double SeqOp::gccontent(const std::string &seq)
{
	string tmpseq = SeqOp::Letter2Digital(seq);
	double gc = 0, total = 0;
	for(int i=0; i<(int)tmpseq.size(); ++i)
	{
		++ total;
		if(tmpseq[i] == '1' || tmpseq[i] == '2')
			++ gc;
	}
	return gc / total;
}

double SeqOp::theta(const std::string &seq)
{
	string tmpseq = SeqOp::Letter2Digital(seq);
	double gc13 = 0, gc = 0.00001;
	for(int i=0; i<(int)tmpseq.size()-2; i+=3)
	{
		if(tmpseq[i] == '1' || tmpseq[i] == '2')
		{
			++ gc13;
			++ gc;
		}
		if(tmpseq[i+1] == '1' || tmpseq[i+1] == '2')
			++ gc;
		if(tmpseq[i+2] == '1' || tmpseq[i+2] == '2')
		{
			++ gc13;
			++ gc;
		}
	}
	return gc13/gc;
}
