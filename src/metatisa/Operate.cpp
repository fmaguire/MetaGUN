
#include "Operate.h"
#include <math.h>
#include <fstream>
#include <iostream>

string ReadFNA(const string& fna)
{
	std::ifstream fin;
	fin.open(fna.c_str());
	if( !fin.good() ) {
		cout << "can't open " << fna << endl;
		exit(0);
	}

	char c[200];
	fin.getline(c,200);	
	std::string seq="";
	seq.reserve(10000000);
	string s;
	while (!fin.eof())
	{
		fin>>s;
		seq+=s;
	}
	fin.close();
	return seq;
}

string drawseq(const string& genome,int start,unsigned length) 
{ 
	if (start<0) start+=(int)genome.size();
	if (start>(int)genome.size()) start-=(int)genome.size();
	if (start+length<=genome.size()) return genome.substr(start,length);
	else return genome.substr(start,genome.size()-start)+genome.substr(0,start+length-genome.size());
}

string drawORF(const string& genome, int left, int right, char strand)
{
	bool pos = (strand == '+' ? true : false);
	if(pos) return drawseq(genome, left, right-left+1);
	else return revs(drawseq(genome, left, right-left+1));
}

vector<double> AAEDP(const char* seq, int startPos, int endPosition)
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


vector<double> CodonEDP(const char* seq, int startPos, int endPosition)
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
