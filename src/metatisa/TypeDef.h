#ifndef TYPEDEF_H
#define TYPEDEF_H

#include"TypeDefBase.h"

#include"Matrix.h"
typedef Matrix_T<double> M_D;
typedef std::vector<Matrix_T<double> > Ve_M_D;
typedef std::pair< Ve_D, M_D > Pa_Ve_D_M_D;

class Location_T								
{
public:
	Location_T() : isPositive( false ), ORFLength( 0 ),LongestORFLength(0),
		tag("A")
	{
	}
	Location_T( Pa_I_I& rhs, bool positive = false ) 
		: isPositive( positive ),LongestORFLength(0), tag("A")
	{
		location = rhs;
	}

	bool operator <( const Location_T& rhs )const
	{
		return location.second < rhs.location.second;
	}

	bool operator ==( const Location_T& rhs )const
	{
		if( isPositive != rhs.isPositive )
			return false;
		return isPositive ? location.second == rhs.location.second
			: location.first == rhs.location.first;
	}
	
	void location2FileForm( int seqLen )
	{
		if( isPositive )
		{
			location.first += 1;
			location.second += 3;
		}
		else
		{
			int tmp = location.first;
			location.first = abs( seqLen - location.second ) - 2;
			location.second = abs( seqLen - tmp );
		}
	}
	void fileForm2Location( int seqLen )
	{
		if( isPositive )
		{
			location.first -= 1;
			location.second -= 3;
		}
		else
		{
			int tmp = location.first;
			location.first = abs( seqLen - location.second );
			location.second = abs( seqLen - tmp ) - 2;
		}
	}

	Pa_I_I location;//|A|TG........|T|TA。
	bool isPositive;

	int start;
	int trueTISPos; //测试用 记录真实的起始位点

	int ORFLength;
	int seqLen;
	int LongestORFLength;
	Str tag;

	
};

class ORF_T : public Location_T
{
public:
	ORF_T() {
		trueScore = -1.0, upFalseScore = -1.0, downFalseScore = -1.0;
		leftmost = -1, haveTIS = -10;
		have_XTG = false, have_upstream_stop = false, is_leftmost = false, enough_upbp_1stXTG = false;
	}

	Str* codingSeq;
	double trueScore, upFalseScore, downFalseScore;
	Pa_I_I origLocation; //储存Location 对于location.first>=3 与校正后的Location比较 会有更新
						//对于location.first<3 不更新

	int leftmost; //最长ORF的location.first

	Str id; //记录其在哪个片段上
	int haveTIS; // =-1时表示这段ORF找不到TIS

	// ycliu
	bool have_XTG;					// candidate TIS
	bool have_upstream_stop;		// upstream same frame stop
	bool enough_upbp_1stXTG;		// the first candidate of start have enough upstream sequence for scoring
	bool is_leftmost;				
};


typedef std::vector< Location_T > Ve_Location; 
typedef const std::vector< Location_T > con_Ve_Location; 
typedef std::vector< Location_T > Ve_Loca;  
typedef const std::vector< Location_T > con_Ve_Loca; 
typedef std::list< Location_T > Li_Loca;  
typedef std::set< Location_T > Set_Loca;
using std::cout;
using std::endl;
using std::setw;
using std::make_pair;

#endif //TYPEDEF_H
