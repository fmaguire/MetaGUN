#ifndef METAGENOME_H
#define METAGENOME_H

#include <math.h>
#include "TypeDef.h"
#include "SequenceTransform.h"
#include "OftenUsedOperatLib.h"

class Metagenome_T
{
public:
	std::map< Str, std::pair<Str,Str> > IdSeqMap;
	std::map< Str, std::vector<Location_T> > IdLocationMap;
	std::multimap< Str, Str > GenusIdMap;
	
	std::map< Str, std::vector<ORF_T> > GenusORFSetMap;

	void getSeqs( Str& seqsFilename );
	void getLocations( Str& LocationsFilename );
	void getGenusId(Str &seqsFilename, Str &bin_model_file, Str &taxon_map_file);
	void getGenusId(const std::string& binfile);
	void getGenusORFSetMap ();
	
	void outSeqs( Str& seqsFilename );
	void outLocations( Str& LocationsFilename );
	void outGenusId ( Str& GenusIdFilename );
	void outORFSet( Str& ORFSetFilename );
	void outGenusID(const std::string& bin_result_file);

	M_D getMatrix( Str& MatrixFilename);
	void displayMatrix( const Matrix_T< double >& ma, std::ostream& out );
};

class TIS_T:public Metagenome_T
{
public:
	//Settings
	int UPBP, DOWNBP, MAXORDER, ORDER;
	int HaveHeadNumCutoff;
	Set_I STARTS, STOPS;
	
	//����������
	static const int maximalInterationNum = 15;
	int minimum_upstream;
//	static const int minimum_upstream = 8;

	//Paramters to be calculated
	M_D trueTIS, upFalse, downFalse;                              
	int maxCandidateN;
	double pFU, pFD, pT;      
	Ve_D qi;

	static const bool filter_on = true;
	double cutoff;   //�ж�û��ͷ��ORF��UPBP�Ժ��һ��XTG�Ƿ���coding������ֵ
	std::string prioriparameters;

	//Initiate parameters
	void initiateSettings(const Str& infile);
	void getprioriParameters( Str& genusName );
	void initiatePWMs( std::vector<ORF_T>& ORFSet ); //ģ������Ԥ��ʱ�򿪹�����
	void getcutoff( std::vector<ORF_T>& ORFSet );

	void x2LongestORF( const char* seq, Pa_I_I& location );
	void x2LongestORF();
	
	//Make prediction
	void reviseTIS_HaveHead_Enough( std::vector<ORF_T>& ORFSet,int HaveHeadNum );     //У����ͷ��ORF �����ͷ��Ŀ�㹻 ��ѧϰPWM
	void reviseTIS_HaveHead_NotEnough( std::vector<ORF_T>& ORFSet );				  //�����ͷ��Ŀ���� �����Ѿ�ѧ�õ�genus�������
	void reviseTIS_NoHead( std::vector<ORF_T>& ORFSet );
	void reviseTIS();

	// ycliu: new revise function 
	bool is_windows;
	void new_reviseTIS();
	void revise_head_enough(std::vector<ORF_T>& ORFSet, int head_enough_number);
	void revise_nohead(std::vector<ORF_T>& ORFSet);
	void revise_not_enough(std::vector<ORF_T>& ORFSet);
	void scoring_XTG(const char* ssseq, double& t, double& d, double& u);
	// scoring for candidate TISs with insufficient upstream sequences, start means the when to start for scroing by PWMs
	void scoring_XTG(const char* ssseq, double& t, double& d, double& u, int start);
	void set_ordinary_starts();
	void set_ordinary_stops();
	void set_abnormal_stops();
	std::string generate_seq(const std::vector<double>& gc, int length);
	std::vector<double> gc_content(const char* seq);
	void read_prior_tismodel(const std::string& tismodel_file);
							 	
	//evaluate coverage
	//void getCoverage( Str& resultFilename); //�������ݵ�Ԥ�⾫�� ��Ҫ��LocationFile�ж�ȡstart��trueTISPos��Ϣ
	//void getCoverageForWholeGenome(Str& AnnotationFilename,Str& resultFilename); //ȫ���������ʱ��coverage �൱��TriTISA
	
	//load and save results                                                                                            
	void resultToFile( Str& resultFile, Str& resultFormat );                                                                                                                              
	//void parametersOut(Str& parameters); //ѧϰ��ÿ��genus�Ĳ��� ���һ����������ļ� ���������ļ���2��markovѧϰ�����
};

#endif
