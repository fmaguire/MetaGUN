//	created:	2008-11-17
//  author:		Jiangtao Guo
//	Metagenome TIS prediction

#include <ctime>
#include "TypeDef.h"
#include "SequenceTransform.h"
#include "OftenUsedOperatLib.h"
#include "Metagenome.h"
#include "OftenUsedFun.h"
using   namespace   std;

string shotgunSeqFile, shotgunLocationsFile, bin_model_file, taxon_map_file, prioriparameters;
string bin_result_file;
string resultFile = "metagenome.med";
string resultFormat = "GFF";
string settingsFile = "metatisa-settings.txt";
bool windows_platform = true;
int minimum_upstream_length = 8;
 
void exit_with_help()
{
	cout<< "MetaTISA     2012-10-10" << endl
		<< "MetaTISA [-s meta-sequence] [-l meta-location] [-b bin-result-file]\n"
		<< "         [-p path-of-prior-parameters] [-f output-format] [-o output-file]\n"
		<< "Options:\n"
		<< "  -s FILE        metagenome sequence file in FASTA format\n"
		<< "  -l FILE        gene location file in MED format\n"
		<< "  -b FILE        file stored the binning result\n"
		<< "  -p STRING      path of directory stores the priori TIS modelsparameters for PWM\n"
		<< "  -o FILE        output filename, 'metagenome.gff' by default\n"
		<< "  -f STRING      Format of result file, currently MED and GFF available, 'GFF' by default\n"
		<< "  -c FILE        settings file, 'metatisa-settings.txt' by default\n"
		<< "  -w INT         1: windows platform, 0: linux platform; default, 1, windows\n"
		<< "  -m INT         minimum upstream scoring length, default 8,\n";
	exit(1);
}

int main(int argc, char* argv[])
{
	//parse options
	if(argc == 1) {
		exit_with_help();
	}
	int i = 1;
	for(; i<argc; ++i)
	{
		if(argv[i][0] != '-') continue;
		if(++i >= argc) 
		{
			cout << "Not specified option " << argv[i-1] << endl;
			exit_with_help();
		}
		switch(argv[i-1][1])
		{
		case 's' : 
			shotgunSeqFile = string(argv[i]);
			break;
		case 'l' :
			shotgunLocationsFile = string(argv[i]);
			break;
		case 'b':
			bin_result_file = string(argv[i]);
			break;
		case 'M' :
			bin_model_file = string(argv[i]);
			break;
		case 't' :
			taxon_map_file = string(argv[i]);
			break;
		case 'p' :
			prioriparameters = string(argv[i]);
			break;
		case 'c':
			settingsFile = string(argv[i]);
			break;
		case 'w':
			if(atoi(argv[i]) == 0) {
				windows_platform = false;
			} else {
				windows_platform = true;
			}
			break;
		case 'f' :
			resultFormat = string(argv[i]);
			if (resultFile == "metagenome.med")
			{
				if (resultFormat == "MED" || resultFormat == "med")
					resultFile = "metagenome.med";
				else
					resultFile = "metagenome.gff";
			}
			break;
		case 'o' :
			resultFile = string(argv[i]);
			break;
		case 'm':
			minimum_upstream_length = atoi(argv[i]);
			break;
		default :
			cout << "Unknown Option!" << endl;
			exit_with_help();
		}
	}
	if(i>argc) exit_with_help();

	if(shotgunSeqFile.empty())
	{
		cout<< "Not enough arguments:\n"
			<< "     shotgunSeqFile is not defined\n";
		exit_with_help();
	}
	if(shotgunLocationsFile.empty())
	{
		cout<< "Not enough arguments:\n"
			<< "     shotgunLocationsFile is not defined\n";
		exit_with_help();
	}

	if(prioriparameters.empty())
	{
		cout<< "Not enough arguments:\n"
			<< "     priori parameters is not defined\n";
		exit_with_help();
	}
	if(resultFormat!="MED" && resultFormat!="GFF")
	{
		cout << "Wrong Output Format, currently only MED and GFF format is available\n\n";
		exit_with_help();
	}

	TIS_T MetagenomeTIS;
	MetagenomeTIS.minimum_upstream = minimum_upstream_length;
	cout << "Minimum upstream scoring length is: " << MetagenomeTIS.minimum_upstream << endl;
	MetagenomeTIS.prioriparameters = prioriparameters;
	MetagenomeTIS.is_windows = windows_platform;
	if(MetagenomeTIS.is_windows) {
		cout << "MetaTISA is running on Windows platform." << endl;
	} else {
		cout << "MetaTISA is running on Linux platform." << endl;
	}

	cout << "Getting sequences and input locations ..." << endl;
	MetagenomeTIS.getSeqs(shotgunSeqFile);

	cout << "Getting sequence binning informations ..." << endl;
	if(bin_result_file.empty()) {
		MetagenomeTIS.getGenusId(shotgunSeqFile, bin_model_file, taxon_map_file);
	} else {
		MetagenomeTIS.getGenusId(bin_result_file); 
	} 
	
	MetagenomeTIS.getLocations(shotgunLocationsFile);
	
	MetagenomeTIS.getGenusORFSetMap();

	MetagenomeTIS.initiateSettings(settingsFile);
	MetagenomeTIS.x2LongestORF();

	cout << "Revising gene starts ... " << endl;
	MetagenomeTIS.new_reviseTIS();
	MetagenomeTIS.resultToFile(resultFile, resultFormat);
	return 0;
}
