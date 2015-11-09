#include <math.h>
#include <sstream>
#include <fstream>
#include "TypeDef.h"
#include "SequenceTransform.h"
#include "OftenUsedOperatLib.h"
#include "Metagenome.h"
#include "Bin.h"
#include "OftenUsedFun.h"

M_D Metagenome_T::getMatrix(Str &MatrixFilename)
{
	std::ifstream in( MatrixFilename.data() );
	if( !in )
	{ cout << "can't open file " << MatrixFilename << endl; exit(0); }

	int lin = 0;
	int col = 0;
	std::string line;

	in >> lin >> col;
	getline( in, line );

	Matrix_T< double > ma( lin, col, 0.0 );
	for( int j=0;  j<col; ++j )
	{
		for( int i=0; i<lin; ++i )
			in >> ma( i, j );
		getline( in, line );
	}
	in.close();
	return ma;
}

void Metagenome_T::displayMatrix( const Matrix_T< double >& ma, std::ostream& out )
{
	out << std::setw(10) << ma.getLine() << std::setw(10) << ma.getColum() << endl;
	for( int i=0; i<ma.getColum(); ++i )
	{
		for( int j=0; j<ma.getLine(); ++j )
			out << ma(j,i) << "\t";
		out << endl;
	}
}

void TIS_T::x2LongestORF( const char* seq, Pa_I_I& location )
{
	int STPPos = location.second;
	int ATGPos = -1;
	int tmp = ((seq[STPPos] - 48)<<4) + ((seq[STPPos+1] - 48)<<2) + seq[STPPos + 2] - 48;

	for( STPPos -= 3; ; STPPos -= 3 )
	{
		if( STPPos < 0 )
		{
			location.first = ATGPos;
			return;
		}
		int subStr = ((seq[STPPos] - 48)<<4) + ((seq[STPPos+1] - 48)<<2) + seq[STPPos + 2] - 48;
		//			TAA-300			TGA-320			TAG-302
		if( STOPS.find(subStr) != STOPS.end() ){
			location.first = ATGPos;
			return;
		}
		//			ATG-032			GTG-232          TTG-332        CTG-132
		if( STARTS.find(subStr) != STARTS.end() )
			ATGPos = STPPos;
	}
	location.first = ATGPos;
}

void TIS_T::x2LongestORF()
{
	set_ordinary_starts();

	std::map<Str, std::vector<ORF_T> >::iterator iter;
	iter = GenusORFSetMap.begin();
	for(; iter != GenusORFSetMap.end(); ++iter) 
	{
		Str genus = iter->first;
		if ( genus == "Mycoplasma" || genus == "Acholeplasma" || genus == "Aster" || genus == "Onion" || genus == "Ureaplasma" ) {
			set_abnormal_stops();
		} else {
			set_ordinary_stops();
		}

		std::vector<ORF_T>& ORFSet = (iter->second);
		for ( int i=0 ; i < (int)ORFSet.size(); i++ )
		{
			Pa_I_I tmp  = ORFSet[i].location;
			if (tmp.first > 2)
			{
				x2LongestORF(ORFSet[i].codingSeq->data(), tmp);
			}
			ORFSet[i].leftmost = tmp.first; //记录最开始的location.first 后面迭代时要用
			ORFSet[i].location.first = tmp.first;  // relocated the start to the longest ORF -- by yc.liu
			ORFSet[i].origLocation = ORFSet[i].location;
			ORFSet[i].is_leftmost = true;
			if(ORFSet[i].location.first != -1)
				ORFSet[i].have_XTG = true;
			if(ORFSet[i].location.first > 2)
				ORFSet[i].have_upstream_stop = true;
			if(ORFSet[i].location.first >= UPBP)
				ORFSet[i].enough_upbp_1stXTG = true;
		}
	}
}

void Metagenome_T::getSeqs(Str& seqsFilename)
{
	std::ifstream seqsFile( seqsFilename.data() );
	if( !seqsFile.good() ){
		std::cout<<"file "<<seqsFilename<<" not found!"<<std::endl;
		exit(1);
	}

	Str line;
	Str seq,negativeseq;
	Str id;

	std::getline( seqsFile, line );
	while(!seqsFile.eof()){
		
		id = line.substr(1);
		getline( seqsFile,line );
		while ( line.find(">") == std::string::npos )
		{
			if ( seqsFile.eof())
				break;
			seq+=SequenceTransform_T::char2DigitalSeq(line);
			getline( seqsFile, line);
		}

		negativeseq = seq;
		std::reverse( negativeseq.begin(), negativeseq.end() );	//将序列翻转。
		std::for_each( negativeseq.begin(), negativeseq.end(),	//转换为互补链。
			SequenceTransform_T::ToOppRule_T() );

		IdSeqMap.insert( make_pair(id, make_pair(seq,negativeseq)) );

		seq.erase( seq.begin(), seq.end() );
		negativeseq.erase( negativeseq.begin(), negativeseq.end() );
	}	
	seqsFile.close();
}

void Metagenome_T::outSeqs( Str& seqsFilename )
{
	std::ofstream seqsFile( seqsFilename.data() );
	std::map< Str, std::pair<Str,Str> >::iterator iter = IdSeqMap.begin();
	for( ; iter!=IdSeqMap.end() ; iter++)
	{
		seqsFile << ">" << iter->first << endl;

		Str seq = SequenceTransform_T::digital2CharSeq(iter->second.first);
		int k;
		for( k = 0; k < (int)seq.size()/70; k++ )
			seqsFile << seq.substr(70*k,70) << endl;
		if ( seq.substr(70*k).size() != 0 )    //有可能最后还有一小段序列
			seqsFile << seq.substr(70*k) << endl;
	}
	seqsFile.close();
}

void Metagenome_T::getLocations(Str &LocationsFilename)
{
	std::ifstream LocationsFile( LocationsFilename.data() );
	if( !LocationsFile.good() ){
		std::cout<<"file "<<LocationsFilename<<" not found!"<<std::endl;
		exit(1);
	}

	Str line;
	Str id;
	Str pos1,pos2,tmp;
	bool isPositive;

	std::vector<Location_T> subLocations;

	std::getline( LocationsFile, line );
	
	
	while(!LocationsFile.eof()){
		id = line.substr(1);
		getline( LocationsFile, line );
		while( line.find(">") == std::string::npos )
		{

			std::istringstream in ( line );
			in >> pos1 >> pos2 >> tmp;

			isPositive = (tmp=="+");

			Pa_I_I location = make_pair(atoi(pos1.data()),atoi(pos2.data()));
			Location_T Location(location,isPositive);

			subLocations.push_back(Location);

			std::getline( LocationsFile, line );
			
			if(LocationsFile.eof())
				break;
		}

		if (subLocations.size()!=0)
		{
			IdLocationMap.insert( make_pair(id,subLocations) );
			subLocations.erase( subLocations.begin(), subLocations.end() );
		}
	}
	
	LocationsFile.close();
}

void Metagenome_T::outLocations( Str& LocationsFilename )
{
	std::ofstream LocationsFile( LocationsFilename.data() );
	std::map< Str, std::vector<Location_T> >::iterator iter = IdLocationMap.begin();
	for( ; iter!=IdLocationMap.end(); iter++ )
	{
		LocationsFile << ">" << iter->first << endl ;
		std::vector<Location_T> subLocations = iter->second;
		for( int j = 0; j < (int)subLocations.size(); j++ )
		{
			LocationsFile << subLocations[j].location.first << "\t" << subLocations[j].location.second << "\t";
			if ( subLocations[j].isPositive )
				LocationsFile << "+";
			else
				LocationsFile << "-";
			LocationsFile << endl;
		}
	}
	LocationsFile.close();
}

void Metagenome_T::getGenusId(Str &seqsFilename, Str &bin_model_file, Str &taxon_map_file)
{
	std::vector<string> genome_list;
	std::vector<std::vector<float> > table;
	int start = (int)time(NULL);
	get_Freqs_table(bin_model_file, genome_list, table);
	cout << time(NULL)-start << " seconds costs" << endl;

	std::vector<string> genomes;
	std::vector<string> groups;
	get_group_map(taxon_map_file, genomes, groups);
	Metagenome_T::GenusIdMap = bin_fasta(table, IdSeqMap, genomes, groups);
}

void Metagenome_T::getGenusId(const std::string& binfile)
{
	vector<string> lines = readlines(binfile);
	std::multimap<std::string, std::string> genus_id_map;
	for(int i=0; i<(int)lines.size(); ++i) {
		vector<string> vv = split_string(lines[i]);
		genus_id_map.insert(std::make_pair(vv[1], vv[0].substr(1)));
	}
	Metagenome_T::GenusIdMap = genus_id_map;
}

void Metagenome_T::outGenusID(const string& bin_result_file)
{
	std::ofstream fout(bin_result_file.data());
	if(!fout) { cout << "can't open " << bin_result_file << endl; exit(-1); }
	std::multimap<string, string>::iterator itergi = Metagenome_T::GenusIdMap.begin();
	for(; itergi != Metagenome_T::GenusIdMap.end(); ++ itergi)
		fout << itergi->first << '\t' << itergi->second << endl;
	fout.close();
}

void Metagenome_T::outGenusId(Str &GenusIdFilename)
{
	std::ofstream GenusIdFile( GenusIdFilename.data() );
	std::multimap< Str, Str >::iterator iter =GenusIdMap.begin();
	for( ; iter!=GenusIdMap.end(); iter++ )
		GenusIdFile << ">" << iter->second << "\t" << iter->first << endl ;

	GenusIdFile.close();
}

void Metagenome_T::getGenusORFSetMap()
{
	std::multimap< Str, Str >::iterator iter1 = GenusIdMap.begin();
	std::map< Str, std::pair<Str,Str> >::iterator iter2 = IdSeqMap.begin();
	std::map< Str, std::vector<Location_T> >::iterator iter3 = IdLocationMap.begin();

	for ( ; iter1!=GenusIdMap.end(); iter1++ )
	{
		Str id = iter1->second;
		if ( IdSeqMap.find( id )!= IdSeqMap.end() ) 
			iter2 = IdSeqMap.find( id ) ;
		else 
			continue;
		if ( IdLocationMap.find(id) != IdLocationMap.end() )
			iter3 = IdLocationMap.find(id);
		else
			continue;

		for ( int j=0; j<(int)iter3->second.size(); j++ )
		{
			ORF_T ORF;
			ORF.id = iter1->second;
			ORF.isPositive = (iter3->second)[j].isPositive;
			ORF.location = (iter3->second)[j].location;
			ORF.fileForm2Location((int)iter2->second.first.size()); //在这里即将location转换为后面计算要用的形式

			ORF.start = (iter3->second)[j].start;
			ORF.trueTISPos = (iter3->second)[j].trueTISPos;

			if (ORF.isPositive)
				ORF.codingSeq = &iter2->second.first;
			else
				ORF.codingSeq = &iter2->second.second;

			ORF.ORFLength =  ORF.location.second - ORF.location.first ; //此时的location.second指向的是终止子的第一位

			ORF.seqLen = (int)iter2->second.first.size();

			GenusORFSetMap[iter1->first].push_back(ORF);
		}
	}
}

void TIS_T::initiateSettings(const string& infile)
{
	std::ifstream in(infile.data());
	if( !in.good() ){
		cout<<"The settings file " << infile << " is missing in this folder."<<endl;
		exit(0);
	}
	while( !in.eof() ){
		Str line;
		std::getline(in,line);
		if( !line.empty() && (line.find("#") == std::string::npos) ){
			int pos = (int)line.find("=");
			Str key = line.substr(0,pos), value = line.substr(pos+1);
			if( key == "MAXORDER" ){
				MAXORDER = atoi( value.data() );
			}
			if( key == "UPSTREAM" ){
				UPBP = atoi( value.data() );
			}
			if( key == "DOWNSTREAM" ){
				DOWNBP = atoi( value.data() );
			}
			if( key == "HaveHeadNumCutoff" ){
				HaveHeadNumCutoff = atoi(value.data());
			}			
		}
	}
}

void TIS_T::read_prior_tismodel(const std::string &tismodel_file)
{
	std::ifstream fin(tismodel_file.data());
	if(!fin) { cout << "can't open " << tismodel_file << endl; exit(-1); }
	std::string line;
	while(!fin.eof()) {
		line.clear();
		getline(fin, line);
		if(line.empty())
			continue;
		if(line == "-------Prior probabilites------") {
			fin >> pT >> pFU >> pFD >> cutoff;
			getline(fin, line);
		} 
		else if(line == "------True TIS Markov model------") {
			int lin=0, col=0;
			fin >> lin >> col;
			getline(fin, line);
			Matrix_T<double> mm(lin, col, 0.0);
			for(int c=0; c<col; ++c) {
				for(int l=0; l<lin; ++l)
					fin >> mm(l, c);
				getline(fin, line);
			}
			trueTIS = mm;
		}
		else if(line == "------Upfalse TIS Markov model------") {
			int lin=0, col=0;
			fin >> lin >> col;
			getline(fin, line);
			Matrix_T<double> mm(lin, col, 0.0);
			for(int c=0; c<col; ++c) {
				for(int l=0; l<lin; ++l)
					fin >> mm(l, c);
				getline(fin, line);
			}
			upFalse = mm;
		}
		else if(line == "------Downfalse TIS Markov model------") {
			int lin=0, col=0;
			fin >> lin >> col;
			getline(fin, line);
			Matrix_T<double> mm(lin, col, 0.0);
			for(int c=0; c<col; ++c) {
				for(int l=0; l<lin; ++l)
					fin >> mm(l, c);
				getline(fin, line);
			}
			downFalse = mm;
		}
	}
	fin.close();
}

void TIS_T::getprioriParameters( Str& genusName )
{
	//测试用
	//genusName = "Ecoli";

	Str path = prioriparameters + "\\" + genusName + "\\";

	Str prioriProbFilename = path+genusName+".prioriProb";
	Str trueMaFilename = path+genusName+".trueMa";
	Str ufMaFilename = path+genusName+".ufMa";
	Str dfMaFilename = path+genusName+".dfMa";

	std::ifstream prioriProbFile ( prioriProbFilename.data() );
	prioriProbFile >> pT >> pFU >> pFD >> cutoff;
	prioriProbFile.close();

	trueTIS = getMatrix(trueMaFilename);
	upFalse = getMatrix(ufMaFilename);
	downFalse = getMatrix(dfMaFilename);
}
void TIS_T::initiatePWMs( std::vector<ORF_T>& ORFSet )
{
	//ORDER = 0; // 暂时只考虑0阶马尔科夫
	int sz = (int)pow(4.0,ORDER+1);
	upFalse = M_D(sz, UPBP+3+DOWNBP,1); //用1做初始化 防止后面出现log0的情况
	downFalse = M_D(sz, UPBP+3+DOWNBP,1);
	trueTIS = M_D( sz, UPBP+3+DOWNBP,1);

	Set_Str upSet,downSet,trueSet; //后面ORF去重时使用

	Ve_D bkg(4,0);
	int C = 0;
	int upCount=0, downCount=0, trueCount=0;

	std::vector< ORF_T >::iterator iter = ORFSet.begin();
	for( ; iter != ORFSet.end(); ++iter )
	{
		const char* seq = iter->codingSeq->data();
		int seqLen = iter->seqLen;
		
		//如果location.first小于三，则该ORF几乎是没有头的ORF，真实TIS在该片段中的概率很小
		//不用此类ORF学习PWM，消除噪音影响
		if ( iter->location.first < 3 )
			continue;
//		if ( iter->ORFLength < 300 )
//				continue;
		
		int initTIS=-1; //记录初始TIS位置 按location.first向下走，找第一个同相位的TIS

		//按location.first向下走，找第一个同相位的TIS 作为初始真正TIS的位置
		int hint = iter->location.first ;

		for( ; hint < iter->location.second; hint += 3 )
		{
			int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
			if(STARTS.find(subStr) != STARTS.end())
			{
				initTIS = hint;
				break;
			}
		}

		//如果initTIS仍为-1 即未找到同相位XTG 
		if ( initTIS == -1 )
		{
			iter->haveTIS = -1;
			continue;
		}

		iter->location.first = initTIS; //把location.first设置为initTIS 
		
		//无过滤器
		if ( filter_on == false)
		{
			//up false TIS, regardless of frame
			hint = 0 > (initTIS - 200) ? 0 : (initTIS - 200);
			for( ; hint < initTIS - 3; ++hint )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end())
				{
					int beg = hint - UPBP;
					if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
					{
						upCount++;
						int i = 0;
						for( ; i < upFalse.getColum(); ++i ){
							int index = 0;
							for( int f = 0; f < ORDER + 1; ++f ){
								index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
							}
							++upFalse(index,i);							
							++bkg[seq[beg+i]-'0'];
							++C;
						}
					}
				}
			}
			//down false TIS
			hint = initTIS + 3;
			for( ; hint < iter->location.second; hint += 3 ){
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if(STARTS.find(subStr) != STARTS.end()){//)
					int beg = hint - UPBP;
					if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
					{
						downCount++;
						int i = 0;
						for( ; i < upFalse.getColum(); ++i ){
							int index = 0;
							for( int f = 0; f < ORDER + 1; ++f ){
								index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
							}
							++downFalse(index,i);							
						}
					}
				}
			}
			//true TIS
			int beg = initTIS - UPBP;
			if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
			{
				trueCount++;
				int i = 0;
				for( ; i < trueTIS.getColum(); ++i ){
					int index = 0;
					for( int f = 0; f < ORDER + 1; ++f ){
						index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
					}
					++trueTIS(index,i);
				}		
			}
		}
		//有过滤器
		else
		{
			Str nter;
			//up false TIS, regardless of frame
			hint = (0 > (initTIS - 200) ? 0 : (initTIS - 200));
			for( ; hint < initTIS - 3; ++hint )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( subStr == 14 || subStr == 46 || subStr == 62 || subStr == 30)
				{
					int beg = hint - UPBP;
					if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1  + ORDER < seqLen  )
					{
						upCount++;
						nter = Str(seq+beg,UPBP + 3 + DOWNBP) ;
						if( upSet.find(nter) == upSet.end() )
						{
							int i = 0;
							for( ; i < upFalse.getColum(); ++i ){
								int index = 0;
								for( int f = 0; f < ORDER + 1; ++f ){
									index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
								}
								++upFalse(index,i);							
								++bkg[seq[beg+i]-'0'];
								++C;
							}
							upSet.insert(nter);
						}
						
					}
				}
			}
			//down false TIS
			hint = initTIS + 3;
			for( ; hint < iter->location.second; hint += 3 ){
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if(STARTS.find(subStr) != STARTS.end()){
					int beg = hint - UPBP;
					if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1  + ORDER < seqLen  )
					{
						downCount++;
						nter = Str(seq+beg,UPBP + 3 + DOWNBP) ;
						if ( downSet.find(nter)==downSet.end() )
						{
							int i = 0;
							for( ; i < downFalse.getColum(); ++i ){
								int index = 0;
								for( int f = 0; f < ORDER + 1; ++f ){
									index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
								}
								++downFalse(index,i);							
							}
							downSet.insert(nter);
						}
					}
				}
			}
			//true TIS
			int beg = initTIS - UPBP;
			if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
			{
				trueCount++;
				nter = Str(seq+beg,UPBP + 3 + DOWNBP) ;
				if( trueSet.find(nter)==trueSet.end() )
				{
					int i = 0;
					for( ; i < trueTIS.getColum(); ++i ){
						int index = 0;
						for( int f = 0; f < ORDER + 1; ++f ){
							index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
						}
						++trueTIS(index,i);
					}
					trueSet.insert(nter);
				}
			}
		
		}
	}

	//cout << "UsedORF:\t" << countUsedORF << endl;
	/*cout << "ORFSet size:\t" << ORFSet.size() << endl;
	cout << "upCount:\t" << upCount << "\tdownCount:\t" << downCount << "\ttrueCount:\t" << trueCount << endl;
	if (filter_on)
		cout << "upAfterFiltered:\t" << upSet.size() << "\tdownAfterFiltered:\t" << downSet.size() << "\ttrueAfterFiltered:\t" << trueSet.size() << endl;*/

	//Normalize PWM
	trueTIS.toBeAveraged();
	upFalse.toBeAveraged();
	downFalse.toBeAveraged();
	trueTIS.toBeLoged();
	upFalse.toBeLoged();
	downFalse.toBeLoged();

	/*displayMatrix(upFalse,cout);
	displayMatrix(downFalse,cout);
	displayMatrix(trueTIS,cout);*/

	int i = 0;
	for( ; i < (int)bkg.size(); ++i )
		bkg[i] /= C;
	
	qi = Ve_D();
	double P = 0;
	double ufN = 0;
	double a = bkg[0], c = bkg[1], g = bkg[2], t = bkg[3];
	for( i = 0; ; ++i ){
		double p = (a*t*g+g*t*g+t*t*g)/(a*t*g+g*t*g+t*t*g+t*g*a+t*a*a+t*a*g);
		double q = (t*g*a+t*a*a+t*a*g)/(a*t*g+g*t*g+t*t*g+t*g*a+t*a*a+t*a*g);
		P += q*pow(p,i);
		qi.push_back(q*pow(p,i));
		if(P>0.9999) {
			maxCandidateN = i+1;
			break;
		}
		else
			ufN += i * q*pow(p,i);
	}

	pT = 1;
	pFU = ufN;
	pFD = maxCandidateN - 1 - pFU;
	//cout << "pFU:\t" << pFU << "\tpT:\t" << pT << "\tpFD:\t" << pFD << endl; 
}

void TIS_T::scoring_XTG(const char* ssseq, double& t, double& d, double& u)
{
	// set to zero
	int beg = 0;
	t = 0, d = 0, u = 0;

	int index = 0;
	int f = 0;
	// find the index of first (ORDER-1) nucleotide(s)
	for( ; f < ORDER ; ++f ){
		index += (int)pow(4.0,ORDER-f)*(ssseq[beg+f]-'0');
	}
	// get the score of the frist (ORDER-1) nucleotides(s)
	for( f = 0; f < 4; ++f ){
		t += exp(trueTIS(index+f,0)); 
		d += exp(downFalse(index+f,0));
		u += exp(upFalse(index+f,0));
	}
	// multipe the prior probabilities, as plus for log
	t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);

	// formulating the ORDER-th Markov model
	int i = 0;
	for( ; i < trueTIS.getColum(); ++i ){
		int index = 0;
		int f = 0;
		for( ; f < ORDER ; ++f ){
			index += (int)pow(4.0,ORDER-f)*(ssseq[beg+i+f]-'0');
		}
		double pt = 0, pt2 = 0, pd = 0, pu = 0;
		for( f = 0; f < 4; ++f ){
			pt += exp(trueTIS(index+f,i)); 
			pd += exp(downFalse(index+f,i));
			pu += exp(upFalse(index+f,i));
		}
		t += trueTIS(index+(ssseq[beg+i+ORDER]-'0'),i) - log(pt);
		d += downFalse(index+(ssseq[beg+i+ORDER]-'0'),i) - log(pd);
		u += upFalse(index+(ssseq[beg+i+ORDER]-'0'),i) - log(pu);
	}
	t = exp(t), d = exp(d); u = exp(u);
}

void TIS_T::scoring_XTG(const char* ssseq, double& t, double& d, double& u, int start)
{
	// set to zero
	int beg = 0;
	t = 0, d = 0, u = 0;

	int index = 0;
	int f = 0;
	// find the index of first (ORDER-1) nucleotide(s)
	for( ; f < ORDER ; ++f ){
		index += (int)pow(4.0,ORDER-f)*(ssseq[beg+f]-'0');
	}
	// get the score of the frist (ORDER-1) nucleotides(s)
	for( f = 0; f < 4; ++f ){
		t += exp(trueTIS(index+f, start)); 
		d += exp(downFalse(index+f, start));
		u += exp(upFalse(index+f, start));
	}
	// multipe the prior probabilities, as plus for log
	t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);

	// formulating the ORDER-th Markov model
	int i = 0;
	for( ; i+start < trueTIS.getColum(); ++i){
		int index = 0;
		int f = 0;
		for( ; f < ORDER ; ++f ){
			index += (int)pow(4.0,ORDER-f)*(ssseq[beg+i+f]-'0');
		}
		double pt = 0, pt2 = 0, pd = 0, pu = 0;
		for( f = 0; f < 4; ++f ){
			pt += exp(trueTIS(index+f, i+start)); 
			pd += exp(downFalse(index+f, i+start));
			pu += exp(upFalse(index+f, i+start));
		}
		t += trueTIS(index+(ssseq[beg+i+ORDER]-'0'), i+start) - log(pt);
		d += downFalse(index+(ssseq[beg+i+ORDER]-'0'), i+start) - log(pd);
		u += upFalse(index+(ssseq[beg+i+ORDER]-'0'), i+start) - log(pu);
	}
	t = exp(t), d = exp(d); u = exp(u);
}

void TIS_T::revise_head_enough(std::vector<ORF_T>& ORFSet, int head_enough_number)
{
	for(int num=0; num < maximalInterationNum ; num++)
	{
		initiatePWMs(ORFSet);
		std::vector< ORF_T >::iterator iter = ORFSet.begin();
		for (; iter != ORFSet.end(); iter++)
		{
			const char* seq = iter->codingSeq->data();
			int seqLen = iter->seqLen;

			if(!iter->enough_upbp_1stXTG)
				continue;

			iter->origLocation = iter->location;
			
			double maxScore = -1;
			iter->location.first = iter->leftmost;//每次都要扩到最长ORF再预测
			int hint = iter->location.first ;
			int indexc = 0;
		
			for( ; hint < iter->location.second -30 && indexc < maxCandidateN; hint += 3 )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end() && hint-UPBP >= 0 && hint+3+DOWNBP+ORDER <= seqLen )
				{
					++indexc;
					int beg = hint - UPBP;
					double t = 0, d = 0, u = 0;

					string ssseq = string(seq).substr(beg, UPBP+3+DOWNBP+ORDER);
					scoring_XTG(ssseq.data(), t, d, u);

					double score = (t)/(t+u+d);
					if( score > maxScore )
					{
						maxScore = score;
						iter->trueScore = score;
						iter->upFalseScore = (u)/(t+u+d);
						iter->downFalseScore = (d)/(t+u+d);
						iter->location.first = hint;
						iter->ORFLength =  iter->location.second - iter->location.first ;
					}
				}
			}
		}

		double unchanged=0;
		for (iter = ORFSet.begin(); iter != ORFSet.end(); iter++)
		{
			if(iter->enough_upbp_1stXTG && iter->origLocation == iter->location)
				++ unchanged;		
		}

		int size=(int)ORFSet.size();
		double percent=unchanged/head_enough_number*100;
//		cout << "\t\tunchanged percent:\t" << percent << " <" << unchanged << '/' << head_enough_number << '>' << endl;
		if( percent > 99.99 ){
			break;
		}
	}
}

void TIS_T::revise_not_enough( std::vector<ORF_T>& ORFSet )
{
	std::vector< ORF_T >::iterator iter = ORFSet.begin();
	for ( ; iter != ORFSet.end(); iter++)
	{
		const char* seq = iter->codingSeq->data();
		int seqLen = iter->seqLen;

		// if the first XTG upstream sequence is enough for scoring, then scoring them to find the one
		// with the highest true score
		if(iter->location.first >= minimum_upstream) {
			iter->origLocation = iter->location;

			double maxScore = -1;
			int hint = iter->location.first ;
			int indexc = 0;
			if(maxCandidateN<14)
				maxCandidateN=14;
			for( ; hint < iter->location.second -30 && indexc < maxCandidateN; hint += 3 )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end() && /*hint-UPBP >= 0 && */hint+3+DOWNBP+ORDER <= seqLen )
				{
					int beg = hint - UPBP;
					double t = 0, d = 0, u = 0;
					if(beg >= 0) {
						string ssseq = string(seq).substr(beg, UPBP+3+DOWNBP+ORDER);
						scoring_XTG(ssseq.data(), t, d, u);
					} else {
						string ssseq = string(seq).substr(0, UPBP+beg+3+DOWNBP+ORDER);
						scoring_XTG(ssseq.data(), t, d, u, (0-beg));
					}
					double score = (t)/(t+u+d);
					if( score > maxScore )
					{
						maxScore = iter->trueScore = score;
						iter->upFalseScore = (u)/(t+u+d);
						iter->downFalseScore = (d)/(t+u+d);
						iter->location.first = hint;
					}
				}
			}
		} else {	// else: find out if the first XTG is in the coding area
			//按location.first向下走，找第一个minimum_upstream后的同相位的XTG
			int hint = iter->location.first - 3;
			int initTIS = -1;
			for( ; hint < iter->location.second; hint += 3 )
			{
				if(hint >= minimum_upstream) 
				{
					int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
					if(STARTS.find(subStr) != STARTS.end())
					{
						initTIS = hint;
						break;
					}
				}
			}

			//如果initTIS仍为-1 即未找到同相位XTG 
			if ( initTIS == -1 )
			{
				iter->haveTIS = -1;
				continue;
			}
			
			double score = -1;
			int beg = initTIS - UPBP;
			// if upstream sequence is enough for scoring, use all PWM
			if( beg >= 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
			{
				double t = 0, d = 0, u = 0;
				string ssseq = string(seq).substr(beg, UPBP+3+DOWNBP+ORDER);
				scoring_XTG(ssseq.data(), t, d, u);
				score = d/(t+d+u);
			}
			// else, if upstream sequence is not enough, use as much as the PWM which is larger than minimum_upstream
			if( beg < 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
			{
				double t = 0, d = 0, u = 0;
				string ssseq = string(seq).substr(0, UPBP+beg+3+DOWNBP+ORDER);
				scoring_XTG(ssseq.data(), t, d, u, (0-beg));
				score = d/(t+d+u);
			}
 
			//大于cutoff 判定这个XTG属于coding区 此类ORF的TIS不在此片段中
			//可能在UPBP内还有XTG 但没有足够信息量 无法判断
			if (score >= cutoff) 
				continue;
			//小于cutoff 认为第一个XTG在noncoding区 则其下游的XTG中得分最高的判定为TIS
			else
			{
				double maxScore = -1;
			
				for( ; hint < iter->location.second -30; hint += 3)
				{
					int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
					if( STARTS.find(subStr) != STARTS.end() && hint-UPBP >= 0 && hint+3+DOWNBP+1+ORDER <= seqLen )
					{
						int beg = hint - UPBP;
						double t = 0, d = 0, u = 0;
						if(beg>=0) {
							string ssseq = string(seq).substr(beg, UPBP+3+DOWNBP+ORDER);
							scoring_XTG(ssseq.data(), t, d, u);
						} else {
							string ssseq = string(seq).substr(0, UPBP+beg+3+DOWNBP+ORDER);
							scoring_XTG(ssseq.data(), t, d, u, (0-beg));
						}
						double score = t/(t+d+u);
						if( score > maxScore )
						{
							maxScore = iter->trueScore = score;
							iter->upFalseScore = (u)/(t+u+d);
							iter->downFalseScore = (d)/(t+u+d);
							iter->location.first = hint;
						}
					}
				}
			}
		}
	}
}

void TIS_T::reviseTIS_HaveHead_NotEnough( std::vector<ORF_T>& ORFSet )
{
	//ORDER = MAXORDER;

	std::vector< ORF_T >::iterator iter = ORFSet.begin();
	for ( ; iter != ORFSet.end(); iter++)
	{
		const char* seq = iter->codingSeq->data();
		int seqLen = iter->seqLen;

		//如果location.first>=3，选取最高分
		if ( iter->location.first >= 3 && iter->location.first >= UPBP)
		{
			iter->origLocation = iter->location;

			if( iter->haveTIS == -1 )
				continue;

			double maxScore = -1;
			int hint = iter->location.first ;
			int indexc = 0;

				
			for( ; hint < iter->location.second -30 && indexc < maxCandidateN; hint += 3 )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end() && hint-UPBP >= 0 && hint+3+DOWNBP+ORDER <= seqLen )
				{
					int beg = hint - UPBP;
					double t = 0, d = 0, u = 0;
					int index = 0;
					int f = 0;
					for( ; f < ORDER ; ++f ){
						index += (int)pow(4.0,ORDER-f)*(seq[beg+f]-'0');
					}
					for( f = 0; f < 4; ++f ){
						t += exp(trueTIS(index+f,0)); 
						d += exp(downFalse(index+f,0));
						u += exp(upFalse(index+f,0));
					}

					t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
					int i = 0;
					for( ; i < trueTIS.getColum(); ++i ){
						int index = 0;
						int f = 0;
						for( ; f < ORDER ; ++f ){
							index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
						}
						double pt = 0, pt2 = 0, pd = 0, pu = 0;
						for( f = 0; f < 4; ++f ){
							pt += exp(trueTIS(index+f,i)); 
							pd += exp(downFalse(index+f,i));
							pu += exp(upFalse(index+f,i));
						}
						t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
						d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
						u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
					}
					t = exp(t), d = exp(d); u = exp(u);
					double score = (t)/(t+u+d);
					if( score > maxScore )
					{
						maxScore = iter->trueScore = score;
						iter->upFalseScore = (u)/(t+u+d);
						iter->downFalseScore = (d)/(t+u+d);
						iter->location.first = hint;
					}
				}
			}
		}
		//location.first<3 用cutoff
		else
		{
			//按location.first向下走，找第一个UPBP后的同相位的XTG
			int hint = iter->location.first ;
			int initTIS = -1;
			for( ; hint < iter->location.second; hint += 3 )
			{
				if ( hint >= UPBP )
				{
					int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
					if(STARTS.find(subStr) != STARTS.end())
					{
						initTIS = hint;
						break;
					}
				}
			}

			//如果initTIS仍为-1 即未找到同相位XTG 
			if ( initTIS == -1 )
			{
				iter->haveTIS = -1;
				continue;
			}
			
			double score = -1;
			int beg = initTIS - UPBP;
			if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
			{
				int beg = hint - UPBP;
				double t = 0, d = 0, u = 0;
				int index = 0;
				int f = 0;
				for( ; f < ORDER ; ++f ){
					index += (int)pow(4.0,ORDER-f)*(seq[beg+f]-'0');
				}
				for( f = 0; f < 4; ++f ){
					t += exp(trueTIS(index+f,0)); 
					d += exp(downFalse(index+f,0));
					u += exp(upFalse(index+f,0));
				}

				t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
				int i = 0;
				for( ; i < trueTIS.getColum(); ++i ){
					int index = 0;
					int f = 0;
					for( ; f < ORDER ; ++f ){
						index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
					}
					double pt = 0, pt2 = 0, pd = 0, pu = 0;
					for( f = 0; f < 4; ++f ){
						pt += exp(trueTIS(index+f,i)); 
						pd += exp(downFalse(index+f,i));
						pu += exp(upFalse(index+f,i));
					}
					t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
					d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
					u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
				}
				t = exp(t), d = exp(d); u = exp(u);
				score = d/(t+d+u);
			}
 
			//大于cutoff 判定这个XTG属于coding区 此类ORF的TIS不在此片段中
			//可能在UPBP内还有XTG 但没有足够信息量 无法判断
//			cutoff = 0;
			if (score >= cutoff) 
				continue;
			//小于cutoff 认为第一个XTG在noncoding区 则其下游的XTG中得分最高的判定为TIS
			else
			{
				double maxScore = -1;
			
				for( ; hint < iter->location.second -30; hint += 3)
				{
					int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
					if( STARTS.find(subStr) != STARTS.end() && hint-UPBP >= 0 && hint+3+DOWNBP+1+ORDER <= seqLen )
					{
						int beg = hint - UPBP;
						double t = 0, d = 0, u = 0;
						int index = 0;
						int f = 0;
						for( ; f < ORDER ; ++f ){
							index += (int)pow(4.0,ORDER-f)*(seq[beg+f]-'0');
						}
						for( f = 0; f < 4; ++f ){
							t += exp(trueTIS(index+f,0)); 
							d += exp(downFalse(index+f,0));
							u += exp(upFalse(index+f,0));
						}

						t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
						int i = 0;
						for( ; i < trueTIS.getColum(); ++i ){
							int index = 0;
							int f = 0;
							for( ; f < ORDER ; ++f ){
								index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
							}
							double pt = 0, pt2 = 0, pd = 0, pu = 0;
							for( f = 0; f < 4; ++f ){
								pt += exp(trueTIS(index+f,i)); 
								pd += exp(downFalse(index+f,i));
								pu += exp(upFalse(index+f,i));
							}
							t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
							d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
							u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
						}
						t = exp(t), d = exp(d); u = exp(u);
						double score = t/(t+d+u);
						if( score > maxScore )
						{
							maxScore = iter->trueScore = score;
							iter->upFalseScore = (u)/(t+u+d);
							iter->downFalseScore = (d)/(t+u+d);
							iter->location.first = hint;
						}
					}
				}
			}
		}
	}
}

void TIS_T::getcutoff( std::vector<ORF_T>& ORFSet )
{
	Set_Str downSet;
	std::list<double> downFalseScores;

	std::vector< ORF_T >::iterator iter = ORFSet.begin();
	for( ; iter != ORFSet.end(); ++iter )
	{
		const char* seq = iter->codingSeq->data();
		int seqLen = iter->seqLen;

		//如果location.first小于三，则该ORF几乎是没有头的ORF，真实TIS在该片段中的概率很小
		//不用此类ORF学习PWM，消除噪音影响
		if ( iter->location.first < 3 )
			continue;

		Str nter;
		//down false TIS学习coding区的得分分布情况
		int hint = iter->location.first + 3;
		for( ; hint < iter->location.second; hint += 3 )
		{
			int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
			if(STARTS.find(subStr) != STARTS.end())
			{
				int beg = hint - UPBP;
				if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1  + ORDER < seqLen  )
				{
					nter = Str(seq+beg,UPBP + 3 + DOWNBP) ;
					if ( downSet.find(nter)==downSet.end() )
					{
						int beg = hint - UPBP;
						double t = 0, d = 0, u = 0;
						int index = 0;
						int f = 0;
						for( ; f < ORDER ; ++f ){
							index += (int)pow(4.0,ORDER-f)*(seq[beg+f]-'0');
						}
						for( f = 0; f < 4; ++f ){
							t += exp(trueTIS(index+f,0)); 
							d += exp(downFalse(index+f,0));
							u += exp(upFalse(index+f,0));
						}

						t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
						int i = 0;
						for( ; i < trueTIS.getColum(); ++i ){
							int index = 0;
							int f = 0;
							for( ; f < ORDER ; ++f ){
								index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
							}
							double pt = 0, pt2 = 0, pd = 0, pu = 0;
							for( f = 0; f < 4; ++f ){
								pt += exp(trueTIS(index+f,i)); 
								pd += exp(downFalse(index+f,i));
								pu += exp(upFalse(index+f,i));
							}
							t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
							d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
							u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
						}
						t = exp(t), d = exp(d); u = exp(u);
						double score = d/(t+d+u);
						
						downSet.insert(nter);
						downFalseScores.push_back(score);

					}
				}
			}
		}
	}
	downFalseScores.sort();
	int size = (int)downFalseScores.size();
	std::list<double>::iterator iter2 = downFalseScores.begin();
	int i=1;
	for ( ; iter2 != downFalseScores.end(); iter2++ )
	{
		if ( double(i)/size > 0.05 )
		{
			this->cutoff = *iter2;
			break;
		}
		i++;
	}
}

std::string TIS_T::generate_seq(const std::vector<double>& gc, int length)
{ 
	string seq;
	srand(unsigned(time(0)));
	for(int i=0; i<length; ++i) {
		double rr = rand()/(RAND_MAX+1.0);
		double sum = 0;
		for(int f=0; f<(int)gc.size(); ++f) {
			sum += gc[f];
			if(sum-gc[f]<rr && rr<=sum) {
				seq += (f+48);
				break;
			}
		}
	}
	return seq;
}

std::vector<double> TIS_T::gc_content(const char* seq)
{
	int size = (int)string(seq).size();
	vector<double> gc(4, 0.001);
	double sum = 0.004;
	for(int i=0; i<size; ++i) {
		++ sum;
		++ gc[seq[i]-48];
	}
	for(int f=0; f<(int)gc.size(); ++f)
		gc[f] /= sum;
	return gc;
}

void TIS_T::revise_nohead(std::vector<ORF_T>& ORFSet)
{
//	int minimum_upstream = 8;
	int number = 0;

	std::vector< ORF_T >::iterator iter = ORFSet.begin();
	for( ; iter != ORFSet.end(); ++iter )
	{
		const char* seq = iter->codingSeq->data();
		int seqLen = iter->seqLen;

		// if have head and enough upstream sequences, continue
		if(iter->enough_upbp_1stXTG)
			continue;

		//按location.first向下走，找第一个UPBP后的同相位的XTG
		int hint = iter->location.first-3;	// ycliu: find XTG from the beginning of the sequence
		int initTIS = -1;
		for( ; hint < iter->location.second; hint += 3 )
		{
			if(hint >= minimum_upstream)	// ycliu: should be larger than the minimum upstream, not the UPBP
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if(STARTS.find(subStr) != STARTS.end())
				{
					initTIS = hint;
					break;
				}
			}
		}

		//如果initTIS仍为-1 即未找到同相位XTG 
		if ( initTIS == -1 )
		{
			iter->haveTIS = -1;
			continue;
		}
		
		double score = -1;
		int beg = initTIS - UPBP;
		// ycliu: beg>=0 means the upstream sequence is larger than UPBP, so use the sequence to calculate the PWM score
		if( beg >= 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
		{
			string ssseq = string(seq).substr(beg, UPBP+3+DOWNBP/*+1*/+ORDER);
			double t = 0, d = 0, u = 0;
			scoring_XTG(ssseq.data(), t, d, u);			
			score = d/(t+d+u);
		}

		// ycliu: beg<0 means upstream sequences is shorter than UPBP, so use as much sequence as possible 
		//        and then use background gc content to replace the exceeding ones
		if(beg < 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen)
		{
			int length = (0-beg);
			std::string ssseq = string(seq).substr(0, UPBP-length+3+DOWNBP+ORDER);
			double t = 0, d = 0, u = 0;
			scoring_XTG(ssseq.data(), t, d, u, length);
			score = d/(t+d+u);
		}

		//大于cutoff 判定这个XTG属于coding区 此类ORF的TIS不在此片段中
		//可能在UPBP内还有XTG 但没有足够信息量 无法判断

		if (score >= cutoff) 
			continue;

		//小于cutoff 认为第一个XTG在noncoding区 则其下游的XTG中得分最高的判定为TIS
		else
		{
			++ number;
			double maxScore = -1;
		
			for( ; hint < iter->location.second -60; hint += 3 )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end() && hint+3+DOWNBP+1+ORDER <= seqLen )
				{
					int beg = hint - UPBP;
					string ssseq;
					double t = 0, d = 0, u = 0;
					if(beg >= 0) {
						ssseq = string(seq).substr(beg, UPBP+3+DOWNBP+ORDER);
						TIS_T::scoring_XTG(ssseq.data(), t, d, u);
					} else {
						ssseq = string(seq).substr(0, UPBP+beg+3+DOWNBP+ORDER);
						TIS_T::scoring_XTG(ssseq.data(), t, d, u, (0-beg));
					}
					double score = t/(t+d+u);
					if( score > maxScore )
					{
						maxScore = iter->trueScore = score;
						iter->upFalseScore = (u)/(t+u+d);
						iter->downFalseScore = (d)/(t+u+d);
						iter->location.first = hint;
					}
				}
			}
		}
	}
}

void TIS_T::reviseTIS_NoHead( std::vector<ORF_T>& ORFSet )
{
	//this->cutoff = 0.9;

	int minimum_upstream = 50;
	int number = 0;

	std::vector< ORF_T >::iterator iter = ORFSet.begin();
	for( ; iter != ORFSet.end(); ++iter )
	{
		const char* seq = iter->codingSeq->data();
		int seqLen = iter->seqLen;

		//对没有头的ORF打分
		if ( iter->location.first >= 3 && iter->location.first >= UPBP)
			continue;
		++ number;

		// calculate the GC content of the sequence
		vector<double> gc(4,0.01);
		double sssize = 0;
		for(int ss=0; ss<(int)iter->codingSeq->size(); ++ss) {
			++ sssize;
			++ gc[iter->codingSeq->at(ss)-48];
		}
		for(int gg=0; gg<(int)gc.size(); ++gg) {
			gc[gg] /= sssize;
		}

		//按location.first向下走，找第一个UPBP后的同相位的XTG
//		int hint = iter->location.first;
		int hint = iter->location.first-3;	// ycliu: find XTG from the beginning of the sequence
		int initTIS = -1;
		for( ; hint < iter->location.second; hint += 3 )
		{
//			if ( hint >= UPBP )
			if(hint >= minimum_upstream)	// ycliu: should be larger than the minimum upstream, not the UPBP
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if(STARTS.find(subStr) != STARTS.end())
				{
					initTIS = hint;
					break;
				}
			}
		}

		//如果initTIS仍为-1 即未找到同相位XTG 
		if ( initTIS == -1 )
		{
			iter->haveTIS = -1;
			continue;
		}
		
		double score = -1;
		int beg = initTIS - UPBP;
		// ycliu: beg>=0 means the upstream sequence is larger than UPBP, so use the sequence to calculate the PWM score
		if( beg >= 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
		{
//			int beg = hint - UPBP;
			double t = 0, d = 0, u = 0;
			int index = 0;
			int f = 0;
			for( ; f < ORDER ; ++f ){
				index += (int)pow(4.0,ORDER-f)*(seq[beg+f]-'0');
			}
			for( f = 0; f < 4; ++f ){
				t += exp(trueTIS(index+f,0)); 
				d += exp(downFalse(index+f,0));
				u += exp(upFalse(index+f,0));
			}

			t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
			int i = 0;
			for( ; i < trueTIS.getColum(); ++i ){
				int index = 0;
				int f = 0;
				for( ; f < ORDER ; ++f ){
					index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
				}
				double pt = 0, pt2 = 0, pd = 0, pu = 0;
				for( f = 0; f < 4; ++f ){
					pt += exp(trueTIS(index+f,i)); 
					pd += exp(downFalse(index+f,i));
					pu += exp(upFalse(index+f,i));
				}
				t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
				d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
				u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
			}
			t = exp(t), d = exp(d); u = exp(u);
			score = d/(t+d+u);
			//cout << score << endl;
		}

		// ycliu: beg<0 means upstream sequences is shorter than UPBP, so use as much sequence as possible 
		//        and then use background gc content to replace the exceeding ones
		if( beg < 0 && beg + UPBP + 3 + DOWNBP + 1 + ORDER < seqLen  )
		{
			continue;
//			int beg = hint - UPBP;
			cout << "beg: " << beg << endl;
			string ssseq;
			srand(unsigned(time(0)));
			int ssbeg = beg;
			for(; ssbeg<0; ++ssbeg) {
				char c = ' ';
				double rr = rand()/(RAND_MAX+1.0);
				if(rr>0 && rr<=gc[0])
					c = '0';
				else if(rr>gc[0] && rr<=(gc[0]+gc[1]))
					c = '1';
				else if(rr>(gc[0]+gc[1]) && rr<=(gc[0]+gc[1]+gc[2]))
					c = '2';
				else if(rr>(gc[0]+gc[1]+gc[2]) && rr<=1)
					c = '3';
				else 
					c= '$';
				ssseq += c;
			}
			cout << ssseq << endl;
			cout << "ssbeg: " << ssbeg << endl;
			for(; ssbeg<beg+UPBP+3+DOWNBP+1+ORDER; ++ssbeg)
				ssseq += seq[ssbeg];
			cout << ssseq << '\t' << ssseq.size() << '\t' << trueTIS.getColum() << endl;
			/*for(; ssbeg<beg+UPBP+3; ++ssbeg)
				ssseq += seq[ssbeg];
			cout << ssseq << '\t' << ssseq.size() << '\t' << trueTIS.getColum() << endl;*/
			
//			system("pause");
//			continue;

			// ycliu: calculating the start probabilities of true, downfalse and upfalse
			ssbeg = 0;
			double t = 0, d = 0, u = 0;
			int index = 0;
			int f = 0;
			for( ; f < ORDER ; ++f ){
				index += (int)pow(4.0,ORDER-f)*(ssseq[ssbeg+f]-'0');
			}
			for( f = 0; f < 4; ++f ){
				t += exp(trueTIS(index+f,0)); 
				d += exp(downFalse(index+f,0));
				u += exp(upFalse(index+f,0));
			}
			// ycliu: multiple the prior probabilities
			t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);

			int i = 0;
			for( ; i < trueTIS.getColum(); ++i ){
				int index = 0;
				int f = 0;
				for( ; f < ORDER ; ++f ){
					index += (int)pow(4.0,ORDER-f)*(ssseq[ssbeg+i+f]-'0');
				}
				double pt = 0, pt2 = 0, pd = 0, pu = 0;
				for( f = 0; f < 4; ++f ){
					pt += exp(trueTIS(index+f,i)); 
					pd += exp(downFalse(index+f,i));
					pu += exp(upFalse(index+f,i));
				}
				t += trueTIS(index+(ssseq[ssbeg+i+ORDER]-'0'),i) - log(pt);
				d += downFalse(index+(ssseq[ssbeg+i+ORDER]-'0'),i) - log(pd);
				u += upFalse(index+(ssseq[ssbeg+i+ORDER]-'0'),i) - log(pu);
			}
			t = exp(t), d = exp(d); u = exp(u);
			score = d/(t+d+u);
			//cout << score << endl;
		}

		//大于cutoff 判定这个XTG属于coding区 此类ORF的TIS不在此片段中
		//可能在UPBP内还有XTG 但没有足够信息量 无法判断

		if (score >= cutoff) 
			continue;
		//小于cutoff 认为第一个XTG在noncoding区 则其下游的XTG中得分最高的判定为TIS
		else
		{
			double maxScore = -1;
		
			for( ; hint < iter->location.second -30; hint += 3 )
			{
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end() && hint-UPBP >= 0 && hint+3+DOWNBP+1+ORDER <= seqLen )
				{
					int beg = hint - UPBP;
					double t = 0, d = 0, u = 0;
					int index = 0;
					int f = 0;
					for( ; f < ORDER ; ++f ){
						index += (int)pow(4.0,ORDER-f)*(seq[beg+f]-'0');
					}
					for( f = 0; f < 4; ++f ){
						t += exp(trueTIS(index+f,0)); 
						d += exp(downFalse(index+f,0));
						u += exp(upFalse(index+f,0));
					}

					t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
					int i = 0;
					for( ; i < trueTIS.getColum(); ++i ){
						int index = 0;
						int f = 0;
						for( ; f < ORDER ; ++f ){
							index += (int)pow(4.0,ORDER-f)*(seq[beg+i+f]-'0');
						}
						double pt = 0, pt2 = 0, pd = 0, pu = 0;
						for( f = 0; f < 4; ++f ){
							pt += exp(trueTIS(index+f,i)); 
							pd += exp(downFalse(index+f,i));
							pu += exp(upFalse(index+f,i));
						}
						t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
						d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
						u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
					}
					t = exp(t), d = exp(d); u = exp(u);
					double score = t/(t+d+u);
					if( score > maxScore )
					{
						maxScore = iter->trueScore = score;
						iter->upFalseScore = (u)/(t+u+d);
						iter->downFalseScore = (d)/(t+u+d);
						iter->location.first = hint;
					}
				}
			}
		}
	}
	/*cout << "No head: " << number << endl;
	system("pause");*/
}

void TIS_T::set_ordinary_starts()
{
	STARTS.clear();
	STARTS.insert(14); //ATG
	STARTS.insert(30); //CTG
	STARTS.insert(62); //TTG
	STARTS.insert(46); //GTG
}

void TIS_T::set_ordinary_stops()
{
	STOPS.clear();
	STOPS.insert(48); //TAA
	STOPS.insert(56); //TGA
	STOPS.insert(50); //TAG
}

void TIS_T::set_abnormal_stops()
{
	STOPS.clear();
	STOPS.insert(48); //TAA
	//STOPS.insert(56); //TGA
	STOPS.insert(50); //TAG
}

void TIS_T::new_reviseTIS( )
{
	set_ordinary_starts();

	std::map<std::string, std::vector<ORF_T> >::iterator iter = GenusORFSetMap.begin();	
	for( ; iter!=GenusORFSetMap.end(); ++iter)
	{
		std::vector<ORF_T>& ORFSet = iter->second;
		std::string genus = iter->first;

		int head_enough_number = 0 ;
		for( int i=0; i<(int)ORFSet.size(); i++ ) {
			if(ORFSet[i].enough_upbp_1stXTG)
				++ head_enough_number;
		}

		if ( genus == "Mycoplasma" || genus == "Acholeplasma" || genus == "Aster" 
			|| genus == "Onion" || genus == "Ureaplasma" ) {
				set_abnormal_stops();
		} else {
			set_ordinary_stops();
		}

		if ( head_enough_number > HaveHeadNumCutoff )
		{
			cout << '.';
//			cout << "Genus: " << setw(20) << iter->first << '\t' << "Have head No.: " << setw(10) << head_enough_number 
//				<< " > " << HaveHeadNumCutoff << endl;
//			cout<<"  CASCADE COMBINATION OF DIFFERENT MARKOV MODELS..." << endl;
			for ( ORDER=0; ORDER < MAXORDER+0.1; ORDER++)
			{
//				if( ORDER < 1 ) 
//					cout<<"    POST-PROCESS ORIGNAL INPUT TISs BY A 0-TH ORDER MARKOV MODEL:\n";
//				else 
//					cout<<"    POST-PROCESS OUTPUT TISs FROM A "<<(ORDER-1)<<"-TH ORDER MARKOV MODEL BY \n"
//						  "                                  A "<<ORDER<<"-TH ORDER MARKOV MODEL:\n";
				revise_head_enough(ORFSet, head_enough_number);
			}
			ORDER = ORDER-1; //revise的时候 ORDER多加了1 此处恢复到MAXORDER
			getcutoff(ORFSet);
			revise_nohead(ORFSet);
		}
		else
		{
			cout << '.';
//			cout << "Genus: " << setw(20) << iter->first << '\t' << "Have head No.: " << setw(10) << head_enough_number 
//				<< " < " << HaveHeadNumCutoff << endl;
			//读取先验参数			
			string os_sep;
			if(is_windows) {
				os_sep = "\\";
			} else {
				os_sep = "/";
			}

			std::string tismodel_file = prioriparameters + os_sep + iter->first + ".tismodel";
			read_prior_tismodel(tismodel_file);

			// ORDER=2; //只用2阶
			ORDER=MAXORDER;
			revise_not_enough(ORFSet);
			ORDER=0;
		}
	}
	cout << endl;
}

void TIS_T::resultToFile( Str& resultFile, Str& resultFormat )
{
	std::ofstream result(resultFile.data());
	if(!result.good()) {
		cout << "can't open " << resultFile << "." << endl;
		exit(-1);
	}

	Str tmp="";

	std::map< Str, std::vector<ORF_T> >::iterator iter = GenusORFSetMap.begin();

	if( resultFormat == "MED")
	{
		for( ; iter!=GenusORFSetMap.end(); iter++)
		{
			std::vector<ORF_T>& ORFSet = (iter->second);
			for( int i=0; i<(int)ORFSet.size(); i++ )
			{
				ORFSet[i].location2FileForm((int)ORFSet[i].codingSeq->size());
				if(ORFSet[i].id!=tmp)
					result << ">" << ORFSet[i].id << endl;
				result << setw(8) << ORFSet[i].location.first << setw(8);
				result << ORFSet[i].location.second << setw(4);
				if (ORFSet[i].isPositive)
					result << "+" << endl;
				else
					result << "-" << endl;				
				tmp = ORFSet[i].id;
			}
		}
	}
	else if ( resultFormat == "GFF")
	{
		for( ; iter!=GenusORFSetMap.end(); iter++)
		{
			std::vector<ORF_T>& ORFSet = (iter->second);
			for( int i=0; i<(int)ORFSet.size(); i++ )
			{
				ORFSet[i].location2FileForm((int)ORFSet[i].codingSeq->size());
				if(ORFSet[i].id!=tmp)
				{
					result << "##gff-version  3" << endl;
					result << "##MetaTISA" << endl;
					result << "##metagenome sequence name: " << ORFSet[i].id << endl;
				}
				result << ORFSet[i].id << "\t";
				result << "MetaTISA" << "\t" << "CDS";
				result << "\t" << ORFSet[i].location.first << "\t";
				result << ORFSet[i].location.second << "\t" << "." << "\t";
				if (ORFSet[i].isPositive)
					result << "+" << "\t" << "." << endl;
				else
					result << "-" << "\t" << "." << endl;
				tmp = ORFSet[i].id;
				
			}
		}
	} else {
		return;
	}
	result.close();
}
