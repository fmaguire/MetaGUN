
#include "meta.h"
#include "OftenUsedFun.h"

vector<MetaORF> find_pos_LORFs(const string& seq, const int add_size, const int step,
							   bool is_pos_strand, const int min_ORF_len)
{
	Codons cc;
	cc.set_ss_codons();

	vector<MetaORF> pos_ORFs;
	int size = (int)seq.size() - add_size;
	int stop = 0;
	for(; stop < size-3; ++stop) 
	{
		int c = ( ((seq[stop] - 48)<<4) + ((seq[stop+1] - 48)<<2) + seq[stop+2] - 48);
		if(cc.stops.find(c) == cc.stops.end())			// stop codon: TAA, TGA, TAG
			continue;

		int begin = 0, end = 0;
		bool is_noncoding = false;
		bool is_leftmost = false;
		bool is_stopin = false;

		int next_stop = stop + step;
		while(next_stop < size - 3) 
		{
			int c = ( ((seq[next_stop]-48)<<4) + ((seq[next_stop+1]-48)<<2) + seq[next_stop+2]-48 );
			if(cc.stops.find(c) != cc.stops.end())
			{
				is_stopin = true;
				break;
			}
			next_stop += step;
		}

		next_stop += step;
		if(stop >= add_size) 
		{
			int start = stop + step;
			while( start < next_stop - 3 ) {
				int c = ( ((seq[start]-48)<<4) + ((seq[start+1]-48)<<2) + seq[start+2]-48 );
				if(cc.starts.find(c) != cc.starts.end())		// start codon: ATG, GTG, TTG
					break;
				start += step;
			}
			if(start < next_stop - 3) {
				begin = start;
				is_leftmost = true;
			}
			else is_noncoding = true;
		}
		if(is_noncoding)
			continue;
		// extract ORF with 3*N length
		int phase = stop % step;
		if(begin < add_size)
			begin = phase;
		else
			begin -= add_size;
		end = next_stop;
		if(next_stop > size) {
			end -= step;
		}
		end -= add_size;
		int ORFSize = end - begin;
		if(ORFSize < min_ORF_len)
			continue;

		if(!is_pos_strand) {
			int seqSize = (int)seq.size() - 2*add_size;
			int tmp = begin;
			begin = seqSize - end;
			end = seqSize - tmp;
		}

		MetaORF orf;
		orf.begin = begin + 1, orf.end = end;
		orf.leftmost_xtg = (orf.ispos() ? orf.begin : orf.end);		// to file format
		orf.phase = phase;
		orf.is_leftmost = is_leftmost, orf.is_stopin = is_stopin;
		orf.strand = (is_pos_strand ? '+' : '-');
		orf.cut_type[0] = (orf.is_leftmost ? '1' : '0'), orf.cut_type[1] = (orf.is_stopin ? '1' : '0');
		
		pos_ORFs.push_back(orf);
	}
	return pos_ORFs;
}

vector<MetaORF> find_LORFs(const string& seq, const int add_size, const int step, const int min_ORF_len)
{
	vector<MetaORF> pos_LORFs = find_pos_LORFs(seq, add_size, step, true, min_ORF_len);
	string negseq = SeqOp::revs(seq);
	vector<MetaORF> neg_LORFs = find_pos_LORFs(negseq, add_size, step, false, min_ORF_len);
	pos_LORFs.insert(pos_LORFs.end(), neg_LORFs.begin(), neg_LORFs.end());
	return pos_LORFs;
}

map<string, vector<MetaORF> > get_metaORF(const string& locfile)
{
	ifstream fin(locfile.data());
	if(!fin) { cout << "can't open " << locfile << endl; exit(-1); }

	map<string, vector<MetaORF> > metaORFs;
	int genenum = 0;
	string line;
	getline(fin, line);
	while(!fin.eof())
	{
		if(line[0] == '>')
		{
			string id;
			if(line.find(' ')==std::string::npos)
				id = line.substr(1);
			else
				id = line.substr(1, line.find(' ')-1);
			getline(fin, line);
			while(line[0] != '>' && !fin.eof())
			{
				MetaORF c;
				istringstream is(line);
				is >> c.begin >> c.end >> c.strand;
				map<string, vector<MetaORF> >::iterator iter = metaORFs.find(id);
				++ genenum;
				if(iter == metaORFs.end())
					metaORFs[id].push_back(c);
				else
					iter->second.push_back(c);

				getline(fin, line);
			}
		}
		else
			getline(fin, line);
	}
	fin.close();
	return metaORFs;
}

XTGScore scoring_xtg(const string& seq, int xtg_site, const TIS& tis)
{
	if(xtg_site - tis.upbps < 0 || xtg_site + 3 + tis.downbps + tis.this_order > (int)seq.size() - 1)
	{
		cout << xtg_site << '\t' << xtg_site + 3 + tis.downbps + tis.this_order << '\t' << seq.size() << endl;
		cout << "given sequence has not enough information for scoring." << endl;
		exit(-1);
	}

	double ts = 0, ds = 0, us = 0;
	int beg = xtg_site - tis.upbps;
	int index = 0;
	for(int f = 0; f < tis.this_order; ++ f)
		index += (int)pow(4.0, tis.this_order-f) * (seq[beg+f] - '0');
	// get initiation probability
	for(int f = 0; f < 4; ++ f)
	{
		ts += exp(tis.trueTIS(index+f, 0));
		ds += exp(tis.downfalse(index+f, 0));
		us += exp(tis.upfalse(index+f, 0));
	}

	ts = log(ts) + log(tis.prob_true);
	ds = log(ds) + log(tis.prob_downfalse);
	us = log(us) + log(tis.prob_upfalse);

	for(int i = 0; i < tis.trueTIS.getColum(); ++ i)
	{
		int index = 0;
		for(int f = 0; f < tis.this_order; ++ f)
			index += ((int)pow(4.0, tis.this_order-f) * (seq[beg+i+f] - '0'));
		double pt = 0, pd = 0, pu = 0;
		for(int f = 0; f < 4; ++ f)
		{
			pt += exp(tis.trueTIS(index+f, i));
			pd += exp(tis.downfalse(index+f, i));
			pu += exp(tis.upfalse(index+f, i));
		}

		int this_codon = index + (seq[beg+i+tis.this_order] - '0');
		ts += (tis.trueTIS(this_codon, i) - log(pt));
		ds += (tis.downfalse(this_codon, i) - log(pd));
		us += (tis.upfalse(this_codon, i) - log(pu));
	}

	ts = exp(ts), ds = exp(ds), us = exp(us);
	double total_prob = ts + ds + us;
	ts /= total_prob, ds /= total_prob, us /= total_prob;
	XTGScore xtg_s;
	xtg_s.true_score = ts, xtg_s.downfalse_score = ds, xtg_s.upfalse_score = us;
	return xtg_s;
}

// // // -------------- functions of MetaSeqs ---------------------// // //
void MetaSeqs::get_metaseqs(const std::string &seqfile)
{
	cout << "get sequences ... " << endl;
	ifstream fin(seqfile.data());
	if(!fin) { cout << "can't open " << seqfile << endl; exit(-1); }

	string line;
	getline(fin, line);
	while(!fin.eof())
	{
		if(line[0] == '>')
		{
			Read rd;
			rd.descript = line.substr(1);
			getline(fin, line);
			while(line[0] != '>' && !fin.eof())
			{
				rd.seq += line;
				getline(fin, line);
			}
			rd.rid = rd.descript.substr(0, rd.descript.find(' '));
			metaseqs[rd.rid] = rd;
		}
		else
			getline(fin, line);
	}
	fin.close();
	cout << '\t' << metaseqs.size() << " sequences read in." << endl;
}

void MetaSeqs::find_all_LORFs()
{
	cout << "find all longest ORFs ... " << endl;

	int step = 3;
	string addseq_head = "TAATTGATTAGT";
	string addseq_tail = "ACTAATCAATTA";
	int add_size = (int)addseq_head.size();
	
	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string posseq = SeqOp::Letter2Digital(addseq_head + iterread->second.seq + addseq_tail);
		vector<MetaORF> LORFs = find_LORFs(posseq, add_size, step, min_ORF_len);

		vector<MetaORF>::iterator iterorf = LORFs.begin();
		for(; iterorf != LORFs.end(); ++ iterorf)
		{
			string type = "00";
			type[0] = ((*iterorf).is_leftmost ? '1' : '0'), type[1] = ((*iterorf).is_stopin ? '1' : '0');
			iterread->second.ORFs[type].push_back(*iterorf);
			++ cut_type_num[type];
			++ all_ORF_num;
		}
	}

	cout << '\t' << all_ORF_num << " longest ORFs found." << endl;
}

void MetaSeqs::find_all_ORFsets()
{
	all_ORF_num = 0;
	cut_type_num["00"] = cut_type_num["01"] = cut_type_num["10"] = cut_type_num["11"] = 0;
	Codons cc;
	cc.set_ss_codons();

	cout << "find all ORFsets ... " << endl;
	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string readseq = SeqOp::Letter2Digital(iterread->second.seq);
		string neg_readseq = SeqOp::revs(readseq);
		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		iterread->second.ORFsets.erase(iterread->second.ORFsets.begin(), iterread->second.ORFsets.end());
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			for(int i = 0; i < (int)itertype->second.size(); ++ i)
			{
				MetaORF mo = itertype->second[i];
				iterread->second.ORFsets[itertype->first].push_back(mo);
				++ all_ORF_num;
				++ cut_type_num[itertype->first];
				if(mo.ispos())
				{
					for(int pos = mo.begin - 1; pos < mo.end - min_ORF_len; pos += 3)
					{
						int codon = ((readseq[pos]-'0')<<4) + ((readseq[pos+1]-'0')<<2) + (readseq[pos+2]-'0');
						if(cc.starts.find(codon) != cc.starts.end())
						{
							if(pos + 1 == mo.begin)
								continue;
							MetaORF mmo = mo;
							mmo.begin = pos + 1;
							string type = itertype->first;
							type[0] = '1';
							iterread->second.ORFsets[type].push_back(mmo);
							++ all_ORF_num;
							++ cut_type_num[type];
						}
					}
				}
				else
				{
					int start = (int)readseq.size() - mo.end + 1;
					int stop = (int)readseq.size() - mo.begin + 1;
					for(int pos = start - 1; pos < stop - min_ORF_len; pos += 3)
					{
						int codon = ((neg_readseq[pos]-'0')<<4) + ((neg_readseq[pos+1]-'0')<<2) + (neg_readseq[pos+2]-'0');
						if(cc.starts.find(codon) != cc.starts.end())
						{
							if((int)readseq.size() - pos == mo.end)
								continue;
							MetaORF mmo = mo;
							mmo.end = (int)readseq.size() - pos;
							string type = itertype->first;
							type[0] = '1';
							iterread->second.ORFsets[type].push_back(mmo);
							++ all_ORF_num;
							++ cut_type_num[type];
						}
					}
				}
			}
		}
	}
	cout << '\t' << all_ORF_num << " ORFs found in a ORFsets way. " << endl;
}

void MetaSeqs::label_ORFs(const string& coding_cds_file)
{
	cout << "label ORFs as coding and non-coding ORFs ... " << endl;
	map<string, vector<MetaORF> > coding_ORFs = get_metaORF(coding_cds_file);

	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		map<string, vector<MetaORF> >::iterator iterorfs = coding_ORFs.find(iterread->first);
		if(iterorfs == coding_ORFs.end())
		{
			map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
			for(; itertype != iterread->second.ORFs.end(); ++ itertype)
			{
				for(int a = 0; a < (int)itertype->second.size(); ++ a)
				{
					MetaORF& ao = itertype->second[a];
					ao.coding_type = "NC";
					++ noncod_ORF_num;
				}
			}
			continue;
		}

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			for(int a = 0; a < (int)itertype->second.size(); ++ a)
			{
				MetaORF& ao = itertype->second[a];
				bool iscoding = false;
				for(int c = 0; c < (int)iterorfs->second.size(); ++ c)
				{
					MetaORF& co = iterorfs->second[c];
					if(ao.stop() == co.stop())				// check
					{
						++ cod_ORF_num;
						if(ao.start() != co.start())
						{
							ao.adjust_start(co.start());
							++ reloc_ORF_num;
						}
						iscoding = true;
						ao.coding_type = "C";
						break;

					}
				}
				if(!iscoding)
				{
					ao.coding_type = "NC";
					++ noncod_ORF_num;
				}
			}
		}
		
	}
	cout << '\t' << cod_ORF_num << " coding ORFs, " << noncod_ORF_num << " non-coding ORFs, "
		<< reloc_ORF_num << " coding ORFs with starts relocated." << endl;
}

void MetaSeqs::label_ORFs_shadow_rule(const string& coding_cds_file)
{
	int overlap_threshold = 90;

	cout << "label ORFs as coding and non-coding in the shadow rule way ... " << endl;
	map<string, vector<MetaORF> > coding_ORFs = get_metaORF(coding_cds_file);

	map<string, Read >::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		map<string, vector<MetaORF> >::iterator iterorf = coding_ORFs.find(iterread->first);
		if(iterorf == coding_ORFs.end())
		{
			map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
			for(; itertype != iterread->second.ORFs.end(); ++ itertype)
				for(int ii = 0; ii < (int) itertype->second.size(); ++ ii)
					itertype->second[ii].coding_type = "OT";
			continue;
		}
		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			for(int aa = 0; aa < (int) itertype->second.size(); ++ aa)
			{
				MetaORF& ao = itertype->second[aa];
				ao.coding_type = "OT";

				for(int cc = 0; cc < (int)iterorf->second.size(); ++ cc)
				{
					MetaORF& co = iterorf->second[cc];
					if(ao.stop() == co.stop())
					{
						ao.coding_type = "C";
						++ cod_ORF_num;
						if(ao.start() != co.start())
						{
							++ reloc_ORF_num;
							ao.adjust_start(co.start());
						}
						break;
					}

					// label in the shadow rule way
					if(ao.coding_type == "NC")
						continue;
					int overlap = 0 - (std::max(ao.begin, co.begin) - std::min(ao.end, co.end));
					if(overlap > overlap_threshold)
					{
						ao.coding_type = "NC";
						++ noncod_ORF_num;
					}
				}
			}
		}
	}
}

void MetaSeqs::get_bin_info(const std::string &binfile)
{
	vector<string> lines = readlines(binfile);
	readbin.erase(readbin.begin(), readbin.end());
	for(int i=0; i<(int)lines.size(); ++i)
	{
		vector<string> vv = parse_string(lines[i], "\t");
		readbin[vv[0].substr(1, vv[0].find(' ')-1)] = std::make_pair(vv[1], vv[2]);
	}

	genus_binreads.erase(genus_binreads.begin(), genus_binreads.end());
	group_binreads.erase(group_binreads.begin(), group_binreads.end());
	map<string, pair<string, string> >::iterator iterread = readbin.begin();
	for(; iterread != readbin.end(); ++ iterread)
	{
		genus_binreads[iterread->second.first].push_back(iterread->first);
		group_binreads[iterread->second.second].push_back(iterread->first);
	}
	cout << '\t' << readbin.size() << " reads are classified into " << genus_binreads.size() 
		<< " genus and " << group_binreads.size() << " groups." << endl;
}

void displayMatrix( const Matrix_T< double >& ma, std::ostream& out )
{
	out << std::setw(10) << ma.getLine() << std::setw(10) << ma.getColum() << endl;
	for( int i=0; i<ma.getColum(); ++i )
	{
		for( int j=0; j<ma.getLine(); ++j )
			out << ma(j,i) << "\t";
		out << endl;
	}
}

void MetaSeqs::tis_scoring(const std::string &prior_para_path)
{
	cout << "get TIS scores ... " << endl;
	map<string, vector<string> >::iterator iter_genus = genus_binreads.begin();
	for(; iter_genus != genus_binreads.end(); ++ iter_genus)
	{
		TIS tis;
		tis.initiate_settings();
		tis.set_codons();

		string genus = iter_genus->first;
		tis.get_para_path(prior_para_path);
		string os_sep = (is_unix ? "/" : "\\");
		string tismodel_file = prior_para_path + os_sep + iter_genus->first + ".tismodel";
		tis.read_prior_tismodel(tismodel_file);

		for(int i = 0; i < (int)iter_genus->second.size(); ++ i)
		{

			string rid = iter_genus->second[i];
			map<string, Read>::iterator iterread = metaseqs.find(rid);
			if(iterread == metaseqs.end())
				continue;

			map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin(); 
			for(; itertype != iterread->second.ORFs.end(); ++ itertype)
			{
				if(itertype->first == "00" || itertype->first == "01")
					continue;
				for(int k = 0; k < (int)itertype->second.size(); ++ k)
				{
					MetaORF& mo = itertype->second[k];
					string readseq = SeqOp::Letter2Digital(iterread->second.seq);
					int seqsize = (int)readseq.size();
					int start = mo.begin, stop = mo.end;
					if(!mo.ispos())
					{
						readseq = SeqOp::revs(SeqOp::Letter2Digital(iterread->second.seq));
						start = seqsize - mo.end + 1;
						stop = seqsize - mo.begin + 1;
					}
 
					mo.true_ts = mo.upfalse_ts = mo.downfalse_ts = 1e-10;
					for(int hint = start - 1; hint <= stop - min_ORF_len; hint += 3)
					{
						int codon = ((readseq[hint] - '0')<<4) + ((readseq[hint+1] - '0')<<2) + (readseq[hint+2] - '0');
						if(tis.Starts.find(codon) != tis.Starts.end() && hint >= tis.minimum_upstream 
							&& hint + 3 + tis.downbps + tis.this_order <= seqsize)
						{
							int beg = hint - tis.upbps;
							double true_score=0, down_score=0, up_score = 0;
							if(beg >= 0 && beg+tis.upbps+3+tis.downbps+1+tis.this_order < seqsize) {
								string ssseq = readseq.substr(beg, tis.upbps+3+tis.downbps+tis.this_order);
								tis.scoring_XTG(ssseq.data(), true_score, down_score, up_score, 0);
							} 
							if(beg < 0 && beg+tis.upbps+3+tis.downbps+1+tis.this_order < seqsize) {
								int length = (0-beg);
								string ssseq = readseq.substr(0, tis.upbps-length+3+tis.downbps+tis.this_order);
								tis.scoring_XTG(ssseq.data(), true_score, down_score, up_score, length);
							}
							if(true_score > mo.true_ts)
							{
								mo.true_ts = true_score;
								mo.downfalse_ts = down_score;
								mo.upfalse_ts = up_score;
								mo.true_tis_site = (mo.ispos() ? (hint + 1) : (seqsize - hint));
							}
						}
					}
				}				
			}
		}
	}
}

void MetaSeqs::tis_scoring_ORFsets(const string& prior_para_path)
{
	cout << "get TIS scores in the ORFsets way ... " << endl;
	map<string, vector<string> >::iterator iter_genus = genus_binreads.begin();
	for(; iter_genus != genus_binreads.end(); ++ iter_genus)
	{
		TIS tis;
		tis.initiate_settings();
		tis.set_codons();

		/*tis.get_para_path(prior_para_path);
		string genus = iter_genus->first;
		tis.get_prior_paras(genus);*/

		string os_sep = (is_unix ? "/" : "\\");
		string tismodel_file = prior_para_path + os_sep + iter_genus->first + ".tismodel";
		tis.read_prior_tismodel(tismodel_file);

		for(int i = 0; i < (int)iter_genus->second.size(); ++ i)
		{

			string rid = iter_genus->second[i];
			map<string, Read>::iterator iterread = metaseqs.find(rid);
			if(iterread == metaseqs.end())
				continue;

			map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFsets.begin(); 
			for(; itertype != iterread->second.ORFsets.end(); ++ itertype)
			{
				if(itertype->first == "00" || itertype->first == "01")
					continue;
				for(int k = 0; k < (int)itertype->second.size(); ++ k)
				{
					MetaORF& mo = itertype->second[k];
					string readseq = SeqOp::Letter2Digital(iterread->second.seq);
					int seqsize = (int)readseq.size();
					int start = mo.begin, stop = mo.end;
					if(!mo.ispos())
					{
						readseq = SeqOp::revs(SeqOp::Letter2Digital(iterread->second.seq));
						start = seqsize - mo.end + 1;
						stop = seqsize - mo.begin + 1;
					}
					if(start - 1< tis.upbps || (start + tis.downbps) >= stop)
						continue;

					XTGScore xtg_s = scoring_xtg(readseq, start - 1, tis);
					mo.true_ts = xtg_s.true_score;
					mo.downfalse_ts = xtg_s.downfalse_score; 
					mo.upfalse_ts = xtg_s.upfalse_score;
					mo.true_tis_site = mo.start();
				}
			}
		}
	}
}

void MetaSeqs::get_allfea(const std::string &name, bool is_ORFsets)
{
	if(is_ORFsets)
		cout << "get feature files for svm predicting in the ORFsets way... " << endl;
	else
		cout << "get feature files for svm predicting ... " << endl;
	ofstream foutfea[2][2], foutorf[2][2];
	map<string, int> ORF_distri;

	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		string seq = iterread->second.seq;

		map<string, vector<MetaORF> >::iterator itertype, endtype;
		if(is_ORFsets)
			itertype = iterread->second.ORFsets.begin(), endtype = iterread->second.ORFsets.end();
		else
			itertype = iterread->second.ORFs.begin(), endtype = iterread->second.ORFs.end();
		for(; itertype != endtype; ++ itertype)
		{
			int head = atoi(itertype->first.substr(0,1).data());
			int tail = atoi(itertype->first.substr(1,1).data());

			if(!foutfea[head][tail].is_open())
			{
				string feafile(name + ".fea_" + itertype->first);
				foutfea[head][tail].open(feafile.data());
				if(!foutfea[head][tail]) { cout << "can't open " << feafile << endl; exit(-1); }
			}
			if(!foutorf[head][tail].is_open())
			{
				string orffile(name + ".orf_" + itertype->first);
				foutorf[head][tail].open(orffile.data());
				if(!foutorf[head][tail]) { cout << "can't open " << orffile << endl; exit(-1); }
			}

			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
					: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
				vector<double> fea = Feature::getfea(orfseq, max_ORF_len, mo.upfalse_ts, mo.true_ts, mo.downfalse_ts, mo.is_leftmost);
				Feature::foutfea(foutfea[head][tail], fea, mo.iscoding());
				int start = mo.begin;
				int stop = mo.end;
				if(mo.cut_type[0] == '1')
				{
					if(mo.true_tis_site != -1)
					{
						mo.ispos() ? start = mo.true_tis_site : 0;
						mo.ispos() ? 0 : stop = mo.true_tis_site;
					}
				}
				foutorf[head][tail] << rid << '|' << start << '|' << stop << '|' << mo.strand 
					<< '|' << mo.phase << '|' << mo.cut_type << endl;
				++ ORF_distri[itertype->first], ++ ORF_distri["all"];
			}
		}
	}

	for(int head = 0; head < 2; ++ head)
		for(int tail = 0; tail < 2; ++ tail)
			foutfea[head][tail].close(), foutorf[head][tail].close();

	map<string, int>::iterator iter = ORF_distri.begin();
	for(; iter != ORF_distri.end(); ++ iter)
		cout << '\t' << std::setw(4) << iter->first << ':' << std::setw(10) << iter->second << endl;
}

void MetaSeqs::get_allfea(const std::string &name, int level, bool is_ORFsets)
{
	if(level != 1 && level != 2)
	{
		cout << "Level setting error! 1 for genus, 2 for groups, other is invalid." << endl;
		exit(-1);
	}

	if(is_ORFsets)
		cout << "get feature files for svm predicting in the ORFsets way ... " << endl;
	else
		cout << "get feature files for svm predicting ... " << endl;
	map<string, int> ORF_distri;

	map<string, vector<string> >::iterator iterbin, enditer;
	if(level == 1) 
		iterbin = genus_binreads.begin(), enditer = genus_binreads.end();
	else 
		iterbin = group_binreads.begin(), enditer = group_binreads.end();

	for(; iterbin != enditer; ++ iterbin)
	{
		if(iterbin->second.empty())
			continue;

		ofstream foutfea[2][2], foutorf[2][2];
		for(int i = 0; i < (int)iterbin->second.size(); ++ i)
		{
			map<string, Read>::iterator iterread = metaseqs.find(iterbin->second[i]);
			if(iterread == metaseqs.end())
				continue;

			string rid = iterread->first;
			string seq = iterread->second.seq;

			map<string, vector<MetaORF> >::iterator itertype, endtype;
			if(is_ORFsets) 
				itertype = iterread->second.ORFsets.begin(), endtype = iterread->second.ORFsets.end();
			else
				itertype = iterread->second.ORFs.begin(), endtype = iterread->second.ORFs.end();
			for(; itertype != endtype; ++ itertype)
			{
				if(itertype->second.empty())
					continue;
				int head = atoi(itertype->first.substr(0,1).data());
				int tail = atoi(itertype->first.substr(1,1).data());

				if(!foutfea[head][tail].is_open())
				{
					string feafile(name + "." + iterbin->first + ".fea_" + itertype->first);
					foutfea[head][tail].open(feafile.data());
					if(!foutfea[head][tail]) { cout << "can't open " << feafile << endl; exit(-1); }
				}
				if(!foutorf[head][tail].is_open())
				{
					string orffile(name + "." + iterbin->first + ".orf_" + itertype->first);
					foutorf[head][tail].open(orffile.data());
					if(!foutorf[head][tail]) { cout << "can't open " << orffile << endl; exit(-1); }
				}

				vector<MetaORF>::iterator iterorf = itertype->second.begin();
				for(; iterorf != itertype->second.end(); ++ iterorf)
				{
					MetaORF mo = (*iterorf);
					string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
						: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
					vector<double> fea = Feature::getfea(orfseq, max_ORF_len, mo.upfalse_ts, mo.true_ts, mo.downfalse_ts, mo.is_leftmost);
					Feature::foutfea(foutfea[head][tail], fea, mo.iscoding());
					int start = mo.begin;
					int stop = mo.end;
					if(mo.cut_type[0] == '1')
					{
						if(mo.true_tis_site != -1)
						{
							mo.ispos() ? start = mo.true_tis_site : 0;
							mo.ispos() ? 0 : stop = mo.true_tis_site;
						}
					}
					foutorf[head][tail] << rid << '|' << start << '|' << stop << '|' << mo.strand 
						<< '|' << mo.phase << '|' << mo.cut_type << endl;
					++ ORF_distri[itertype->first], ++ ORF_distri["all"];
				}
			}
		}

		for(int head = 0; head < 2; ++ head)
			for(int tail = 0; tail < 2; ++ tail)
				foutfea[head][tail].close(), foutorf[head][tail].close(); 
	}

	map<string, int>::iterator iter = ORF_distri.begin();
	for(; iter != ORF_distri.end(); ++ iter)
		cout << '\t' << std::setw(4) << iter->first << ':' << std::setw(10) << iter->second << endl;
}

void MetaSeqs::get_trainfea(const std::string &name)
{
	cout << "get feature files for svm training ... " << endl;
	ofstream fout[2][2];

	map<string, pair<int, int> > ORF_TF_distri;
	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		string seq = iterread->second.seq;

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			int head = atoi(itertype->first.substr(0,1).data());
			int tail = atoi(itertype->first.substr(1,1).data());

			if(!fout[head][tail].is_open())
			{
				string feafile(name + ".trainfea_" + itertype->first);
				fout[head][tail].open(feafile.data());
				if(!fout[head][tail]) { cout << "can't open " << feafile << endl; exit(-1); }
			}

			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				if(mo.coding_type == "OT")
					continue;
				string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
					: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
				vector<double> fea = Feature::getfea(orfseq, max_ORF_len, mo.upfalse_ts, mo.true_ts, mo.downfalse_ts, mo.is_leftmost);
				Feature::foutfea(fout[head][tail], fea, mo.iscoding());
				mo.iscoding() ? ++ ORF_TF_distri[itertype->first].first : ++ ORF_TF_distri[itertype->first].second;
				mo.iscoding() ? ++ ORF_TF_distri["all"].first : ++ ORF_TF_distri["all"].second;
			}
		}
	}

	for(int head = 0; head < 2; ++ head)
		for(int tail = 0; tail < 2; ++ tail)
			fout[head][tail].close();

	map<string, pair<int, int> >::iterator iter = ORF_TF_distri.begin();
	cout << '\t' << std::setw(4) << "  " <<  setw(10) << "C" << setw(10) << "NC" 
		<< setw(10) << "total" << endl;
	for(; iter != ORF_TF_distri.end(); ++ iter)
		cout << '\t' << std::setw(4) << iter->first << ':' << std::setw(10) << iter->second.first 
		<< std::setw(10) << iter->second.second << setw(10) << iter->second.first + iter->second.second << endl;
}

void MetaSeqs::get_trainfea(const std::string &name, int level)
{
	if(level != 1 && level != 2)
	{
		cout << "Level setting error! 1 for genus, 2 for groups, other is invalid." << endl;
		exit(-1);
	}

	cout << "get feature files for svm training according to binning result ... " << endl;
	map<string, vector<string> >::iterator iterbin, enditer;
	if(level == 1) 
		iterbin = genus_binreads.begin(), enditer = genus_binreads.end();
	else 
		iterbin = group_binreads.begin(), enditer = group_binreads.end();

	for(; iterbin != enditer; ++ iterbin)
	{
		cout << iterbin->first << endl;
		if(iterbin->second.empty())
			continue;
		ofstream fout[2][2];
		map<string, pair<int, int> > ORF_TF_distri;
		for(int ii=0; ii<(int)iterbin->second.size(); ++ii)
		{
			map<string, Read>::iterator iterread = metaseqs.find(iterbin->second[ii]);
			if(iterread == metaseqs.end())
				continue;
			string rid = iterread->first;
			string seq = iterread->second.seq;
			map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
			for(; itertype != iterread->second.ORFs.end(); ++ itertype)
			{
				int head = atoi(itertype->first.substr(0,1).data());
				int tail = atoi(itertype->first.substr(1,1).data());

				if(!fout[head][tail].is_open())
				{
					string feafile(name + "." + iterbin->first + ".trainfea_" + itertype->first);
					fout[head][tail].open(feafile.data());
					if(!fout[head][tail]) { cout << "can't open " << feafile << endl; exit(-1); }
				}

				vector<MetaORF>::iterator iterorf = itertype->second.begin();
				for(; iterorf != itertype->second.end(); ++ iterorf)
				{
					MetaORF mo = (*iterorf);
					if(mo.coding_type == "OT")
						continue;
					string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
						: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
					vector<double> fea = Feature::getfea(orfseq, max_ORF_len, mo.upfalse_ts, mo.true_ts, mo.downfalse_ts, mo.is_leftmost);
					Feature::foutfea(fout[head][tail], fea, mo.iscoding());
					mo.iscoding() ? ++ ORF_TF_distri[itertype->first].first : ++ ORF_TF_distri[itertype->first].second;
					mo.iscoding() ? ++ ORF_TF_distri["all"].first : ++ ORF_TF_distri["all"].second;
				}
			}
		}
		for(int head = 0; head < 2; ++ head)
			for(int tail = 0; tail < 2; ++ tail)
				fout[head][tail].close();

		map<string, pair<int, int> >::iterator iter = ORF_TF_distri.begin();
		cout << '\t' << std::setw(4) << "  " <<  setw(10) << "C" << setw(10) << "NC" 
			<< setw(10) << "total" << endl;
		for(; iter != ORF_TF_distri.end(); ++ iter)
			cout << '\t' << std::setw(4) << iter->first << ':' << std::setw(10) << iter->second.first 
			<< std::setw(10) << iter->second.second << setw(10) << iter->second.first + iter->second.second << endl;
	}
}

//*
void MetaSeqs::get_trainfea_group(const std::string &name)
{
	cout << "get feature files for svm training ... " << endl;
	ofstream fout[2][2];

	map<string, pair<int, int> > ORF_TF_distri;
	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		string seq = iterread->second.seq;

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			int head = atoi(itertype->first.substr(0,1).data());
			int tail = atoi(itertype->first.substr(1,1).data());

			if(!fout[head][tail].is_open())
			{
				string feafile(name + ".trainfea_" + itertype->first);
				fout[head][tail].open(feafile.data());
				if(!fout[head][tail]) { cout << "can't open " << feafile << endl; exit(-1); }
			}

			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				if(mo.coding_type == "OT")
					continue;
				string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
					: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
				vector<double> fea = Feature::getfea(orfseq, max_ORF_len, mo.upfalse_ts, mo.true_ts, mo.downfalse_ts, mo.is_leftmost);
				Feature::foutfea(fout[head][tail], fea, mo.iscoding());
				mo.iscoding() ? ++ ORF_TF_distri[itertype->first].first : ++ ORF_TF_distri[itertype->first].second;
				mo.iscoding() ? ++ ORF_TF_distri["all"].first : ++ ORF_TF_distri["all"].second;
			}
		}
	}

	for(int head = 0; head < 2; ++ head)
		for(int tail = 0; tail < 2; ++ tail)
			fout[head][tail].close();

	map<string, pair<int, int> >::iterator iter = ORF_TF_distri.begin();
	cout << '\t' << std::setw(4) << "  " <<  setw(10) << "C" << setw(10) << "NC" 
		<< setw(10) << "total" << endl;
	for(; iter != ORF_TF_distri.end(); ++ iter)
		cout << '\t' << std::setw(4) << iter->first << ':' << std::setw(10) << iter->second.first 
		<< std::setw(10) << iter->second.second << setw(10) << iter->second.first + iter->second.second << endl;
}
//*/

void MetaSeqs::fout_ORF_seqs(const std::string &outfile, bool is_aa)
{
	ofstream foutseq(outfile.data());
	if(!foutseq) { cout << "can't open " << outfile << endl; exit(-1); }

	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		string seq = iterread->second.seq;

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
					: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
				orfseq = is_aa ? SeqOp::toAminoSeq(orfseq.data(), 0, (int)orfseq.size()) : SeqOp::Digital2Letter(orfseq);
				foutseq << '>' << rid << '|' << mo.begin << '|' << mo.end << '|' << mo.strand 
					<< '|' << mo.phase << '|' << mo.cut_type << endl;
				for(int k = 0; k <= (int)orfseq.size() / 70; ++k)
					foutseq << orfseq.substr(70*k, 70) << endl;
			}
		}
	}
	foutseq.close();
}

bool is_low_complexity(const string& seq, int repeat_thres = 10)
{
	int max_repeat_num = 0;
	int repeat_num = 1;
	char pre = ' ';
	for(int ii = 0; ii < (int)seq.size(); ++ ii)
	{
		char now = seq[ii];
		if(now == pre)
			++ repeat_num;
		else
			repeat_num = 1;
		pre = now;
		if(repeat_num > max_repeat_num)
			max_repeat_num = repeat_num;
	}
	if(max_repeat_num > repeat_thres)
		return true;
	else
		return false;
}

void MetaSeqs::fout_cod_seqs(const string& outfile, bool is_amino)
{
	ofstream foutseq(outfile.data());
	if(!foutseq) { cout << "can't open " << outfile << endl; exit(-1); }

	int coding_orfseq_num = 0;
	int low_complexity_num = 0, wrong_location_num = 0, short_length_num = 0;
	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		string seq = iterread->second.seq;

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				if(!mo.iscoding())
					continue;
				if(mo.start() < 0 || mo.start() > (int)seq.size() || mo.stop() < 0 || mo.stop() > (int)seq.size())
				{
					++ wrong_location_num;
					continue;
				}
				++ coding_orfseq_num;
				string orfseq = (mo.ispos() ? SeqOp::Letter2Digital(seq.substr(mo.begin-1, mo.size()))
					: SeqOp::Letter2Digital(SeqOp::revs(seq.substr(mo.begin-1, mo.size()))));
				if(is_amino)
					orfseq = SeqOp::toAminoSeq(orfseq.data(), 0, (int)orfseq.size());
				else 
					orfseq = SeqOp::Digital2Letter(orfseq);
				if(is_low_complexity(orfseq))
				{
					++ low_complexity_num;
					continue;
				}
				if((int)orfseq.size() < min_ORF_len / 3 - 1)
				{
					++ short_length_num;
					continue;
				}
//				orfseq = is_aa ? SeqOp::toAminoSeq(orfseq.data(), 0, (int)orfseq.size()) : SeqOp::Digital2Letter(orfseq);
				foutseq << '>' << rid << '|' << mo.begin << '|' << mo.end << '|' << mo.strand 
					<< '|' << mo.phase << '|' << mo.cut_type << endl;
				for(int k = 0; k <= (int)orfseq.size() / 70; ++k)
					foutseq << orfseq.substr(70*k, 70) << endl;
			}
		}
	}
	foutseq.close();
	cout << coding_orfseq_num << " coding sequences output to file " << outfile << endl;
	cout<< "wrong location: " << wrong_location_num << endl
		<< "short   length: " << short_length_num << endl
		<< "low complexity: " << low_complexity_num << endl;
}

void MetaSeqs::fout_ORF_locs(const std::string &outfile, bool is_coding_only)
{
	ofstream foutorf(outfile.data());
	if(!foutorf) { cout << "can't open " << outfile << endl; exit(-1); }

	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		foutorf << '>' << rid << endl;

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				if(is_coding_only)
				{
					if(!mo.iscoding())
						continue;
					foutorf << setw(8) << mo.begin << setw(8) << mo.end << setw(4) << mo.strand 
						<< setw(4) << mo.phase << setw(4) << itertype->first << endl;
				}
				else
					foutorf << setw(8) << mo.begin << setw(8) << mo.end << setw(4) << mo.strand 
					<< setw(4) << mo.phase << setw(4) << itertype->first << endl;
			}
		}
	}
	foutorf.close();
}

void MetaSeqs::fout_shd_ORF_locs(const string& outfile_basename)
{
	string shd_c = outfile_basename + ".shd_c.locs";
	string shd_nc = outfile_basename + ".shd_nc.locs";
	string shd_ot = outfile_basename + ".shd_ot.locs";

	ofstream fout_c(shd_c.data());
	if(!fout_c) { cout << "can't open " << shd_c << endl; exit(-1); }
	ofstream fout_nc(shd_nc.data());
	if(!fout_nc) { cout << "can't open " << shd_nc << endl; exit(-1); }
	ofstream fout_ot(shd_ot.data());
	if(!fout_ot) { cout << "can't open " << shd_ot << endl; exit(-1); }

	map<string, Read>::iterator iterread = metaseqs.begin();
	for(; iterread != metaseqs.end(); ++ iterread)
	{
		string rid = iterread->first;
		fout_c << '>' << rid << endl;
		fout_nc << '>' << rid << endl;
		fout_ot << '>' << rid << endl;

		map<string, vector<MetaORF> >::iterator itertype = iterread->second.ORFs.begin();
		for(; itertype != iterread->second.ORFs.end(); ++ itertype)
		{
			vector<MetaORF>::iterator iterorf = itertype->second.begin();
			for(; iterorf != itertype->second.end(); ++ iterorf)
			{
				MetaORF mo = (*iterorf);
				if(mo.coding_type == "C")
					fout_c << setw(8) << mo.begin << setw(8) << mo.end << setw(4) << mo.strand 
						<< setw(4) << mo.phase << setw(4) << itertype->first << endl;
				else if(mo.coding_type == "NC")
					fout_nc << setw(8) << mo.begin << setw(8) << mo.end << setw(4) << mo.strand 
						<< setw(4) << mo.phase << setw(4) << itertype->first << endl;
				else if(mo.coding_type == "OT")
					fout_ot << setw(8) << mo.begin << setw(8) << mo.end << setw(4) << mo.strand 
						<< setw(4) << mo.phase << setw(4) << itertype->first << endl;
				else
					cout << "No this type of ORFs." << endl;
			}
		}
	}

	fout_c.close(), fout_nc.close(), fout_ot.close();
}

void MetaSeqs::get_seqs_info(const string& name)
{
	map<int, double> seqs_len_distri;
	map<double, double> seqs_gc_distri;
	map<string, Read>::iterator iterread = metaseqs.begin();
	int min_frag_len = 1000000, max_frag_len = -1;
	double average_gc = 0;
	int total_size = 0;

	for(; iterread != metaseqs.end(); ++ iterread)
	{
		int readsize = (int)iterread->second.seq.size();
		++ seqs_len_distri[ (readsize / 50 + 1)* 50];
		double gc = SeqOp::gccontent(iterread->second.seq);
		++ seqs_gc_distri[ (double(int(gc*20) + 1)) / 20.0 ];

		if(min_frag_len > readsize)
			min_frag_len = readsize;
		if(max_frag_len < readsize)
			max_frag_len = readsize;

		average_gc += gc * (double) readsize;
		total_size += readsize;
	}

	string info_file = name + ".info";
	ofstream fout(info_file.data());
	if(!fout) { cout << "can't open " << info_file << endl; exit(-1); }

	fout<< "minimal fragment: " << min_frag_len << " bps." << endl
		<< "maximal fragment: " << max_frag_len << " bps." << endl
		<< "average gc content: " << average_gc / (double) total_size * 100 << endl;
	fout << "fragments length distribution" << endl;
	double seqssize = (double)MetaSeqs::seqs_num();
	map<int, double>::iterator iterlen = seqs_len_distri.begin();
	for(; iterlen != seqs_len_distri.end(); ++ iterlen)
		fout << setw(10) << iterlen->first << setw(15) << iterlen->second / seqssize * 100 << endl;
	fout << endl << "fragments gc content distribution" << endl;
	map<double, double>::iterator itergc = seqs_gc_distri.begin();
	for(; itergc != seqs_gc_distri.end(); ++ itergc)
		fout << setw(10) << itergc->first << setw(15) << itergc->second / seqssize * 100 << endl;
	fout.close();
}
