#ifndef GAURD_TESTS
#define GAURD_TESTS


std::ostream& den_writer(const double P,const EM_data& EM_data_O,const vcf_line_cols& O,const each_line_data& ELDO,std::ostream& out){
	out.precision(20);
	const std::vector<std::string>& split_line(ELDO.split_line);
	const size_t INFO_col(O.INFO_col);
	for(size_t i(0);i<split_line.size()-1;++i){
		out<<split_line[i];
		if(i!=INFO_col){
			out<<"\t";
			continue;
		}
		out<<";P_denovo="<<P<<";rho="<<EM_data_O.rho_new<<";GF="<<EM_data_O.prior_new<<";ExAC_AF_computed="<<ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF;
		out<<";MDC="<<ELDO.INFO_data_obj.CSQ_data_obj.MDQ;
		out<<";MDC_rank="<<ELDO.INFO_data_obj.CSQ_data_obj.MDQ_rank;
		out<<";CSQ_gene="<<ELDO.INFO_data_obj.CSQ_data_obj.gene;
		out<<"\t";
	}
	out<<split_line[split_line.size()-1];
	return out;
}

class calls_record{
	public:
	static std::list<std::string> all_records;
	static size_t INFO_col;
	
	static void updater(const double P,const EM_data& EM_data_O,const vcf_line_cols& O,const each_line_data& ELDO);
	static bool compare_individual_records(const std::string&,const std::string&);
	static void sort_record();
	static void writer(std::ostream& out);
};

void calls_record::writer(std::ostream& out){
	using namespace std;
	list<string>::const_iterator I( all_records.begin() );
	for(;I!=all_records.end();++I)
		out<<(*I)<<"\n";
}


size_t calls_record::INFO_col=0;
std::list<std::string> calls_record::all_records(0);


void calls_record::sort_record(){
	using namespace std;
	all_records.sort(compare_individual_records);
}


bool calls_record::compare_individual_records(const std::string& line1,const std::string& line2){
	using namespace std;
	vector<string> split_line,split_INFO;
	string INFO;
	double P1,P2;
	
	split_line=split(line1);
	INFO=split_line[INFO_col];
	split_INFO=split(INFO,';');
	for(size_t I(0);I<split_INFO.size();++I){
		const string element(split_INFO[I]);
		if( element.size()>9 ){
			const string temp( element.begin(),element.begin()+9 );
			if(temp=="P_denovo="){
				const string val( element.begin()+9,element.end() );
				P1=atof( val.c_str() );
			}
		}
	}
	
	split_line=split(line2);
	INFO=split_line[INFO_col];
	split_INFO=split(INFO,';');
	for(size_t I(0);I<split_INFO.size();++I){
		const string element(split_INFO[I]);
		if( element.size()>9 ){
			const string temp( element.begin(),element.begin()+9 );
			if(temp=="P_denovo="){
				const string val( element.begin()+9,element.end() );
				P2=atof( val.c_str() );
			}
		}
	}
	
	return P1>P2;
}


void calls_record::updater(const double P,const EM_data& EM_data_O,const vcf_line_cols& O,const each_line_data& ELDO){
	using namespace std;
	const std::vector<std::string>& split_line(ELDO.split_line);
	INFO_col=O.INFO_col;
	string temp;
	for(size_t i(0);i<split_line.size()-1;++i){
		temp+=split_line[i];
		if(i!=INFO_col){
			temp+="\t";
			continue;
		}
		ostringstream OS_P;
		OS_P.precision(30);
		OS_P<<P;
		temp+= ";P_denovo="+OS_P.str();
		
		ostringstream OS_rho;
		OS_rho.precision(30);
		OS_rho<<EM_data_O.rho_new;
		temp+= ";rho="+OS_rho.str();
		
		ostringstream OS_GF;
		OS_GF.precision(30);
		OS_GF<<EM_data_O.prior_new;
		temp+= ";GF="+OS_GF.str();
		
		ostringstream OS_EAC;
		OS_EAC.precision(30);
		OS_EAC<<ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF;
		temp+= ";ExAC_AF_computed="+OS_EAC.str();
		
		ostringstream OS_MDC;
		OS_MDC.precision(30);
		OS_MDC<<ELDO.INFO_data_obj.CSQ_data_obj.MDQ;
		temp+= ";MDC="+OS_MDC.str();
		
		ostringstream OS_MDC_rank;
		OS_MDC_rank.precision(30);
		OS_MDC_rank<<ELDO.INFO_data_obj.CSQ_data_obj.MDQ_rank;
		temp+= ";MDC_rank="+OS_MDC_rank.str();
		
		ostringstream OS_CSQ_gene;
		OS_CSQ_gene.precision(30);
		OS_CSQ_gene<<ELDO.INFO_data_obj.CSQ_data_obj.gene;
		temp+= ";CSQ_gene="+OS_CSQ_gene.str();
		
		temp+="\t";
	}
	
	temp+=split_line[split_line.size()-1];
	
	all_records.push_back(temp);
}



void denovo_worker(std::istream& fin,const std::string& trio_ID_filename,std::ostream& fout,const std::string& X_choice,const double PP_thresh,const double ExAC_thresh){
	using namespace std;
	string line;
	getline(fin,line);
	while(! vcf_line_cols::CSQ_line_checker(line) )
		getline(fin,line);
	const string CSQ_line(line);
	while(line[0]=='#' && line[1]=='#')
		getline(fin,line);
	const string header_line(line);
	vcf_line_cols O(header_line,trio_ID_filename,CSQ_line);
	O.tables_memory_obj.set_table(1,1e-8);
	cout<<O;
	
	size_t DD(0);
	
	size_t upper_limit;
	if(X_choice=="1")
		upper_limit=23;
	else
		upper_limit=22;
	
	getline(fin,line);
	each_line_data ELDO(line,O);
	EM_data EM_data_obj(0.8,ELDO.combos);
	const double& ExAC_AF(ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF);
	const int& MDQ_rank(ELDO.INFO_data_obj.CSQ_data_obj.MDQ_rank);
	if(ELDO.LOC.chrom>=1 && ELDO.LOC.chrom<=upper_limit && ExAC_AF<ExAC_thresh && MDQ_rank>=17){
		EM_data_obj.EM_full(ELDO.AD_list,O.parent_cols,O.GT_likelihood_wrt_allele_L_table[0],1);
		ELDO.GL_L_list_updater(EM_data_obj.rho_new,O);
		double PP(  ELDO.denovo_BGM(O,EM_data_obj.prior_L_new)  );
		if(PP>PP_thresh){
			calls_record::updater(PP,EM_data_obj,O,ELDO);
			++DD;
		}
	}
	size_t count(0);
	getline(fin,line);
	while(line.size()>0){
		++count;
		if(count%1000000==0){
			cout<<count<<" ";
			cout.flush();
		}
		
		ELDO.mem_set(line,O);
		EM_data EM_data_obj(0.8,ELDO.combos);
		const double& ExAC_AF(ELDO.INFO_data_obj.CSQ_data_obj.ExAC_AF);
		const int& MDQ_rank(ELDO.INFO_data_obj.CSQ_data_obj.MDQ_rank);
		if(ELDO.LOC.chrom>=1 && ELDO.LOC.chrom<=upper_limit && ExAC_AF<ExAC_thresh && MDQ_rank>=17){
			EM_data_obj.EM_full(ELDO.AD_list,O.parent_cols,O.GT_likelihood_wrt_allele_L_table[0],1);
			ELDO.GL_L_list_updater(EM_data_obj.rho_new,O);
			double PP(  ELDO.denovo_BGM(O,EM_data_obj.prior_L_new) );
			if(PP>PP_thresh){
				calls_record::updater(PP,EM_data_obj,O,ELDO);
				++DD;
			}
		}
		
		getline(fin,line);
	}
	cout<<"\n";
	cout<<"count="<<count<<"\n";
	
	calls_record::sort_record();
	
	calls_record::writer(fout);
}



void tester_1(){
	using namespace std;
	
	const string filename("/net/data/bgm/cases/bgm0058_wes/analysis/upstream/bgm0058_final.vep.vcf");
	const string trio_id_filename("trio_ID_bgm0058.txt");
	ifstream fin(filename.c_str());
	string line;
	getline(fin,line);
	while(! vcf_line_cols::CSQ_line_checker(line) )
		getline(fin,line);
	const string CSQ_line(line);
	while(line[0]=='#' && line[1]=='#')
		getline(fin,line);
	const string header_line(line);
	
	vcf_line_cols O(header_line,trio_id_filename,CSQ_line);
	
	getline(fin,line);
	each_line_data ELDO(line,O);
	
	cout<<O;
	cout<<ELDO;
}




#endif







