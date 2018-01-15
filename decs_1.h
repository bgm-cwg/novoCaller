#ifndef GAURD_DECS_1
#define GAURD_DECS_1

// marked

std::vector<double> range(const size_t N){
	std::vector<double> A(N,0);
	for(size_t i(0);i<N;++i)
		A[i]=i;
	return A;
}


size_t ALT_count_to_combos(const size_t ALT_count){
	return (ALT_count+1)*(ALT_count+2)/2;
}


size_t combos_to_ALT_count(const size_t C){
	return (pow((C-1)*8+9,0.5)-3)/2;
}



std::vector<std::string> split(const std::string& line){
    std::vector<std::string> split_string;
    std::string::const_iterator i(line.begin()),j(line.begin());
    while(i!=line.end() && j!=line.end()){
        while(isspace(*i) && i!=line.end())
            ++i;
        j=i;
        while(!isspace(*j) && j!=line.end())
            ++j;
        std::string word(i,j);
        if(word.length()>0)
            split_string.push_back(word);
        i=j;
    }
    return split_string;
}


std::string shorten_word(const std::string& word){
	std::string::const_iterator i(word.begin()),j;
	while(isspace(*i) && i!=word.end())
		++i;
	j=i;
	while(!isspace(*j) && j!=word.end())
		++j;
	return std::string(i,j);
}


std::vector<std::string> split(const std::string& line,const char sep){
	std::vector<std::string> split_string;
	std::string::const_iterator i(line.begin()),j(line.begin());
	while(i!=line.end() && j!=line.end()){
		while(*j!=sep && j!=line.end())
			++j;
		std::string full_word(std::string(i,j));
		//std::cout<<full_word.size()<<'\n';
		std::string word(shorten_word(full_word));
		split_string.push_back(word);
		if(j!=line.end())
			i=j+1;
		else
			i=j;
		j=i;
	}
	if(line[line.size()-1]==sep)
		split_string.push_back("");
	return split_string;
}






void col_extractor(const std::string vcf_header,size_t& CHROM_col,size_t& POS_col,size_t& ID_col,size_t& REF_col,size_t& ALT_col,size_t& QUAL_col,size_t& FILTER_col,size_t& INFO_col,size_t& FORMAT_col){
	std::vector<std::string> split_line(split(vcf_header));
	CHROM_col=POS_col=ID_col=REF_col=ALT_col=QUAL_col=FILTER_col=INFO_col=FORMAT_col=1000;
	for(size_t i(0);i<split_line.size();++i){
		std::string F(split_line[i]);
		if(F=="#CHROM")
			CHROM_col=i;
		if(F=="POS")
			POS_col=i;
		if(F=="ID")
			ID_col=i;
		if(F=="REF")
			REF_col=i;
		if(F=="ALT")
			ALT_col=i;
		if(F=="QUAL")
			QUAL_col=i;
		if(F=="FILTER")
			FILTER_col=i;
		if(F=="INFO")
			INFO_col=i;
		if(F=="FORMAT")
			FORMAT_col=i;
	}
}





template<class T>
std::vector<T> LOG_1D(const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=std::log(v[i]);
	return ans;
}


template<class T>
std::vector<std::vector<T> > LOG_2D(const std::vector<std::vector<T> >& v){
	size_t d1(v.size()),d2(v[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=log(v[i][j]);
		}
	}
	return ans;
}


template<class T>
std::vector<T> EXP_1D(const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=std::exp(v[i]);
	return ans;
}

template<class T>
std::vector<std::vector<T> > EXP_2D(const std::vector<std::vector<T> >& v){
	size_t d1(v.size()),d2(v[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=std::exp(v[i][j]);
		}
	}
	return ans;
}

// no need for row extract since = operator copies the entire row vector
template<class T>
std::vector<T> column_extract(const std::vector<std::vector<T> >& P,const size_t C){
	if(C>P[0].size()-1)
		throw(std::domain_error("error in column_extract"));
	const size_t states(P.size());
	std::vector<T> col(states,P[0][0]);
	for(size_t i(0);i<states;++i)
		col[i]=P[i][C];
	return col;
}

template<class T>
void column_insert(std::vector<std::vector<T> >& P,const std::vector<T>& Col,const size_t C){
	if(P.size()!=Col.size())
		throw(std::domain_error("error in column_insert"));
	const size_t states(P.size());
	for(size_t i(0);i<states;++i)
		P[i][C]=Col[i];
}



template<class T>
T sum(const std::vector<T>& v){
	T ans(0);
	for(size_t i(0);i<v.size();++i)
		ans+=v[i];
	return ans;
}

template<class T>
T sum_2D(const std::vector<std::vector<T> >& v){
	T ans(0);
	for(size_t i(0);i<v.size();++i)
		ans+=sum(v[i]);
	return ans;
}


template<class T>
T max(const std::vector<T>& v){
	T ans(v[0]);
	for(size_t i(1);i<v.size();++i)
		if(v[i]>ans)
			ans=v[i];
	return ans;
}

template<class T>
T max_abs(const std::vector<T>& v){
	T ans(std::abs(v[0]));
	for(size_t i(1);i<v.size();++i)
		if(std::abs(v[i])>ans)
			ans=std::abs(v[i]);
	return ans;
}


template<class T>
size_t arg_max(const std::vector<T>& v){
	size_t ans(0);
	for(size_t i(1);i<v.size();++i)
		if(v[i]>v[ans])
			ans=i;
	return ans;
}


template<class T>
size_t arg_max_abs(const std::vector<T>& v){
	size_t ans(0);
	for(size_t i(1);i<v.size();++i)
		if(std::abs(v[i])>std::abs(v[ans]))
			ans=i;
	return ans;
}



template<class T>
T min(const std::vector<T>& v){
	T ans(v[0]);
	for(size_t i(1);i<v.size();++i)
		if(v[i]<ans)
			ans=v[i];
	return ans;
}


template<class T>
std::vector<T> dot(const std::vector<std::vector<T> >& v1,const std::vector<T>& v2){
	if(v1[0].size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v1.size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=sum(v1[i]*v2);
	return ans;
}

template<class T>
std::vector<T> dot(const std::vector<T>& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v2[0].size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=sum(column_extract(v2,i)*v1);
	return ans;
}


template<class T>
std::vector<std::vector<T> > dot(const std::vector<std::vector<T> >& v1,const std::vector<std::vector<T> >& v2){
	if(v1[0].size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	std::vector<std::vector<T> > ans(v1.size(),std::vector<T>(v2[0].size()));
	for(size_t i(0);i<ans.size();++i){
		for(size_t j(0);j<ans[0].size();++j){
			T summ(0);
			for(size_t k(0);k<v1[0].size();++k){
				summ+=v1[i][k]*v2[k][j];
			}
			ans[i][j]=summ;
		}
	}
	return ans;
}



template<class T>
std::vector<T> max_sum_dot(const std::vector<std::vector<T> >& v1,const std::vector<T>& v2){
	if(v1[0].size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v1.size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=max(v1[i]+v2);
	return ans;
}

template<class T>
std::vector<T> max_sum_dot(const std::vector<T>& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("not compatible for dot product"));
	
	const size_t l(v2[0].size());
	std::vector<T> ans(l,0);
	for(size_t i(0);i<l;++i)
		ans[i]=max(column_extract(v2,i)+v1);
	return ans;
}




template<class T>
std::vector<T> operator*(const double S,const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<v.size();++i)
		ans[i]=v[i]*S;
	return ans;
}



template<class T>
std::vector<T> operator*(const std::vector<T>& v,const double S){
	return S*v;
}

template<class T>
std::vector<T> operator+(const double S,const std::vector<T>& v){
	size_t l(v.size());
	std::vector<T> ans(l,v[0]);
	for(size_t i(0);i<v.size();++i)
		ans[i]=v[i]+S;
	return ans;
}

template<class T>
std::vector<T> operator+(const std::vector<T>& v,const double S){
	return S+v;
}


template<class T>
std::vector<T> operator*(const std::vector<T>& v1,const std::vector<T>& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("sizes do not match in multiplication"));
	
	size_t l(v1.size());
	std::vector<T> ans(l,v1[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=v1[i]*v2[i];
	return ans;
}


template<class T>
std::vector<std::vector<T> > operator*(const std::vector<std::vector<T> >& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size() || v1[0].size()!=v2[0].size())
		throw(std::domain_error("sizes do not match in multiplication"));
	
	size_t d1(v1.size()),d2(v1[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v1[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=v1[i][j]*v2[i][j];
		}
	}
	return ans;
}


template<class T>
std::vector<T> operator+(const std::vector<T>& v1,const std::vector<T>& v2){
	if(v1.size()!=v2.size())
		throw(std::domain_error("sizes do not match in addition"));
	
	size_t l(v1.size());
	std::vector<T> ans(l,v1[0]);
	for(size_t i(0);i<l;++i)
		ans[i]=v1[i]+v2[i];
	return ans;
}


template<class T>
std::vector<std::vector<T> > operator+(const std::vector<std::vector<T> >& v1,const std::vector<std::vector<T> >& v2){
	if(v1.size()!=v2.size() || v1[0].size()!=v2[0].size())
		throw(std::domain_error("sizes do not match in addition"));
	
	size_t d1(v1.size()),d2(v1[0].size());
	std::vector<std::vector<T> > ans(d1,std::vector<T>(d2,v1[0][0]));
	for(size_t i(0);i<d1;++i){
		for(size_t j(0);j<d2;++j){
			ans[i][j]=v1[i][j]+v2[i][j];
		}
	}
	return ans;
}


template<class T>
std::ostream& operator<<(std::ostream& out,const std::vector<T>& v){
	if(v.size()>0){
		for(size_t i(0);i<v.size();++i)
			out<<v[i]<<':';
	}
	return out;
}

template<class T>
std::ostream& operator<<(std::ostream& out,const std::list<T>& v){
	typename std::list<T>::const_iterator I(v.begin());
	while(I!=v.end()){
		std::cout<<(*I)<<":";
		++I;
	}
	return out;
}


template<class T>
std::ostream& operator<<(std::ostream& out,const std::vector<std::vector<T> >& v){
	if(v.size()>0){
		for(size_t i(0);i<v.size();++i)
			out<<v[i]<<'\n';
	}
	return out;
}





void col_extractor_for_FORMAT(const std::string& FORMAT,size_t& GT_col,size_t& AD_col,size_t& DP_col,size_t& GQ_col,size_t& PGT_col,size_t& PID_col,size_t& PL_col){
	GT_col=AD_col=DP_col=GQ_col=PGT_col=PID_col=PL_col=1000;
	std::vector<std::string> split_FORMAT(split(FORMAT,':'));
	for(size_t i(0);i<split_FORMAT.size();++i){
		if(split_FORMAT[i]=="GT"){
            GT_col=i;
            continue;
		}
        if(split_FORMAT[i]=="AD"){
            AD_col=i;
            continue;
		}
        if(split_FORMAT[i]=="DP"){
            DP_col=i;
            continue;
		}
        if(split_FORMAT[i]=="GQ"){
            GQ_col=i;
            continue;
		}
        if(split_FORMAT[i]=="PGT"){
            PGT_col=i;
            continue;
		}
        if(split_FORMAT[i]=="PID"){
            PID_col=i;
            continue;
		}
        if(split_FORMAT[i]=="PL"){
            PL_col=i;
            continue;
		}
	}
}



void AD_extractor(std::vector<std::vector<int> >& all_cols,const std::vector<std::string>& split_line,const size_t AD_col,const size_t ALT_count,const size_t FORMAT_length,const size_t FORMAT_col){
	using namespace std;
	for(size_t i(FORMAT_col+1);i<split_line.size();++i){
		const size_t II(i-FORMAT_col-1);
		const string& value(split_line[i]);
		const vector<string> split_value(split(value,':'));
		if(split_value.size()!=FORMAT_length){
			vector<int>& row(all_cols[II]);
			for(size_t j(0);j<2;++j)
				row[j]=-1;
			continue;
		}
		const string& AD_string(split_value[AD_col]);
		const vector<string> split_AD_string(split(AD_string,','));
		if(split_AD_string.size()!=ALT_count+1){
			vector<int>& row(all_cols[II]);
			for(size_t j(0);j<2;++j)
				row[j]=-1;
			continue;
		}
		
		
		vector<int>& row(all_cols[II]);
		row[0]=atoi(split_AD_string[0].c_str());
		int DPA(0);
		for(size_t j(1);j<split_AD_string.size();++j)
			DPA+=atoi(split_AD_string[j].c_str());
		row[1]=DPA;
		
	}
}



std::vector<std::vector<size_t> > GT_ordering_vanilla(const size_t ALT_count,std::vector<std::vector<size_t> >& ordering){
	size_t combos((ALT_count+1)*(ALT_count+2)/2);
	//std::vector<std::vector<size_t> > ordering(combos,std::vector<size_t>(2,0));
	size_t count(0);
	for(size_t a1(0);a1<ALT_count+1;++a1){
		for(size_t a2(a1);a2<ALT_count+1;++a2){
			ordering[count][0]=a1;
			ordering[count][1]=a2;
			count+=1;
		}
	}
	return ordering;
}





std::vector<std::vector<double> > GT_likelihood_wrt_allele_calc(const size_t ALT_count){
	
	size_t combos((ALT_count+1)*(ALT_count+2)/2);
	std::vector<std::vector<size_t> > ordering(combos,std::vector<size_t>(2,0));
	GT_ordering_vanilla(ALT_count,ordering);
	std::vector<std::vector<double> > table(combos,std::vector<double>(ALT_count+1,0.0));
	for(size_t i(0);i<combos;++i){
		size_t a1(ordering[i][0]);
		size_t a2(ordering[i][1]);
		table[i][a1]=table[i][a1]+0.5;
		table[i][a2]=table[i][a2]+0.5;
	}
	return table;
	
}

std::vector<std::vector<std::vector<double> > > GT_likelihood_wrt_allele_L_table_gen(const size_t max_ALT_count){
	std::vector<std::vector<std::vector<double> > > GT_likelihood_wrt_allele_L_table;
	for(size_t i(1);i<=max_ALT_count;++i){
		std::vector<std::vector<double> > GT_likelihood_wrt_allele(GT_likelihood_wrt_allele_calc(i));
		std::vector<std::vector<double> > GT_likelihood_wrt_allele_L(LOG_2D(GT_likelihood_wrt_allele));
		GT_likelihood_wrt_allele_L_table.push_back(GT_likelihood_wrt_allele_L);
	}
	return GT_likelihood_wrt_allele_L_table;
}



class Messages{
	public:
	size_t ALT_count,combos;
	
	std::vector<std::vector<double> > M1_L,M2_L,M3_L,M4_L,A_marg_L;
	std::vector<double> GT_marg_L,GL_L,work_vec_1D_1,work_vec_1D_2;
	// T terms
	double T1_term,T2_term,joint_probty_term;
	std::vector<double> T_term_for_prior;
	
	Messages(const size_t);
	
	void M1_L_calc(const double rho);
	
	void M2_L_calc_aux(const size_t& k,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L);
	
	void M2_L_calc(const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L);
	
	void GT_marg_L_calc(const std::vector<int>& AD,const std::vector<double>& prior_L);
	
	
	void M3_L_calc();
	
	void M4_L_calc_aux(const size_t& k,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L);
	
	void M4_L_calc(const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L);
	
	void A_marg_L_calc();
	
	void updater(const std::vector<int>& AD,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const double& rho,const std::vector<double>& prior_L);
	
	void T_term_calc_for_rho(const std::vector<int>& AD);
	
	void T_term_calc_for_prior();
	
	void GL_L_calc(const std::vector<int>& AD,const double& rho,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L);
	
};

Messages::Messages(const size_t A): ALT_count(A),combos((ALT_count+1)*(ALT_count+2)/2) ,  M1_L(std::vector<std::vector<double> >(ALT_count+1,std::vector<double>(ALT_count+1,0))) , M2_L(std::vector<std::vector<double> >(ALT_count+1,std::vector<double>(combos,0))) , M3_L(std::vector<std::vector<double> >(ALT_count+1,std::vector<double>(combos,0))) , M4_L(std::vector<std::vector<double> >(ALT_count+1,std::vector<double>(ALT_count+1,0))) , GT_marg_L(std::vector<double>(combos,0)) , GL_L(std::vector<double>(combos,0)) , A_marg_L(std::vector<std::vector<double> >(ALT_count+1,std::vector<double>(ALT_count+1,0))) , work_vec_1D_1(std::vector<double>(ALT_count+1,0)) , work_vec_1D_2(std::vector<double>(combos,0)) , T_term_for_prior(std::vector<double>(combos,0)) {}

void Messages::M1_L_calc(const double rho){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t i(0);i<ALT_count+1;++i){
		for(size_t j(0);j<ALT_count+1;++j){
			if(i!=j)
				M1_L[i][j]=log((1-rho)/ALT_count);
			else
				M1_L[i][j]=log(rho);
		}
	}
}



void Messages::M2_L_calc_aux(const size_t& k,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t i(0);i<combos;++i){
		//std::vector<double> row(GT_likelihood_wrt_allele_L[i]+M1_L[k]);
		for(size_t j(0);j<ALT_count+1;++j)
			work_vec_1D_1[j]=GT_likelihood_wrt_allele_L[i][j]+M1_L[k][j];
		double row_max(max(work_vec_1D_1));
		//row=row+(-1)*row_max;
		for(size_t j(0);j<ALT_count+1;++j)
			work_vec_1D_1[j]=work_vec_1D_1[j]-row_max;
		//M2_L[k][i]=log(sum(EXP_1D(work_vec_1D_1)))+row_max;
		double temp_var(0);
		for(size_t j(0);j<ALT_count+1;++j)
			temp_var+=exp(work_vec_1D_1[j]);
		temp_var=log(temp_var)+row_max;
		M2_L[k][i]=temp_var;
	}
}


void Messages::M2_L_calc(const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t k(0);k<ALT_count+1;++k){
		Messages::M2_L_calc_aux(k,GT_likelihood_wrt_allele_L);
	}
	
}




void Messages::GT_marg_L_calc(const std::vector<int>& AD,const std::vector<double>& prior_L){
	for(size_t i(0);i<combos;++i)
		GT_marg_L[i]=prior_L[i];
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t k(0);k<ALT_count+1;++k){
		for(size_t i(0);i<combos;++i){
			GT_marg_L[i]=GT_marg_L[i]+AD[k]*M2_L[k][i];
		}
	}
}






void Messages::M3_L_calc(){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t k(0);k<ALT_count+1;++k){
		//std::vector<double> M3_L_k(GT_marg_L+(-1)*M2_L[k]);
		
		for(size_t i(0);i<combos;++i){
			double M3_L_k_i(GT_marg_L[i]-M2_L[k][i]);
			
			M3_L[k][i]=M3_L_k_i;
			
		}
	}
}


void Messages::M4_L_calc_aux(const size_t& k,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t i(0);i<ALT_count+1;++i){
		//std::vector<double> column(column_extract(GT_likelihood_wrt_allele_L,i)+M3_L[k]);
		for(size_t j(0);j<combos;++j)
			work_vec_1D_2[j]=GT_likelihood_wrt_allele_L[j][i]+M3_L[k][j];
		double column_max(max(work_vec_1D_2));
		//column=column+(-1)*column_max;
		for(size_t j(0);j<combos;++j)
			work_vec_1D_2[j]=work_vec_1D_2[j]-column_max;
		//M4_L[k][i]=log(sum(EXP_1D(work_vec_1D_2)))+column_max;
		double temp_var(0);
		for(size_t j(0);j<combos;++j)
			temp_var+=exp(work_vec_1D_2[j]);
		temp_var=log(temp_var)+column_max;
		M4_L[k][i]=temp_var;
	}
}


void Messages::M4_L_calc(const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t k(0);k<ALT_count+1;++k){
		M4_L_calc_aux(k,GT_likelihood_wrt_allele_L);
	}
}



void Messages::A_marg_L_calc(){
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t k(0);k<ALT_count+1;++k){
		for(size_t i(0);i<ALT_count+1;++i){
			A_marg_L[k][i]=M1_L[k][i]+M4_L[k][i];
		}
	}
}

void Messages::T_term_calc_for_rho(const std::vector<int>& AD){
	//const size_t ALT_count(M1_L[0].size()-1);
	T1_term=0;
	T2_term=0;
	for(size_t k(0);k<ALT_count+1;++k){
		double A_marg_L_k_max(max(A_marg_L[k]));
		double summ(0);
		for(size_t i(0);i<ALT_count+1;++i){
			work_vec_1D_1[i]=exp(A_marg_L[k][i]-A_marg_L_k_max);
			summ+=work_vec_1D_1[i];
		}
		for(size_t i(0);i<ALT_count+1;++i)
			work_vec_1D_1[i]=work_vec_1D_1[i]/summ;
		T1_term+=work_vec_1D_1[k]*AD[k];
		T2_term+=(1-work_vec_1D_1[k])*AD[k];
	}
}


void Messages::T_term_calc_for_prior(){
	// this is exactly GT_marg
	double max_val(max(GT_marg_L));
	double summ(0);
	for(size_t i(0);i<combos;++i){
		work_vec_1D_2[i]=exp(GT_marg_L[i]-max_val);
		summ+=work_vec_1D_2[i];
	}
	joint_probty_term=log(summ)+max_val;
	for(size_t i(0);i<combos;++i)
		T_term_for_prior[i]=work_vec_1D_2[i]/summ;
	
}


void Messages::updater(const std::vector<int>& AD,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const double& rho,const std::vector<double>& prior_L){
	M1_L_calc(rho);
	M2_L_calc(GT_likelihood_wrt_allele_L);
	GT_marg_L_calc(AD,prior_L);
	M3_L_calc();
	M4_L_calc(GT_likelihood_wrt_allele_L);
	A_marg_L_calc();
	T_term_calc_for_rho(AD);
	T_term_calc_for_prior();
}

void Messages::GL_L_calc(const std::vector<int>& AD,const double& rho,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L){
	M1_L_calc(rho);
	M2_L_calc(GT_likelihood_wrt_allele_L);
	
	for(size_t i(0);i<combos;++i)
		GL_L[i]=0;
	//const size_t ALT_count(M1_L[0].size()-1);
	for(size_t k(0);k<ALT_count+1;++k){
		for(size_t i(0);i<combos;++i){
			GL_L[i]=GL_L[i]+AD[k]*M2_L[k][i];
		}
	}
	const double max_val(max(GL_L));
	for(size_t i(0);i<combos;++i)
		GL_L[i]=GL_L[i]-max_val;
}


std::ostream& operator<<(std::ostream& out,const class Messages& M){
	out<<"Messages\n**********\n";
	out<<"ALT_count="<<M.ALT_count<<'\n';
	out<<"combos="<<M.combos<<'\n';
	out<<"M1_L=\n"<<M.M1_L<<'\n';
	out<<"M2_L=\n"<<M.M2_L<<'\n';
	out<<"M3_L=\n"<<M.M3_L<<'\n';
	out<<"M4_L=\n"<<M.M4_L<<'\n';
	out<<"GT_marg_L=\n"<<M.GT_marg_L<<'\n';
	out<<"GL_L=\n"<<M.GL_L<<'\n';
	out<<"A_marg_L=\n"<<M.A_marg_L<<'\n';
	out<<"work_vec_1D_1=\n"<<M.work_vec_1D_1<<'\n';
	out<<"work_vec_1D_2=\n"<<M.work_vec_1D_2<<'\n';
	out<<"T1_term="<<M.T1_term<<'\n';
	out<<"T2_term="<<M.T2_term<<'\n';
	
	return out;
}



class EM_data{
	public:
	double rho_old,rho_new,joint_probty_old,joint_probty_new;
	std::vector<double> prior_old,prior_L_old,prior_new,prior_L_new;
	size_t not_there_count,extra_iterations;
	size_t ALT_count,combos;
	Messages M;
	
	EM_data(const double rho_val,const std::vector<double>& prior_val);
	EM_data(const double rho_val,const size_t C);
	void prior_old_initializer();
	
	
	void EM_step(const std::vector<std::vector<int> >& AD_list,const std::vector<size_t>& parent_cols,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const double a,const double b,const std::vector<double>& D);
	void EM_step(const std::vector<std::vector<int> >& AD_list,const std::vector<size_t>& parent_cols,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const double a,const double b,const std::vector<double>& D,double AF);
	
	void EM_full(const std::vector<std::vector<int> >& AD_list,const std::vector<size_t>& parent_cols,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const size_t min_iterations);
	void EM_full(const std::vector<std::vector<int> >& AD_list,const std::vector<size_t>& parent_cols,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const size_t min_iterations,double AF);
};


void EM_data::EM_full(const std::vector<std::vector<int> >& AD_list,const std::vector<size_t>& parent_cols,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const size_t min_iterations){
	double min_iters;
        if(min_iterations<=3)
                min_iters=3;
        else
                min_iters=min_iterations;
	const double a(2),b(2);
	std::vector<double> D(combos,2);
	for(size_t i(0);i<min_iters;++i)
		EM_step(AD_list,parent_cols,GT_likelihood_wrt_allele_L,a,b,D);
	extra_iterations=0;
	while(std::abs(joint_probty_new-joint_probty_old)>pow(10.0,-4)){
		EM_step(AD_list,parent_cols,GT_likelihood_wrt_allele_L,a,b,D);
		extra_iterations+=1;
	}
}



void EM_data::prior_old_initializer(){
	prior_old=std::vector<double>(combos,0);
	prior_new=std::vector<double>(combos,0);
	const double lambda_val(0.5);
	for(int i(0);i<combos;++i){
		prior_old[i]=exp(-i*lambda_val);
		prior_new[i]=exp(-i*lambda_val);
	}
	prior_old=prior_old*(1/sum(prior_old));
	prior_new=prior_new*(1/sum(prior_new));
	
	prior_L_old=LOG_1D(prior_old);
	prior_L_new=LOG_1D(prior_new);
}


EM_data::EM_data(const double rho_val,const std::vector<double>& prior_val): rho_old(rho_val),rho_new(rho_val),prior_old(prior_val),prior_new(prior_val),combos(prior_val.size()) ,ALT_count(combos_to_ALT_count(combos)),M(ALT_count) {
	prior_old=prior_old*(1/sum(prior_old));
	prior_L_old=LOG_1D(prior_old);
	prior_new=prior_new*(1/sum(prior_new));
	prior_L_new=LOG_1D(prior_new);
}

EM_data::EM_data(const double rho_val,const size_t C): rho_old(rho_val),rho_new(rho_val),combos(C), ALT_count(combos_to_ALT_count(C)),M(ALT_count) {
	prior_old_initializer();
}





void EM_data::EM_step(const std::vector<std::vector<int> >& AD_list,const std::vector<size_t>& parent_cols,const std::vector<std::vector<double> >& GT_likelihood_wrt_allele_L,const double a,const double b,const std::vector<double>& D){
	//const size_t ALT_count(AD_list[0].size()-1);
	//const size_t combos((ALT_count+1)*(ALT_count+2)/2);
	const double rho_lower_limit(1.0/(ALT_count+1));
	if(D.size()!=combos || D.size()!=prior_old.size())
		throw(std::domain_error("error in EM_data.EM_step"));
	
	rho_old=rho_new;
	for(size_t i(0);i<combos;++i){
		prior_old[i]=prior_new[i];
		prior_L_old[i]=prior_L_new[i];
	}
	joint_probty_old=joint_probty_new;
	
	joint_probty_new=(a-1.)*log(rho_old)+(b-1.)*log(1.-rho_old);
	for(size_t i(0);i<combos;++i)
		joint_probty_new+=D[i]*prior_L_old[i];
	
	double T1_for_rho(a-1),T2_for_rho(b-1);
	std::vector<double> T_for_prior(D+(-1));
	
	
	
	
	not_there_count=0;
	for(size_t i(0);i<parent_cols.size();++i){
		if(AD_list[parent_cols[i]][0]!=-1){
			M.updater(AD_list[parent_cols[i]],GT_likelihood_wrt_allele_L,rho_old,prior_L_old);
			T1_for_rho+=M.T1_term;
			T2_for_rho+=M.T2_term;
			joint_probty_new+=M.joint_probty_term;
			for(size_t j(0);j<combos;++j){
				T_for_prior[j]=T_for_prior[j]+M.T_term_for_prior[j];
			}
		}
		else{
			not_there_count+=1;
		}
	}
	
	rho_new=1/(1+T2_for_rho/T1_for_rho);
	if(rho_new<rho_lower_limit)
		rho_new=rho_lower_limit;
	
	double temp(sum(T_for_prior));
	for(size_t i(0);i<combos;++i){
		prior_new[i]=T_for_prior[i]/temp;
		prior_L_new[i]=log(prior_new[i]);
	}
}



std::ostream& operator<<(std::ostream& out,const EM_data& ED){
	out<<"EM_data:\n*********\n";
	out<<"ALT_count="<<ED.ALT_count<<'\n';
	out<<"combos="<<ED.combos<<'\n';
	out<<"rho_old="<<ED.rho_old<<'\n';
	out<<"rho_new="<<ED.rho_new<<'\n';
	out<<"joint_probty_old="<<ED.joint_probty_old<<'\n';
	out<<"joint_probty_new="<<ED.joint_probty_new<<'\n';
	out<<"prior_old=\n"<<ED.prior_old<<'\n';
	out<<"prior_L_old=\n"<<ED.prior_L_old<<'\n';
	out<<"prior_new=\n"<<ED.prior_new<<'\n';
	out<<"prior_L_new=\n"<<ED.prior_L_new<<'\n';
	out<<"not_there_count="<<ED.not_there_count<<'\n';
	out<<"extra_iterations="<<ED.extra_iterations<<'\n';
	out<<"sum(prior_old)="<<sum(ED.prior_old)<<'\n';
	out<<"sum(prior_new)="<<sum(ED.prior_new)<<'\n';
	out<<"ExAC_AF_calculated_old="<<0.5*ED.prior_old[1]+ED.prior_old[2]<<"\n";
	out<<"ExAC_AF_calculated_new="<<0.5*ED.prior_new[1]+ED.prior_new[2]<<"\n";
	
	
	return out;
	
}


void table_gen(std::vector<std::vector<double> >& table,const size_t& ALT_count,const double& mut_rate);


class tables_memory{
	public:
	std::vector<std::vector<std::vector<double> > > tables;
	std::vector<std::vector<std::vector<double> > > tables_L;
	size_t max_ALT_count;
	
	
	tables_memory(const size_t& MAC): max_ALT_count(MAC){
		for(size_t i(1);i<=max_ALT_count;++i){
			const size_t combos(ALT_count_to_combos(i));
			std::vector<std::vector<double> > table(combos*combos,std::vector<double>(combos,0));
			std::vector<std::vector<double> > table_L(combos*combos,std::vector<double>(combos,0));
			tables.push_back(table);
			tables_L.push_back(table_L);
		}
	}
	
	void set_table(const size_t& ALT_count,const double& mut_rate){
		if(ALT_count>max_ALT_count)
			throw(std::domain_error("ERROR in tables_memory.set_table"));
		table_gen(tables[ALT_count-1],ALT_count,mut_rate);
		const size_t combos(ALT_count_to_combos(ALT_count));
		for(size_t i(0);i<combos*combos;++i){
			for(size_t j(0);j<combos;++j){
				tables_L[ALT_count-1][i][j]=std::log(tables[ALT_count-1][i][j]);
			}
		}
	}
	
};



std::ostream& operator<<(std::ostream& out,const tables_memory& obj){
	for(size_t i(0);i<obj.tables.size();++i){
		out<<"*********\n"<<obj.tables[i];
		out<<"*********\n"<<obj.tables_L[i];
	}
	return out;
}



std::vector<std::vector<size_t> > trio_set_calc(const std::string& filename){
	std::ifstream fin(filename.c_str());
	std::string line;
	getline(fin,line);
	std::vector<std::vector<size_t> > trio_set;
	while(line.size()>0){
		const std::vector<std::string> split_line(split(line,':'));
		size_t l;
		if((*(line.end()-1))==':')
			l=split_line.size()-1;
		else
			l=split_line.size();
		std::vector<size_t> row(l,0);
		for(size_t i(0);i<l;++i)
			row[i]=atoi(split_line[i].c_str());
		trio_set.push_back(row);
		
		getline(fin,line);
	}
	return trio_set;
}


std::vector<size_t> parent_cols_calc(const std::string& filename){
	std::ifstream fin(filename.c_str());
	std::string line;
	getline(fin,line);
	std::vector<std::string> split_line(split(line,':'));
	size_t l;
	if((*(line.end()-1))==':')
		l=split_line.size()-1;
	else
		l=split_line.size();
	
	std::vector<size_t> parent_cols(l,0);
	for(size_t i(0);i<l;++i){
		const size_t val(atoi(split_line[i].c_str()));
		parent_cols[i]=val;
	}
	return parent_cols;
}



class vcf_line_cols{
	public:
	size_t CHROM_col,POS_col,ID_col,REF_col,ALT_col,QUAL_col,FILTER_col,INFO_col,FORMAT_col;
	std::vector<std::string> split_line;
	size_t total_candidates,end_col;
	std::vector<size_t> parent_cols;
	std::vector<std::vector<size_t> > trio_set;
	std::vector<std::vector<std::vector<double> > > GT_likelihood_wrt_allele_L_table;
	tables_memory tables_memory_obj;
	
	int CSQ_ExAC_AF_col;
	
	
	
	vcf_line_cols(const std::string& vcf_header,const std::string& trio_id_filename,const std::string& CSQ_line): tables_memory_obj(7){
		split_line=split(vcf_header);
		CHROM_col=POS_col=ID_col=REF_col=ALT_col=QUAL_col=FILTER_col=INFO_col=FORMAT_col=1000;
		for(size_t i(0);i<split_line.size();++i){
			std::string F(split_line[i]);
			if(F=="#CHROM")
				CHROM_col=i;
			if(F=="POS")
				POS_col=i;
			if(F=="ID")
				ID_col=i;
			if(F=="REF")
				REF_col=i;
			if(F=="ALT")
				ALT_col=i;
			if(F=="QUAL")
				QUAL_col=i;
			if(F=="FILTER")
				FILTER_col=i;
			if(F=="INFO")
				INFO_col=i;
			if(F=="FORMAT")
				FORMAT_col=i;
		}
		total_candidates=split_line.size()-(FORMAT_col+1);
		end_col=total_candidates-1;
		
		BGM_set_cols(trio_id_filename);
		
		GT_likelihood_wrt_allele_L_table=GT_likelihood_wrt_allele_L_table_gen(7);
		
		CSQ_ExAC_AF_col=CSQ_ExAC_AF_col_extractor(CSQ_line);
	}
	
	void BGM_set_cols(const std::string& trio_id_filename);
	int CSQ_ExAC_AF_col_extractor(const std::string& CSQ_line);
	static bool CSQ_line_checker(const std::string& line);
	
	
};


bool vcf_line_cols::CSQ_line_checker(const std::string& line){
	using namespace std;
	if(line.size()>=14){
		const string temp(line.begin()+11,line.begin()+14);
		if(temp=="CSQ")
			return 1;
		else
			return 0;
	}
	else
		return 0;
}


int vcf_line_cols::CSQ_ExAC_AF_col_extractor(const std::string& CSQ_line){
	using namespace std;
	string::const_iterator I1,I2;
	I1=CSQ_line.begin();
	while(*I1!='"')
		++I1;
	I2=I1;
	++I2;
	while(*I2!='"')
		++I2;
	string temp(I1+1,I2);
	vector<string> split_temp( split(temp) );
	const size_t L( split_temp.size() );
	temp=split_temp[L-1];
	split_temp=split(temp,'|');
	
	for(size_t I(0);I<split_temp.size();++I){
		if(split_temp[I]=="ExAC_AF")
			return I;
	}
	return -1;
}

void vcf_line_cols::BGM_set_cols(const std::string& trio_id_filename){
	using namespace std;
	
	ifstream fin(trio_id_filename.c_str());
	std::string trio_line;
	getline(fin,trio_line);
	vector<string> trio_split_line(split(trio_line));
	if(trio_split_line.size()!=3)
		throw(domain_error("ERROR in vcf_line_cols::BGM_set_cols"));
	
	vector<size_t> trio(3);
	for(size_t i(FORMAT_col+1);i<split_line.size();++i){
		const string ID( split_line[i] );
		bool mark(0);
		for(size_t j(0);j<trio_split_line.size();++j){
			const string trio_ID(trio_split_line[j]);
			if(trio_ID==ID){
				mark=1;
				break;
			}
		}
		
		if(!mark)
			parent_cols.push_back(i-FORMAT_col-1);
		
		if(trio_split_line[0]==ID)
			trio[0]=i-FORMAT_col-1;
		if(trio_split_line[1]==ID)
			trio[1]=i-FORMAT_col-1;
		if(trio_split_line[2]==ID)
			trio[2]=i-FORMAT_col-1;
		
		
	}
	
	trio_set.push_back(trio);
}


std::ostream& operator<<(std::ostream& out,const vcf_line_cols& obj){
	out<<"vcf_line_cols:\n**************\n";
	out<<obj.CHROM_col<<' '<<obj.POS_col<<' '<<obj.ID_col<<' '<<obj.REF_col<<' '<<obj.ALT_col<<' '<<obj.QUAL_col<<' '<<obj.FILTER_col<<' '<<obj.INFO_col<<' '<<obj.FORMAT_col<<'\n';
	out<<obj.split_line[obj.CHROM_col]<<' '<<obj.split_line[obj.POS_col]<<' '<<obj.split_line[obj.ID_col]<<' '<<obj.split_line[obj.REF_col]<<' '<<obj.split_line[obj.ALT_col]<<' '<<obj.split_line[obj.QUAL_col]<<' '<<obj.split_line[obj.FILTER_col]<<' '<<obj.split_line[obj.INFO_col]<<' '<<obj.split_line[obj.FORMAT_col]<<'\n';
	out<<"total_candidates="<<obj.total_candidates<<'\n';
	out<<"end_col="<<obj.end_col<<'\n';
	out<<"number of parents = "<<obj.parent_cols.size()<<"\n";
	out<<"number of children = "<<obj.total_candidates-obj.parent_cols.size()<<"\n";
	out<<"parent_cols=\n"<<obj.parent_cols<<'\n';
	out<<"trio_set=\n"<<obj.trio_set;
	out<<"CSQ_ExAC_AF_col="<<obj.CSQ_ExAC_AF_col<<"\n";
	
	return out;
}





class max_sum_data{
	public:
	size_t ALT_count,combos;
	std::vector<double> C_mess,M_mess,D_mess;
	std::vector<double> C_mess_inwards,M_mess_inwards,D_mess_inwards;
	std::vector<double> C_marg,M_marg,D_marg;
	std::vector<std::vector<double> > GT_posterior_probty;
	double max_probty;
	size_t C_index,M_index,D_index;
	std::vector<std::vector<size_t> > ordering;
	bool MVC_1,MVC_2,MVC_GoNL;
	size_t max_pos_1,max_pos_2;
	
	max_sum_data(const size_t AC): ALT_count(AC),combos(ALT_count_to_combos(ALT_count)), C_mess(std::vector<double>(combos,0)) , M_mess(std::vector<double>(combos,0)) , D_mess(std::vector<double>(combos,0)) , C_mess_inwards(std::vector<double>(combos,0)) , M_mess_inwards(std::vector<double>(combos,0)) , D_mess_inwards(std::vector<double>(combos,0)) , C_marg(std::vector<double>(combos,0)) , M_marg(std::vector<double>(combos,0)) , D_marg(std::vector<double>(combos,0)) , ordering(std::vector<std::vector<size_t> >(combos,std::vector<size_t>(2,0))) , GT_posterior_probty(std::vector<std::vector<double> >(combos*combos,std::vector<double>(combos,0))) , max_probty(0) , max_pos_1(0) , max_pos_2(0) {GT_ordering_vanilla(ALT_count,ordering);}
	//max_sum_data(){}
	
	void update(const std::vector<double>& C_GL_L,const std::vector<double>& M_GL_L,const std::vector<double>& D_GL_L,const std::vector<double>& prior_L,const std::vector<std::vector<double> >& table_L);
	
	// these below functions are used by the update function
	void mendelian_viol_check_1();
	void mendelian_viol_check_2();
	void mendelian_viol_check_GoNL();
	void GT_posterior_probty_update(const std::vector<double>& C_GL_L,const std::vector<double>& M_GL_L,const std::vector<double>& D_GL_L,const std::vector<double>& prior_L,const std::vector<std::vector<double> >& table_L);
	
};

std::ostream& operator<<(std::ostream& out,const max_sum_data& O);

const char *vinit[]={"intergenic_variant","feature_truncation","feature_elongation","regulatory_region_variant","regulatory_region_amplification","regulatory_region_ablation","TF_binding_site_variant","TFBS_amplification","TFBS_ablation","downstream_gene_variant","upstream_gene_variant","non_coding_transcript_variant","NMD_transcript_variant","intron_variant","non_coding_transcript_exon_variant","3_prime_UTR_variant","5_prime_UTR_variant","mature_miRNA_variant","coding_sequence_variant","synonymous_variant","stop_retained_variant","incomplete_terminal_codon_variant","splice_region_variant","missense_variant","inframe_deletion","inframe_insertion","initiator_codon_variant","stop_lost","frameshift_variant","stop_gained","splice_donor_variant","splice_acceptor_variant","transcript_ablation"};

class CSQ_data{
	public:
	std::map<std::string,double> vals;
	double ExAC_AF;
	std::string MDQ;
	int MDQ_rank;
	std::string gene;
	
	static const std::vector<std::string> ordering;
	
	CSQ_data() {}
	CSQ_data(const std::string& CSQ_line,const int ExAC_AF_col);
	
	void MDQ_update(const std::string& CON,const std::string& GENE_val);
};

const std::vector<std::string> CSQ_data::ordering(vinit,vinit+33);


void CSQ_data::MDQ_update(const std::string& CON,const std::string& GENE_val){
	using namespace std;
	vector<string> split_CON(split(CON,'&'));
	for(size_t i(0);i<split_CON.size();++i){
		const string& conseq(split_CON[i]);
		int val1(0),val2(0);
		for(int j(0);j<ordering.size();++j){
			if(MDQ==ordering[j])
				val1=j;
			if(conseq==ordering[j])
				val2=j;
		}
		if(val2>val1){
			MDQ=conseq;
			MDQ_rank=val2;
			if(GENE_val.size()>0)
				gene=GENE_val;
		}
	}
}


CSQ_data::CSQ_data(const std::string& CSQ_line,const int ExAC_AF_col): MDQ("intergenic_variant"),MDQ_rank(0) {
	using namespace std;
	const vector<string> split_line(split(CSQ_line,','));
	for(size_t i(0);i<split_line.size();++i){
		const string& each_val(split_line[i]);
		const vector<string> split_each_val(split(each_val,'|'));
		/*
		if(split_each_val[3].size()>0){
			bool mark(0);
			for(size_t II(0);II<gene.size();++II){
				if(split_each_val[3]==gene[II]){
					mark=1;
					break;
				}
			}
			if(!mark)
				gene.push_back(split_each_val[3]);
		}
		*/
		string freq_val;
		if(split_each_val.size()>ExAC_AF_col && ExAC_AF_col>0)
			freq_val=split_each_val[ExAC_AF_col];
		
		const string Consequence_val(split_each_val[1]);
		MDQ_update(Consequence_val,split_each_val[3]);
		if(freq_val.size()>0)
			vals[split_each_val[0]]=atof(freq_val.c_str());
		else
			vals[split_each_val[0]]=-1;
	}
	
	map<string,double>::const_iterator II(vals.begin());
	bool mark(1);
	while(II!=vals.end()){
		if(II->second>=0){
			mark=0;
			break;
		}
		++II;
	}
	if(mark){
		ExAC_AF=-1;
	}
	else{
		ExAC_AF=0;
		II=vals.begin();
		while(II!=vals.end()){
			if(II->second>=0)
				ExAC_AF+=II->second;
			++II;
		}
	}
	if(ExAC_AF>1)
		ExAC_AF=1;
}


std::ostream& operator<<(std::ostream& out,const CSQ_data& O){
	using namespace std;
	out<<"CSQ_data:"<<"\n";
	out<<"*********"<<"\n";
	map<string,double>::const_iterator O_I(O.vals.begin());
	while(O_I!=O.vals.end()){
		out<<O_I->first<<"\t"<<O_I->second<<"\n";
		++O_I;
	}
	out<<"CSQ_data::ExAC_AF="<<O.ExAC_AF<<"\n";
	out<<"CSQ_data::MDQ="<<O.MDQ<<"\n";
	out<<"CSQ_data::MDQ_rank="<<O.MDQ_rank<<"\n";
	out<<"CSQ_data::gene="<<O.gene<<"\n";
	return out;
}


class INFO_data{
	public:
	std::vector<std::string> split_INFO;
	double MQ,FS,DP,AF,ExAC_AF;
	CSQ_data CSQ_data_obj;
	
	INFO_data(){}
	INFO_data(const std::string& INFO,const int ExAC_AF_col);
	
};

INFO_data::INFO_data(const std::string& INFO,const int ExAC_AF_col): split_INFO(split(INFO,';')) ,MQ(-1),FS(-1),DP(-1),AF(-1),ExAC_AF(-1) {
	std::vector<std::string>::const_iterator i(split_INFO.begin());
	for(;i!=split_INFO.end();++i){
		const std::string& temp(*i);
		if(temp[0]=='M' && temp[1]=='Q' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			MQ=atof(new_temp.c_str());
			continue;
		}
		if(temp[0]=='F' && temp[1]=='S' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			FS=atof(new_temp.c_str());
			continue;
		}
		if(temp[0]=='D' && temp[1]=='P' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			DP=atof(new_temp.c_str());
			continue;
		}
		if(temp[0]=='A' && temp[1]=='F' && temp[2]=='='){
			const std::string new_temp(temp.begin()+3,temp.end());
			AF=atof(new_temp.c_str());
			continue;
		}
		if(std::string(temp.begin(),temp.begin()+8)=="ExAC_AF="){
			const std::string new_temp(temp.begin()+8,temp.end());
			ExAC_AF=atof(new_temp.c_str());
		}
		if(std::string(temp.begin(),temp.begin()+4)=="CSQ="){
			const std::string new_temp(temp.begin()+4,temp.end());
			CSQ_data_obj=CSQ_data(new_temp,ExAC_AF_col);
		}
	}
}
std::ostream& operator<<(std::ostream& out,const INFO_data& O){
	out<<"INFO_data:\n**********\n";
	out<<"MQ="<<O.MQ<<'\n';
	out<<"FS="<<O.FS<<'\n';
	out<<"DP="<<O.DP<<'\n';
	out<<"AF="<<O.AF<<'\n';
	out<<"ExAC_AF="<<O.ExAC_AF<<"\n";
	out<<O.CSQ_data_obj;
	return out;
}


std::ostream& operator<<(std::ostream& out,const INFO_data& O);

class location{
	public:
	size_t chrom,pos;
	
	location(const size_t& C,const size_t& P): chrom(C),pos(P) {}
	location(){}
	
};
bool operator<(const location& loc1,const location& loc2);
bool operator==(const location& loc1,const location& loc2);
bool operator>(const location& loc1,const location& loc2);
bool operator<=(const location& loc1,const location& loc2);
bool operator>=(const location& loc1,const location& loc2);
bool operator!=(const location& loc1,const location& loc2);
long operator-(const location& L1,const location& L2);


std::ostream& operator<<(std::ostream& out,const location& loc){
	out<<loc.chrom<<"\t"<<loc.pos;
	return out;
}

size_t str_to_size_t(const std::string& S);

class HG19_query{
	public:
	
	location LOC;
	char BASE;
	size_t line_start;
	
	std::ifstream fin;
	std::string line;
	std::string prev_line;
	
	
	HG19_query(): LOC(0,1) , fin("/net/data/hg19/hg19.fasta") , line_start(1) {
		getline(fin,line);
		prev_line=line;
		getline(fin,line);
		BASE=line[0];
		if(BASE>=97)
			BASE-=32;
	}
	
	HG19_query(const std::string& fasta_filename): LOC(0,1) , fin(fasta_filename.c_str()) , line_start(1) {
		getline(fin,line);
		prev_line=line;
		getline(fin,line);
		BASE=line[0];
		if(BASE>=97)
			BASE-=32;
	}
	
	char base_extractor(const location& loc);
	
	// below here, the functions are not used directly
	void func_1(const location& loc);
	
};


std::ostream& operator<<(std::ostream& out,const HG19_query& O){
	out<<"HG19_query\n*********\n";
	out<<"LOC="<<O.LOC<<"\n";
	out<<"BASE="<<O.BASE<<"\n";
	out<<"line_start="<<O.line_start<<"\n";
	out<<"prev_line="<<O.prev_line<<"\t:"<<O.prev_line.size()<<"\n";
	out<<"line     ="<<O.line<<"\t:"<<O.line.size()<<"\n";
}

class mut_rates_per_line{
	// ATGC
	public:
	std::vector<std::vector<double> > SNPs;
	double CGtoTG,CGtoCA;
	
	mut_rates_per_line(): SNPs(std::vector<std::vector<double> >(4,std::vector<double>(4,-1))) {}
	mut_rates_per_line(const std::vector<std::string>& split_line,const std::vector<std::string>& header_line): SNPs(std::vector<std::vector<double> >(4,std::vector<double>(4,-1))) {
		set_vals(split_line,header_line);
	}
	
	void set_vals(const std::vector<std::string>& split_line,const std::vector<std::string>& header_line);
	
	size_t base_I(const char& B);
	
};


class mut_rate_tracker{
	public:
	std::ifstream fin;
	std::string line;
	std::vector<std::string> header_line;
	std::vector<std::string> split_line;
	location loc_start,loc_end;
	mut_rates_per_line O;
	const double default_mut_rate;
	
	mut_rate_tracker(): fin("/net/home/amohanty/projects/pgcs/denovo/all_cpp/GoNL_denovo_study_data/local_mutation_rate.bias_corrected.SEXAVG.bed") , default_mut_rate(1.2e-8) {
		getline(fin,line);
		header_line=split(line);
		getline(fin,line);
		split_line=split(line);
		const size_t C(str_to_size_t(split_line[0])),P1(atoi(split_line[1].c_str())),P2(atoi(split_line[2].c_str()));
		loc_start=location(C,P1);
		loc_end=location(C,P2);
		O=mut_rates_per_line(split_line,header_line);
	}
	
	// set_vals need not be used directly, but through compute_mut_rate
	void set_vals(const location& loc){
		if(loc<=loc_end)
			return;
		
		//std::cout<<"GOING INTO WHILE LOOP\n";
		while(loc>loc_end){
			if(line.size()==0){
				//std::cout<<"BREAKING_1\n";
				break;
			}
			getline(fin,line);
			if(line.size()==0){
				//std::cout<<"BREAKING_2\n";
				break;
			}
			split_line=split(line);
			const size_t C(str_to_size_t(split_line[0])),P1(atoi(split_line[1].c_str())),P2(atoi(split_line[2].c_str()));
			loc_start=location(C,P1);
			loc_end=location(C,P2);
			O=mut_rates_per_line(split_line,header_line);
		}
		
	}
	
	
	double compute_mut_rate(const location& loc,const char& B_l,const char& B_m,const char& B_r,const char& ALT);
	
	
};



class const_mut_rate_per_line;
class const_mut_rate_handler{
	public:
	std::vector<const_mut_rate_per_line> all;
	
	const_mut_rate_handler();
	double mut_rate_calc(const char B_l,const char B_m,const char B_r,const char ALT) const;
	
	//not used directly, whats below
	size_t index_calc(const char B_l,const char B_m,const char B_r,const char ALT) const;
};


class each_line_data{
	// an object of this class is declred once and mem_set is used after that
	public:
	std::vector<std::string> split_line;
	std::string ALT;
	std::vector<std::string> split_ALT;
	const size_t combos;
	size_t ALT_count,combos_actual;
	size_t GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col;
	std::string FORMAT;
	std::string INFO;
	size_t FORMAT_length;
	std::vector<std::vector<int> > AD_list;
	std::vector<std::vector<double> > GL_L_list;
	max_sum_data max_sum_data_obj;
	INFO_data INFO_data_obj;
	
	std::string REF;
	location LOC;
	
	
	
	
	
	
	
	each_line_data(const std::string& line,const vcf_line_cols& obj);
	
	void mem_set(const std::string& line,const vcf_line_cols& obj);
	
	void GL_L_list_updater(const double& rho,const vcf_line_cols& obj);
	
	
	
	double denovo_BGM(const vcf_line_cols& vcf_line_cols_obj,const std::vector<double>& prior_L);
};


double each_line_data::denovo_BGM(const vcf_line_cols& O,const std::vector<double>& prior_L){
	
	//const double mut_rate( (1.2e-8)*SCALE);
	//O.tables_memory_obj.set_table(1,mut_rate);
	const std::vector<std::vector<double> >& table_L(O.tables_memory_obj.tables_L[0]);
	
	const std::vector<size_t>& trio(O.trio_set[0]);
	const size_t& M_I(trio[0]);
	const size_t& D_I(trio[1]);
	const size_t& C_I(trio[2]);
	const std::vector<int>& M_AD(AD_list[M_I]);
	const std::vector<int>& D_AD(AD_list[D_I]);
	const std::vector<int>& C_AD(AD_list[C_I]);
	const std::vector<double>& M_GL_L(GL_L_list[M_I]);
	const std::vector<double>& D_GL_L(GL_L_list[D_I]);
	const std::vector<double>& C_GL_L(GL_L_list[C_I]);
	double P(-1);
	if( M_AD[0]>-1 && D_AD[0]>-1 && C_AD[0]>-1 && M_GL_L[0]<1 && D_GL_L[0]<1 && C_GL_L[0]<1 ){
		max_sum_data_obj.update(C_GL_L,M_GL_L,D_GL_L,prior_L,table_L);
		P= max_sum_data_obj.GT_posterior_probty[0][1] ;
	}
	return P;
	
}



each_line_data::each_line_data(const std::string& line,const vcf_line_cols& obj): split_line(split(line)) , ALT(split_line[obj.ALT_col]),split_ALT(split(ALT,',')) , ALT_count(split_ALT.size()) , combos(3) , combos_actual(ALT_count_to_combos(ALT_count)) , FORMAT(split_line[obj.FORMAT_col]) , INFO(split_line[obj.INFO_col]) , FORMAT_length(split(FORMAT,':').size()) , AD_list(obj.total_candidates,std::vector<int>(2,100)) , GL_L_list(AD_list.size(),std::vector<double>(combos,100)) , max_sum_data_obj(1) , INFO_data_obj(INFO,obj.CSQ_ExAC_AF_col) , REF(split_line[obj.REF_col]) , LOC(str_to_size_t(split_line[obj.CHROM_col]),atoi(split_line[obj.POS_col].c_str())) {
	//ALT=split_line[obj.ALT_col];
	//split_ALT=split(ALT,',');
	//ALT_count=split_ALT.size();
	//combos=(ALT_count+1)*(ALT_count+2)/2;
	//FORMAT=split_line[obj.FORMAT_col];
	//FORMAT_length=split(FORMAT,':').size();
	col_extractor_for_FORMAT(FORMAT,GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col);
	//AD_list=std::vector<std::vector<int> >(obj.total_candidates,std::vector<int>(ALT_count+1,100));
	AD_extractor(AD_list,split_line,AD_col,ALT_count,FORMAT_length,obj.FORMAT_col);
	//GL_L_list=std::vector<std::vector<double> >(AD_list.size(),std::vector<double>(combos,100));
	
	
}




void each_line_data::mem_set(const std::string& line,const vcf_line_cols& obj){
	split_line=split(line);
	ALT=split_line[obj.ALT_col];
	split_ALT=split(ALT,',');
	ALT_count=split_ALT.size();
	combos_actual=(ALT_count+1)*(ALT_count+2)/2;
	FORMAT=split_line[obj.FORMAT_col];
	INFO=split_line[obj.INFO_col];
	FORMAT_length=split(FORMAT,':').size();
	col_extractor_for_FORMAT(FORMAT,GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col);
	/*
	if(AD_list[0].size()!=ALT_count+1){
		AD_list=std::vector<std::vector<int> >(obj.total_candidates,std::vector<int>(ALT_count+1,100));
		GL_L_list=std::vector<std::vector<double> >(AD_list.size(),std::vector<double>(combos,100));
		//max_sum_data_obj=max_sum_data(ALT_count);
	}
	*/
	AD_extractor(AD_list,split_line,AD_col,ALT_count,FORMAT_length,obj.FORMAT_col);
	INFO_data_obj=INFO_data(INFO,obj.CSQ_ExAC_AF_col);
	REF=split_line[obj.REF_col];
	LOC=location(str_to_size_t(split_line[obj.CHROM_col]),atoi(split_line[obj.POS_col].c_str()));
}



void each_line_data::GL_L_list_updater(const double& rho,const vcf_line_cols& obj){
	Messages M(1);
	for(size_t i(0);i<GL_L_list.size();++i){
		if(AD_list[i][0]>-1){
			M.GL_L_calc(AD_list[i],rho,obj.GT_likelihood_wrt_allele_L_table[0]);
			for(size_t j(0);j<combos;++j){
				GL_L_list[i][j]=M.GL_L[j];
			}
		}
	}
}


std::ostream& operator<<(std::ostream& out,const each_line_data& obj){
	out<<"each_line_data:\n***************\n";
	out<<"LOC="<<obj.LOC<<"\n";
	out<<"REF="<<obj.REF<<"\n";
	out<<"ALT="<<obj.ALT<<'\n';
	out<<"split_ALT="<<obj.split_ALT<<'\n';
	out<<"ALT_count, combos,combos_actual = "<<obj.ALT_count<<", "<<obj.combos<<", "<<obj.combos_actual<<'\n';
	out<<"GT_col,AD_col,DP_col,GQ_col,PGT_col,PID_col,PL_col="<<obj.GT_col<<' '<<obj.AD_col<<' '<<obj.DP_col<<' '<<obj.GQ_col<<' '<<obj.PGT_col<<' '<<obj.PID_col<<' '<<obj.PL_col<<'\n';
	out<<"FORMAT="<<obj.FORMAT<<'\n';
	out<<"FORMAT_length="<<obj.FORMAT_length<<'\n';
	out<<"INFO=  "<<obj.INFO<<'\n';
	out<<"INFO_data_obj=\n"<<obj.INFO_data_obj;
	out<<"AD_list=\n";
	for(size_t i(0);i<obj.AD_list.size();++i)
		out<<i<<")\t"<<obj.AD_list[i]<<"\t"<<obj.GL_L_list[i]<<'\n';
	out<<"AD_list.size()="<<obj.AD_list.size()<<'\n';
	//out<<"max_sum_data_obj=\n"<<obj.max_sum_data_obj;
	
	return out;
}



// note to self: stay the f#@* away from iterators in nested loops. they somehow did not work in this function. Spent 3 f#@*ing hours on this POS function before getting it to work.
// also uses vanilla ordering
void row_gen(std::vector<double>& row,const std::vector<size_t>& GT1,const std::vector<size_t>& GT2,const size_t& ALT_count,const double& mut_rate){
	const size_t N(ALT_count);
	const size_t combos(ALT_count_to_combos(N));
	if(row.size()!=combos)
		throw(std::domain_error("error in row_gen"));
	
	for(size_t i(0);i<combos;++i)
		row[i]=0;
	
	size_t count(0);
	for(size_t a1(0);a1<N+1;++a1){
		for(size_t a2(0);a2<N+1;++a2){
			for(size_t a3(0);a3<N+1;++a3){
				for(size_t a4(0);a4<N+1;++a4){
					double P(1);
					if(a1==GT1[0])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					if(a2==GT1[1])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					if(a3==GT2[0])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					if(a4==GT2[1])
						P=P*(1-mut_rate);
					else
						P=P*mut_rate/N;
					
					std::vector<size_t> B1(2,0);B1[0]=a1;B1[1]=a2;
					std::vector<size_t> B2(2,0);B2[0]=a3;B2[1]=a4;
					++count;
					
					for(size_t b1(0);b1<2;++b1){
						for(size_t b2(0);b2<2;++b2){
							std::vector<size_t> gt_work(2,0);
							gt_work[0]=B1[b1];
							gt_work[1]=B2[b2];
							sort(gt_work.begin(),gt_work.end());
							const size_t index(  N*gt_work[0]+gt_work[1] - int(gt_work[0]*gt_work[0]-gt_work[0])/2  );
							row[index]=row[index]+0.25*P;
						}
					}
				}
			}
		}
	}
	
	const double summ(sum(row));
	for(size_t i(0);i<row.size();++i)
		row[i]=row[i]/summ;
}



void table_gen(std::vector<std::vector<double> >& table,const size_t& ALT_count,const double& mut_rate){
	const size_t N(ALT_count);
	const size_t combos(ALT_count_to_combos(N));
	if(table.size()!=combos*combos || table[0].size()!=combos)
		throw(std::domain_error("error in table_gen"));
		
	for(size_t a1(0);a1<N+1;++a1){
		for(size_t a2(a1);a2<N+1;++a2){
			for(size_t a3(0);a3<N+1;++a3){
				for(size_t a4(a3);a4<N+1;++a4){
					std::vector<size_t> GT1(2,0);GT1[0]=a1;GT1[1]=a2;
					std::vector<size_t> GT2(2,0);GT2[0]=a3;GT2[1]=a4;
					const size_t I1(  N*GT1[0]+GT1[1] - int(GT1[0]*GT1[0]-GT1[0])/2  );
					const size_t I2(  N*GT2[0]+GT2[1] - int(GT2[0]*GT2[0]-GT2[0])/2  );
					const size_t II(I1*combos+I2);
					
					row_gen(table[II],GT1,GT2,ALT_count,mut_rate);
					
				}
			}
		}
	}
}





void child_mess_calc(std::vector<double>& M3,const std::vector<double>& M1,const std::vector<double>& M2,const std::vector<std::vector<double> >& table_L,const size_t& ALT_count){
	const size_t combos(ALT_count_to_combos(ALT_count));
	if(table_L.size()!=combos*combos || table_L[0].size()!=combos)
		throw(std::domain_error("ERROR in child_mess_calc"));
	
	for(size_t i(0);i<combos;++i){
		double maximum(-std::numeric_limits<double>::infinity());
		for(size_t I1(0);I1<combos;++I1){
			for(size_t I2(0);I2<combos;++I2){
				const size_t II(I1*combos+I2);
				if(table_L[II][i]+M1[I1]+M2[I2]>maximum){
					maximum=table_L[II][i]+M1[I1]+M2[I2];
				}
			}
		}
		M3[i]=maximum;
	}
}




void parent_mess_calc(std::vector<double>& M_Mom,const std::vector<double>& M_C,const std::vector<double>& M_Dad,const std::vector<std::vector<double> >& table_L,const size_t& ALT_count){
	const size_t combos(ALT_count_to_combos(ALT_count));
	for(size_t i(0);i<combos;++i){
		double maximum(-std::numeric_limits<double>::infinity());
		for(size_t j(0);j<combos;++j){
			for(size_t k(0);k<combos;++k){
				const double temp( table_L[j+i*combos][k]+M_C[k]+M_Dad[j] );
				if(temp>maximum)
					maximum=temp;
			}
		}
		M_Mom[i]=maximum;
	}
}





void max_sum_data::update(const std::vector<double>& C_GL_L,const std::vector<double>& M_GL_L,const std::vector<double>& D_GL_L,const std::vector<double>& prior_L,const std::vector<std::vector<double> >& table_L){
	for(size_t i(0);i<combos;++i){
		C_mess_inwards[i]=C_GL_L[i];
		M_mess_inwards[i]=M_GL_L[i]+prior_L[i];
		D_mess_inwards[i]=D_GL_L[i]+prior_L[i];
	}
	child_mess_calc(C_mess,M_mess_inwards,D_mess_inwards,table_L,ALT_count);
	parent_mess_calc(M_mess,C_mess_inwards,D_mess_inwards,table_L,ALT_count);
	parent_mess_calc(D_mess,C_mess_inwards,M_mess_inwards,table_L,ALT_count);
	for(size_t i(0);i<combos;++i){
		C_marg[i]=C_mess[i]+C_mess_inwards[i];
		M_marg[i]=M_mess[i]+M_mess_inwards[i];
		D_marg[i]=D_mess[i]+D_mess_inwards[i];
	}
	C_index=arg_max(C_marg);
	M_index=arg_max(M_marg);
	D_index=arg_max(D_marg);
	
	mendelian_viol_check_1();
	mendelian_viol_check_2();
	mendelian_viol_check_GoNL();
	
	if(1){
		GT_posterior_probty_update(C_GL_L,M_GL_L,D_GL_L,prior_L,table_L);
	}
}


std::ostream& operator<<(std::ostream& out,const max_sum_data& O){
	out<<"max_sum_data:\n************\n";
	out<<"ALT_count,combos = "<<O.ALT_count<<", "<<O.combos<<'\n';
	out<<"C_marg="<<O.C_marg<<'\n';
	out<<"M_marg="<<O.M_marg<<'\n';
	out<<"D_marg="<<O.D_marg<<'\n';
	out<<"C_index,M_index,D_index = "<<O.C_index<<", "<<O.M_index<<", "<<O.D_index<<'\n';
	out<<"MVC_GoNL="<<O.MVC_GoNL<<'\n';
	out<<"MVC_2="<<O.MVC_2<<'\n';
	out<<"MVC_1="<<O.MVC_1<<'\n';
	out<<"ordering=\n"<<O.ordering;
	out<<"max_probty="<<O.max_probty<<'\n';
	out<<"max_pos_1,max_pos_2 = "<<O.max_pos_1<<", "<<O.max_pos_2<<'\n';
	out<<"GT_posterior_probty=\n"<<O.GT_posterior_probty;
	out<<"sum_2D(GT_posterior_probty)="<<sum_2D(O.GT_posterior_probty)<<'\n';
}


void max_sum_data::mendelian_viol_check_2(){
	const std::vector<size_t>& GT_c(ordering[C_index]);
	const std::vector<size_t>& GT_m(ordering[M_index]);
	const std::vector<size_t>& GT_f(ordering[D_index]);
	
	for(size_t i(0);i<2;++i){
		for(size_t j(0);j<2;++j){
			const size_t a1(GT_m[i]),a2(GT_f[j]);
			std::vector<size_t> GT_temp(2,0);
			GT_temp[0]=a1;
			GT_temp[1]=a2;
			sort(GT_temp.begin(),GT_temp.end());
			if(GT_c[0]==GT_temp[0] && GT_c[1]==GT_temp[1]){
				//return 0;
				MVC_2=0;
				return;
			}
			else{
				size_t summ_1(0),summ_2(0);
				for(size_t I(0);I<2;++I){
					if(GT_c[I]!=GT_temp[I] && GT_c[I]==0)
						++summ_1;
					if(GT_c[I]!=GT_temp[I])
						++summ_2;
				}
				if(summ_1==summ_2){
					//return 0;
					MVC_2=0;
					return;
				}
			}
		}
	}
	//return 1;
	MVC_2=1;
	return;
}


void max_sum_data::mendelian_viol_check_1(){
	const std::vector<size_t>& GT_c(ordering[C_index]);
	const std::vector<size_t>& GT_m(ordering[M_index]);
	const std::vector<size_t>& GT_f(ordering[D_index]);
	for(size_t i(0);i<2;++i){
		for(size_t j(0);j<2;++j){
			const size_t a1(GT_m[i]),a2(GT_f[j]);
			std::vector<size_t> GT_temp(2,0);
			GT_temp[0]=a1;
			GT_temp[1]=a2;
			sort(GT_temp.begin(),GT_temp.end());
			if(GT_c[0]==GT_temp[0] && GT_c[1]==GT_temp[1]){
				MVC_1=0;
				return;
			}
		}
	}	
	MVC_1=1;
	return;
}

void max_sum_data::mendelian_viol_check_GoNL(){
	const std::vector<size_t>& GT_c(ordering[C_index]);
	const std::vector<size_t>& GT_m(ordering[M_index]);
	const std::vector<size_t>& GT_f(ordering[D_index]);
	if(GT_c[0]==0 && GT_c[1]!=0 && GT_m[0]==0 && GT_m[1]==0 && GT_f[0]==0 && GT_f[1]==0){
		MVC_GoNL=1;
		return;
	}
	else{
		MVC_GoNL=0;
		return;
	}
}


void max_sum_data::GT_posterior_probty_update(const std::vector<double>& C_GL_L,const std::vector<double>& M_GL_L,const std::vector<double>& D_GL_L,const std::vector<double>& prior_L,const std::vector<std::vector<double> >& table_L){
	double maximum(-std::numeric_limits<double>::infinity());
	max_pos_1=0;
	max_pos_2=0;
	for(size_t I1(0);I1<combos;++I1){
		for(size_t I2(0);I2<combos;++I2){
			const size_t II(I1*combos+I2);
			for(size_t k(0);k<combos;++k){
				GT_posterior_probty[II][k]=table_L[II][k]+C_GL_L[k]+prior_L[I1]+prior_L[I2]+M_GL_L[I1]+D_GL_L[I2];
				if(GT_posterior_probty[II][k]>maximum){
					maximum=GT_posterior_probty[II][k];
					max_pos_1=II;
					max_pos_2=k;
				}
			}
		}
	}
	
	double summ(0);
	for(size_t i(0);i<combos*combos;++i){
		for(size_t j(0);j<combos;++j){
			GT_posterior_probty[i][j]=exp(GT_posterior_probty[i][j]-maximum);
			summ+=GT_posterior_probty[i][j];
		}
	}
	max_probty=GT_posterior_probty[max_pos_1][max_pos_2]/summ;
	
	for(size_t i(0);i<combos*combos;++i){
		for(size_t j(0);j<combos;++j){
			GT_posterior_probty[i][j]=GT_posterior_probty[i][j]/summ;
		}
	}
	
	
}






bool operator<(const location& loc1,const location& loc2){
	if(loc1.chrom<loc2.chrom)
		return 1;
	else{
		if(loc1.chrom==loc2.chrom){
			if(loc1.pos<loc2.pos)
				return 1;
			else
				return 0;
		}
		else
			return 0;
	}
}



bool operator==(const location& loc1,const location& loc2){
	if(loc1.chrom==loc2.chrom && loc1.pos==loc2.pos)
		return 1;
	else
		return 0;
}

bool operator!=(const location& loc1,const location& loc2){
	return !(loc1==loc2);
}



bool operator>(const location& loc1,const location& loc2){
	if( !(loc1<loc2) && !(loc1==loc2) )
		return 1;
	else
		return 0;
}

bool operator<=(const location& loc1,const location& loc2){
	if(loc1<loc2 || loc1==loc2)
		return 1;
	else
		return 0;
}

bool operator>=(const location& loc1,const location& loc2){
	if(loc1>loc2 || loc1==loc2)
		return 1;
	else
		return 0;
}





char HG19_query::base_extractor(const location& loc){
	if(loc==LOC)
		return BASE;
	
	if(loc.chrom<LOC.chrom)
		throw(std::domain_error("ERROR_1 in HG19_query::base_extractor"));
	
	
	if(loc.chrom==LOC.chrom && loc.pos<line_start){
		const size_t diff(line_start-loc.pos);
		const size_t l(prev_line.size());
		if(diff>l)
			throw(std::domain_error("ERROR_2 in HG19_query::base_extractor"));
		
		char temp(prev_line[l-diff]);
		
		if(temp>=97)
			temp-=32;
		
		return temp;
	}
	
	
	if(loc.chrom==LOC.chrom)
		func_1(loc);
	else{
		while(line[0]!='>' || str_to_size_t(std::string(line.begin()+1,line.end()))!=loc.chrom )
			getline(fin,line);
		prev_line=line;
		getline(fin,line);
		
		LOC=location(loc.chrom,1);
		BASE=line[0];
		if(BASE>=97)
			BASE-=32;
		line_start=1;
		
		
		func_1(loc);
	}
	
	return BASE;
}



void HG19_query::func_1(const location& loc){
	while(line_start+line.size()-1<loc.pos){
		line_start+=line.size();
		prev_line=line;
		getline(fin,line);
		if(line[0]=='>')
			throw(std::domain_error("ERROR in HG19_query.func_1"));
	}
	BASE=line[loc.pos-line_start];
	if(BASE>=97)
		BASE-=32;
	LOC=loc;
}




size_t str_to_size_t(const std::string& S){
	std::string S2,S_temp;
	if(S.size()>3)
		S_temp=std::string(S.begin(),S.begin()+3);
	else
		S_temp=S;
	
	if(S_temp=="chr")
		S2=std::string(S.begin()+3,S.end());
	else
		S2=S;
	
	//std::cout<<"S2="<<S2<<'\n';
	
	const size_t num(atoi(S2.c_str()));
	if(num>=1 && num<=22)
		return num;
	else
		if(S2==std::string("M") || S2==std::string("m"))
			return 0;
		else
			if(S2==std::string("X") || S2==std::string("x"))
				return 23;
			else
				if(S2==std::string("Y") || S2==std::string("y"))
					return 24;
				else
					return 25;
}



char base_extractor(std::istream& in,const size_t pos){
        // we use 1 indexing instead of zero indexing
        // this function is just a tester, not to be used for the real thing
        std::string line;
        std::getline(in,line);
        std::getline(in,line);
        //std::cout<<line<<":\t"<<line.size()<<"\n";
        size_t count(0);
        if(line.size()+count>pos){
                return line[pos-int(count)-1];
        }
        while(line.size()>0){
                count+=line.size();
                std::getline(in,line);
                //std::cout<<line<<":\t"<<line.size()<<"\n";
                if(line.size()+count>pos){
                        return line[pos-int(count)-1];
                }
        }
        throw(std::domain_error("pos>chrom length"));
}








void mut_rates_per_line::set_vals(const std::vector<std::string>& split_line,const std::vector<std::string>& header_line){
	CGtoTG=atof(split_line[3].c_str());
	CGtoCA=atof(split_line[4].c_str());
	for(size_t i(5);i<header_line.size();++i){
		const size_t l(header_line[i].size());
		const char B1(header_line[i][0]);
		const char B2(header_line[i][l-1]);
		const size_t I1(base_I(B1)),I2(base_I(B2));
		if(I1==I2)
			throw(std::domain_error("ERROR in mut_rates_per_line::set_vals"));
		SNPs[I1][I2]=atof(split_line[i].c_str());
	}
}



std::ostream& operator<<(std::ostream& out,const mut_rates_per_line& O){
	out<<"CGtoTG="<<O.CGtoTG<<"\n";
	out<<"CGtoCA="<<O.CGtoCA<<"\n";
	out<<"ordering is ATGC\n";
	out<<"SNPs=\n";
	out<<O.SNPs;
}



size_t mut_rates_per_line::base_I(const char& B){
		if(B=='A')
			return 0;
		if(B=='T')
			return 1;
		if(B=='G')
			return 2;
		if(B=='C')
			return 3;
}








double mut_rate_tracker::compute_mut_rate(const location& loc,const char& B_l,const char& B_m,const char& B_r,const char& ALT){
	set_vals(loc);
	if(loc<loc_start || loc>loc_end)
		return default_mut_rate;
	
	if(B_m=='C' && B_r=='G' && ALT=='T')
		return O.CGtoTG;
	if(B_m=='G' && B_l=='C' && ALT=='A')
		return O.CGtoCA;
	
	const size_t I1(O.base_I(B_m)),I2(O.base_I(ALT));
	return O.SNPs[I1][I2];
}


std::ostream& operator<<(std::ostream& out,const mut_rate_tracker& O){
	out<<"loc_start="<<O.loc_start<<"\n";
	out<<"loc_end="<<O.loc_end<<"\n";
	out<<"mut_rates_per_line=\n"<<O.O;
	out<<"default_mut_rate="<<O.default_mut_rate<<"\n";
}



long operator-(const location& L1,const location& L2){
	if(L1.chrom==L2.chrom)
		return long(L1.pos)-long(L2.pos);
	else{
		if(L1.chrom>L2.chrom)
			return std::numeric_limits<long>::max();
		else
			return -std::numeric_limits<long>::max();
	}
}




class const_mut_rate_per_line{
	public:
	std::string triplet;
	char ALT;
	double mut_rate;
	const_mut_rate_per_line(const std::string&);
};

const_mut_rate_per_line::const_mut_rate_per_line(const std::string& line){
	const std::vector<std::string> split_line(split(line));
	triplet=split_line[0];
	const std::string temp(split_line[3]);
	ALT=(*(temp.end()-1));
	mut_rate=atof(split_line[14].c_str());
	if(triplet[1]==ALT)
		throw(std::domain_error("ERROR in const_mut_rate_per_line::const_mut_rate_per_line"));
}



std::ostream& operator<<(std::ostream& out,const const_mut_rate_per_line& O){
	out<<O.triplet<<"\t"<<O.ALT<<"\t"<<O.mut_rate;
	return out;
}



double const_mut_rate_handler::mut_rate_calc(const char B_l,const char B_m,const char B_r,const char ALT) const{
	const size_t I( index_calc(B_l,B_m,B_r,ALT) );
	const double M( all[I].mut_rate );
	return M;
}


const_mut_rate_handler::const_mut_rate_handler(){
	std::ifstream fin("/net/home/amohanty/projects/pgcs/denovo/all_cpp/GoNL_denovo_study_data/triplets-2.txt");
	std::string line;
	getline(fin,line);
	
	getline(fin,line);
	while(line.size()>0){
		const const_mut_rate_per_line O(line);
		all.push_back(O);
		
		getline(fin,line);
	}
	
}


size_t const_mut_rate_handler::index_calc(const char B_l,const char B_m,const char B_r,const char ALT) const{
	size_t count(0);
	char temp[]={'A','C','G','T'};
	std::vector<char> T(temp,temp+4);
	for(size_t i(0);i<4;++i){
		if(T[i]!=B_m)
			++count;
		if(T[i]==ALT)
			break;
	}
	--count;
	for(size_t i(0);i<4;++i){
		if(T[i]==B_l)
			count+=i*48;
		if(T[i]==B_m)
			count+=i*12;
		if(T[i]==B_r)
			count+=i*3;
		
	}
	return count;
}





std::ostream& operator<<(std::ostream& out,const const_mut_rate_handler& O){
	for(size_t i(0);i<O.all.size();++i)
		out<<O.all[i]<<"\n";
	return out;
}


double lbeta(const double a,const double b){
	return boost::math::lgamma(a)+boost::math::lgamma(b)-boost::math::lgamma(a+b);
}








#endif















