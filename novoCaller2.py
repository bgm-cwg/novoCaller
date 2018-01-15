import pysam
import numpy as np
import sys
import os


def GT_ordering_alternate(ALT_count):
	combos=(ALT_count+1)*(ALT_count+2)/2
	ordering=np.empty([combos,2])
	count=0
	for a1 in range(0,ALT_count+1):
		for a2 in range(a1,ALT_count+1):
			ordering[count,0]=a1
			ordering[count,1]=a2
			count=count+1
	return ordering


def row_gen(GT1,GT2,alt_count,mut_rate):
	N=alt_count
	combos=(N+1)*(N+2)/2
	row=np.zeros(combos)
	count=0
	for a1 in range(N+1):
		for a2 in range(N+1):
			for a3 in range(N+1):
				for a4 in range(N+1):
					P=1.0
					if a1==GT1[0]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					if a2==GT1[1]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					if a3==GT2[0]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					if a4==GT2[1]:
						P=P*(1-mut_rate)
					else:
						P=P*mut_rate/N
					count+=1
					
					for b1 in [a1,a2]:
						for b2 in [a3,a4]:
							gt_work=np.sort([b1,b2])
							index=(2*N+3-gt_work[0])*gt_work[0]/2+gt_work[1]-gt_work[0]
							
							row[index]=row[index]+0.25*P
	
	return row



def table_gen(alt_count,mut_rate):
	N=alt_count
	
	
	II_prev=-1
	
	combos=(N+1)*(N+2)/2
	table=np.zeros([combos**2,combos])
	for a1 in range(N+1):
		for a2 in range(a1,N+1):
			for a3 in range(N+1):
				for a4 in range(a3,N+1):
					GT1=[a1,a2]
					GT2=[a3,a4]
					I1=(2*N+3-GT1[0])*GT1[0]/2+GT1[1]-GT1[0]
					I2=(2*N+3-GT2[0])*GT2[0]/2+GT2[1]-GT2[0]
					II=I1*combos+I2
					
					if II<=II_prev:
						print "error in II calc!!!"
						print II_prev,II
					
					
					row=row_gen(GT1,GT2,alt_count,mut_rate)
					table[II,:]=row
	
	
	return table


def GT_likelihood_wrt_allele_calc(ALT_count):
	ordering=GT_ordering_alternate(ALT_count)
	
	combos=(ALT_count+1)*(ALT_count+2)/2
	table=np.zeros([combos,ALT_count+1])*1.
	
	for i in range(combos):
		a1=ordering[i,0]
		a2=ordering[i,1]
		table[i,a1]=table[i,a1]+0.5
		table[i,a2]=table[i,a2]+0.5
	return table


def chr_calc(str_val):
    if str_val[0:3]=="chr":
        str_val=str_val[3:len(str_val)]
    chrom=-1
    rec={"M":0,"X":23,"Y":24}
    try:
        chrom=rec[str_val]
    except:
        try:
            chrom=int(str_val)
        except:
            chrom=25
    return chrom





def get_sam_file_list(list_of_filenames):
	buff=open(list_of_filenames)
	line=buff.readline()
	samfile_list=[]
	while len(line)>1:
		split_line=line.split()
		filename=split_line[0]
		samfile = pysam.AlignmentFile(filename, "rb" )
		samfile_list.append(samfile)
		
		line=buff.readline()
	
	return samfile_list




def get_ADs(samfile,chrom,position_actual,REF,MQ_thresh,BQ_thresh):
	position=position_actual-1
	ADf=np.array([0.,0])
	ADr=np.array([0.,0])
	if chrom==0:
		CC="M"
	if chrom>0 and chrom<=22:
		CC=str(chrom)
	if chrom==23:
		CC="X"
	if chrom==24:
		CC="Y"
	SP=samfile.pileup("chr"+CC, position, position+1)
	try:
		if chrom==0:
			CC="M"
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		if chrom==23:
			CC="X"
		if chrom==24:
			CC="Y"
		SP=samfile.pileup("chr"+CC, position, position+1)
	except ValueError:
		if chrom==0:
			CC="MT"
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		if chrom==23:
			CC="X"
		if chrom==24:
			CC="Y"
		SP=samfile.pileup(CC, position, position+1)
	
	for pileupcolumn in SP:
		if(pileupcolumn.pos==position):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip :
					MQ=pileupread.alignment.mapping_quality
					BQ=pileupread.alignment.query_qualities[pileupread.query_position]
					if MQ>=MQ_thresh and BQ>=BQ_thresh:
						if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper():
							if pileupread.alignment.is_reverse:
								ADr[0]=ADr[0]+1
							else:
								ADf[0]=ADf[0]+1
						else:
							if pileupread.alignment.is_reverse:
								ADr[1]=ADr[1]+1
							else:
								ADf[1]=ADf[1]+1
	return ADf,ADr



def get_ADs_deletion(samfile,chrom,position_actual,REF_0,MQ_thresh,BQ_thresh):
	REF=REF_0
	position=position_actual-1
	ADf=np.array([0.,0])
	ADr=np.array([0.,0])
	try:
		if chrom==0:
			CC="M"
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		if chrom==23:
			CC="X"
		if chrom==24:
			CC="Y"
		SP=samfile.pileup("chr"+CC, position, position+1)
	except ValueError:
		if chrom==0:
			CC="MT"
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		if chrom==23:
			CC="X"
		if chrom==24:
			CC="Y"
		SP=samfile.pileup(CC, position, position+1)
	for pileupcolumn in SP:
		if(pileupcolumn.pos==position):
			for pileupread in pileupcolumn.pileups:
				
				if not pileupread.is_del and not pileupread.is_refskip:
					MQ=pileupread.alignment.mapping_quality
					BQ=pileupread.alignment.query_qualities[pileupread.query_position]
					if MQ>=MQ_thresh and BQ>=BQ_thresh:
						indel_val=pileupread.indel
						if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val>=0:
							if pileupread.alignment.is_reverse:
								ADr[0]=ADr[0]+1
							else:
								ADf[0]=ADf[0]+1
						else:
							if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val<0:
								if pileupread.alignment.is_reverse:
									ADr[1]=ADr[1]+1
								else:
									ADf[1]=ADf[1]+1
	return ADf,ADr




def get_ADs_insertion(samfile,chrom,position_actual,REF_0,MQ_thresh,BQ_thresh):
	REF=REF_0
	position=position_actual-1
	ADf=np.array([0.,0])
	ADr=np.array([0.,0])
	try:
		if chrom==0:
			CC="M"
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		if chrom==23:
			CC="X"
		if chrom==24:
			CC="Y"
		SP=samfile.pileup("chr"+CC, position, position+1)
	except ValueError:
		if chrom==0:
			CC="MT"
		if chrom>0 and chrom<=22:
			CC=str(chrom)
		if chrom==23:
			CC="X"
		if chrom==24:
			CC="Y"
		SP=samfile.pileup(CC, position, position+1)
	for pileupcolumn in SP:
		if(pileupcolumn.pos==position):
			for pileupread in pileupcolumn.pileups:
				
				if not pileupread.is_del and not pileupread.is_refskip:
					MQ=pileupread.alignment.mapping_quality
					BQ=pileupread.alignment.query_qualities[pileupread.query_position]
					if MQ>=MQ_thresh and BQ>=BQ_thresh:
						indel_val=pileupread.indel
						if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val<=0:
							if pileupread.alignment.is_reverse:
								ADr[0]=ADr[0]+1
							else:
								ADf[0]=ADf[0]+1
						else:
							if pileupread.alignment.query_sequence[pileupread.query_position].upper()==REF.upper() and indel_val>0:
								if pileupread.alignment.is_reverse:
									ADr[1]=ADr[1]+1
								else:
									ADf[1]=ADf[1]+1
	return ADf,ADr




def get_ADs_combined(samfile,chrom,position_actual,REF,ALT,MQ_thresh,BQ_thresh):
	split_ALT=ALT.split(",")
	if len(split_ALT)>1:
		ADf,ADr = get_ADs(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
		return ADf,ADr
	if len(REF)>1:
		ADf,ADr = get_ADs_deletion(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
		return ADf,ADr
	if len(ALT)>1:
		ADf,ADr = get_ADs_insertion(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
		return ADf,ADr
	ADf,ADr = get_ADs(samfile,chrom,position_actual,REF[0],MQ_thresh,BQ_thresh)
	return ADf,ADr



def get_all_ADs_combined(unrelated_samfiles,chrom,position_actual,REF,ALT,MQ_thresh,BQ_thresh):
	ADfs=[]
	ADrs=[]
	for samfile in unrelated_samfiles:
		ADf,ADr = get_ADs_combined(samfile,chrom,position_actual,REF,ALT,MQ_thresh,BQ_thresh)
		ADfs.append(ADf)
		ADrs.append(ADr)
		
	ADfs=np.array(ADfs)
	ADrs=np.array(ADrs)
	
	return ADfs,ADrs





def M1_L_calc_aux(rho,k):
	ALT_count=1
	M1_L_k=np.zeros(ALT_count+1)
	default=np.log((1.-rho)/ALT_count)
	M1_L_k=M1_L_k+default
	M1_L_k[k]=np.log(rho)
	return M1_L_k




def M1_L_calc(AD,rho):
	ALT_count=1
	if AD.size-1 != ALT_count:
		print "ERROR in M1_L_calc"
		sys.exit()
	M1_L=[]
	for k in range(ALT_count+1):
		M1_L_k=M1_L_calc_aux(rho,k)
		M1_L.append(M1_L_k)
	return M1_L



def M2_L_calc_aux(M1_L_k,GT_likelihood_wrt_allele_L):
	ALT_count=1
	if M1_L_k.size-1 != ALT_count:
		print "ERROR in M2_L_calc_aux"
		sys.exit()
	combos=(ALT_count+1)*(ALT_count+2)/2
	temp_table=GT_likelihood_wrt_allele_L+np.tile(M1_L_k.reshape([1,ALT_count+1]),[combos,1])
	M2_L_k=np.zeros(combos)
	for i in range(combos):
		row=temp_table[i,:]
		row_max=np.max(row)
		row=row-row_max
		M2_L_k[i]=np.log(np.sum(np.exp(row)))+row_max
	return M2_L_k



def M2_L_calc(M1_L,GT_likelihood_wrt_allele_L):
	ALT_count=1
	if (M1_L[0]).size-1 != ALT_count:
		print "ERROR in M2_L_calc"
		sys.exit()
	M2_L=[]
	for k in range(ALT_count+1):
		M1_L_k=M1_L[k]
		M2_L_k=M2_L_calc_aux(M1_L_k,GT_likelihood_wrt_allele_L)
		M2_L.append(M2_L_k)
	return M2_L






def GT_marg_L_calc(M2_L_f,M2_L_r,ADf,ADr,prior_L):
	GT_marg_L=prior_L
	ALT_count=1
	if ADf.size-1 != ALT_count or ADr.size-1 != ALT_count:
		print "ERROR in GT_marg_L_calc"
		sys.exit()
	for k in range(ALT_count+1):
		M2_L_k=M2_L_f[k]
		GT_marg_L=GT_marg_L+ADf[k]*M2_L_k
	for k in range(ALT_count+1):
		M2_L_k=M2_L_r[k]
		GT_marg_L=GT_marg_L+ADr[k]*M2_L_k
	return GT_marg_L



def M3_L_calc_aux(GT_marg_L,M2_L_k):
	M3_L_k=GT_marg_L-M2_L_k
	return M3_L_k


def M3_L_calc(GT_marg_L,M2_L):
	ALT_count=1
	if len(M2_L)-1 != ALT_count:
		print "ERROR in M3_L_calc"
		sys.exit()
	M3_L=[]
	for k in range(ALT_count+1):
		M2_L_k=M2_L[k]
		M3_L_k=M3_L_calc_aux(GT_marg_L,M2_L_k)
		M3_L.append(M3_L_k)
	return M3_L



def M4_L_calc_aux(M3_L_k,GT_likelihood_wrt_allele_L):
	ALT_count=1
	if (GT_likelihood_wrt_allele_L.shape)[1]-1 != 1:
		print "ERROR in M4_L_calc_aux"
		sys.exit()
	combos=(ALT_count+1)*(ALT_count+2)/2
	temp_table=GT_likelihood_wrt_allele_L+np.tile(M3_L_k.reshape([combos,1]),[1,ALT_count+1])
	M4_L_k=np.zeros(ALT_count+1)
	for i in range(ALT_count+1):
		column=temp_table[:,i]
		column_max=np.max(column)
		column=column-column_max
		M4_L_k[i]=np.log(np.sum(np.exp(column)))+column_max
	return M4_L_k



def M4_L_calc(M3_L,GT_likelihood_wrt_allele_L):
	ALT_count=1
	if (GT_likelihood_wrt_allele_L.shape)[1]-1 != ALT_count:
		print "ERROR in M4_L_calc"
		sys.exit()
	M4_L=[]
	for k in range(ALT_count+1):
		M3_L_k=M3_L[k]
		M4_L_k=M4_L_calc_aux(M3_L_k,GT_likelihood_wrt_allele_L)
		M4_L.append(M4_L_k)
	return M4_L




def A_marg_L_calc(M1_L,M4_L):
	ALT_count=1
	if len(M1_L)-1 != ALT_count:
		print "ERROR in A_marg_L_calc"
		sys.exit()
	A_marg_L=[]
	for k in range(ALT_count+1):
		M1_L_k=M1_L[k]
		M4_L_k=M4_L[k]
		A_marg_L_k=M1_L_k+M4_L_k
		A_marg_L.append(A_marg_L_k)
	return A_marg_L



def T_term_calc_for_rho(A_marg_L,AD):
	if len(A_marg_L)!=AD.size:
		print "ERROR in T_term_calc"
		sys.exit()
	
	ALT_count=AD.size-1
	T1_term=0.
	T2_term=0.
	for k in range(ALT_count+1):
		A_marg_L_k=A_marg_L[k]
		A_marg_temp=np.exp(A_marg_L_k-np.max(A_marg_L_k))
		A_marg=A_marg_temp/np.sum(A_marg_temp)
		
		T1_term=T1_term+A_marg[k]*AD[k]
		T2_term=T2_term+(1.-A_marg[k])*AD[k]
	
	return T1_term,T2_term

def GT_marg_L_to_GT_marg(GT_marg_L):
	M=np.max(GT_marg_L)
	GT_marg_L=GT_marg_L-M
	GT_marg=np.exp(GT_marg_L)
	S=np.sum(GT_marg)
	GT_marg=GT_marg/S
	joint_probty_term=np.log(S)+M
	return GT_marg,joint_probty_term



def EM_step(ADf_list,ADr_list,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D_original,allele_freq):
	D=np.zeros(3)
	D[0]=D_original[0]
	D[1]=D_original[1]
	D[2]=D_original[2]
	
	
	if allele_freq<=0.:
		AF=0.
	else:
		AF=allele_freq
	
	f0=(1.-AF)**2.
	f2=AF**2.
	f1=1.-f0-f2
	D=np.array([f0,f1,f2])*1000.+2.
	
	T1_f=a-1.
	T2_f=b-1.
	T1_r=a-1.
	T2_r=b-1.
	T_for_prior=D-1.
	joint_probty=             (a-1.)*np.log(rho_f_old)+(b-1.)*np.log(1.-rho_f_old)
	joint_probty=joint_probty+(a-1.)*np.log(rho_r_old)+(b-1.)*np.log(1.-rho_r_old)
	for i in range(3):
		joint_probty=joint_probty+(D[i]-1)*prior_L_old[i]
	
	if len(ADf_list)!=len(ADr_list):
		print "ERROR1 in EM_step"
		sys.exit()
	
	for i in range(len(ADf_list)):
		ADf=ADf_list[i]
		ADr=ADr_list[i]
		M1_L_f = M1_L_calc(ADf,rho_f_old)
		M1_L_r = M1_L_calc(ADr,rho_r_old)
		M2_L_f = M2_L_calc(M1_L_f,GT_likelihood_wrt_allele_L)
		M2_L_r = M2_L_calc(M1_L_r,GT_likelihood_wrt_allele_L)
		GT_marg_L = GT_marg_L_calc(M2_L_f,M2_L_r,ADf,ADr,prior_L_old)
		M3_L_f = M3_L_calc(GT_marg_L,M2_L_f)
		M3_L_r = M3_L_calc(GT_marg_L,M2_L_r)
		M4_L_f = M4_L_calc(M3_L_f,GT_likelihood_wrt_allele_L)
		M4_L_r = M4_L_calc(M3_L_r,GT_likelihood_wrt_allele_L)
		A_marg_L_f = A_marg_L_calc(M1_L_f,M4_L_f)
		A_marg_L_r = A_marg_L_calc(M1_L_r,M4_L_r)
		
		T1_term_f,T2_term_f = T_term_calc_for_rho(A_marg_L_f,ADf)
		T1_term_r,T2_term_r = T_term_calc_for_rho(A_marg_L_r,ADr)
		
		T1_f = T1_f + T1_term_f
		T2_f = T2_f + T2_term_f
		T1_r = T1_r + T1_term_r
		T2_r = T2_r + T2_term_r
		
		GT_marg,joint_probty_term = GT_marg_L_to_GT_marg(GT_marg_L)
		joint_probty = joint_probty + joint_probty_term
		
		T_for_prior = T_for_prior + GT_marg
	
	
	rho_f_new=1./(1.+T2_f/T1_f)
	rho_r_new=1./(1.+T2_r/T1_r)
	prior_new=T_for_prior/np.sum(T_for_prior)
	prior_L_new=np.log(prior_new)
	
	return rho_f_new,rho_r_new,prior_L_new,joint_probty




def EM_full(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq):
	
	joint_probty_s=[]
	joint_probty_new=np.nan
	for i in range(3):
		joint_probty_old=joint_probty_new
		rho_f_new,rho_r_new,prior_L_new,joint_probty_new = EM_step(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq)
		rho_f_old=rho_f_new
		rho_r_old=rho_r_new
		prior_L_old=prior_L_new
		joint_probty_s.append(joint_probty_new)
	while np.abs(joint_probty_old-joint_probty_new)>10**-7:
		joint_probty_old=joint_probty_new
		rho_f_new,rho_r_new,prior_L_new,joint_probty_new = EM_step(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq)
		rho_f_old=rho_f_new
		rho_r_old=rho_r_new
		prior_L_old=prior_L_new
		joint_probty_s.append(joint_probty_new)
	return rho_f_new,rho_r_new,prior_L_new,joint_probty_s
	
	

def GTL_L_calc(ADf,ADr,rho_f,rho_r,GT_likelihood_wrt_allele_L):
	M1_L_f = M1_L_calc(ADf,rho_f)
	M1_L_r = M1_L_calc(ADr,rho_r)
	M2_L_f = M2_L_calc(M1_L_f,GT_likelihood_wrt_allele_L)
	M2_L_r = M2_L_calc(M1_L_r,GT_likelihood_wrt_allele_L)
	prior_L = np.zeros(3)
	GTL_L = GT_marg_L_calc(M2_L_f,M2_L_r,ADf,ADr,prior_L)
	GTL_L=GTL_L-np.max(GTL_L)
	return GTL_L
	

def posterior_probty_calc_exact(prior_L,table_L,C_GL_L,M_GL_L,D_GL_L):
	combos=3
	work_column=np.empty(combos**2)
	for I1 in range(combos):
		for I2 in range(combos):
			II=I1*combos+I2
			work_column[II]=prior_L[I1]+prior_L[I2]+M_GL_L[I1]+D_GL_L[I2]
	
	work_table=table_L+np.tile(C_GL_L,[combos**2,1])+np.tile(np.reshape(work_column,[combos**2,1]),[1,combos])
	work_table=work_table-np.max(work_table)
	work_table=np.exp(work_table)
	work_table=work_table/np.sum(work_table)
	PP=np.max(np.array([work_table[0][1],work_table[0][2]]))
	return PP,work_table




def denovo_P_calc(ADfs,ADrs,rho_f,rho_r,GT_likelihood_wrt_allele_L,table_L,prior_L):
	M_GL_L=GTL_L_calc(ADfs[0],ADrs[0],rho_f,rho_r,GT_likelihood_wrt_allele_L)
	D_GL_L=GTL_L_calc(ADfs[1],ADrs[1],rho_f,rho_r,GT_likelihood_wrt_allele_L)
	C_GL_L=GTL_L_calc(ADfs[2],ADrs[2],rho_f,rho_r,GT_likelihood_wrt_allele_L)
	PP,work_table = posterior_probty_calc_exact(prior_L,table_L,C_GL_L,M_GL_L,D_GL_L)
	return PP,work_table






def PP_calc(trio_samfiles,unrelated_samfiles,chrom,pos,REF,ALT,allele_freq,MQ_thresh,BQ_thresh):
	ADfs,ADrs = get_all_ADs_combined(unrelated_samfiles,chrom,pos,REF,ALT,MQ_thresh,BQ_thresh)
	ADfs_U=ADfs
	ADrs_U=ADrs
	rho_f_old=0.8
	rho_r_old=0.8
	prior_old=np.array([1./3,1./3,1./3])
	prior_old=prior_old/np.sum(prior_old)
	prior_L_old=np.log(prior_old)
	GT_likelihood_wrt_allele = GT_likelihood_wrt_allele_calc(1)
	GT_likelihood_wrt_allele_L=np.log(GT_likelihood_wrt_allele)
	a=2.
	b=2.
	D=np.array([2.,2,2])
	
	rho_f_new,rho_r_new,prior_L_new,joint_probty_s = EM_full(ADfs,ADrs,rho_f_old,rho_r_old,prior_L_old,GT_likelihood_wrt_allele_L,a,b,D,allele_freq)
	
	
	AF_unrel=0.
	for i in range(ADfs.shape[0]):
		temp1=GTL_L_calc(ADfs[i],ADrs[i],rho_f_new,rho_r_new,GT_likelihood_wrt_allele_L)
		temp=temp1+prior_L_new
		temp=temp-np.max(temp)
		temp=np.exp(temp)
		temp=temp/np.sum(temp)
		
		AF_unrel=AF_unrel+temp[1]+temp[2]*2.
		
	AF_unrel=AF_unrel/2./ADfs.shape[0]
	
	ADfs,ADrs = get_all_ADs_combined(trio_samfiles,chrom,pos,REF,ALT,MQ_thresh,BQ_thresh)
	
	table=table_gen(1,1e-8)
	table_L=np.log(table)
	PP,work_table = denovo_P_calc(ADfs,ADrs,rho_f_new,rho_r_new,GT_likelihood_wrt_allele_L,table_L,prior_L_new)
	
	return PP,ADfs,ADrs,ADfs_U,ADrs_U,rho_f_new,rho_r_new,prior_L_new,AF_unrel





def cmp_entry(E1,E2):
	if E1[0]>E2[0]:
		return -1
	elif E1[0]==E2[0]:
		return 0
	else:
		return 1




def ALT_read_checker_in_parents(ADfs,ADrs):
	if(len(ADfs)!=3 or len(ADrs)!=3):
		print "ERROR in ALT_read_checker_in_parents"
		sys.exit()
	summ=ADfs[0][1]+ADfs[1][1]+ADrs[0][1]+ADrs[1][1]
	
	if summ>3:
		return False
	else:
		return True





def runner(outfilename,initial_filename,unrelated_filename,trio_filename):
	outbuff_sorted_simple=open(outfilename,'w')
	buff=open(initial_filename,'r')
	line=buff.readline()
	record=[]
	MQ_thresh=-100.
	BQ_thresh=-100.
	count=0
	
	unrelated_samfiles=get_sam_file_list(unrelated_filename)
	trio_samfiles=get_sam_file_list(trio_filename)
	
	
	while line:
		count=count+1
		
		
		print count,
		sys.stdout.flush()
		split_line=line.split()
		chrom=chr_calc(split_line[0])
		pos=int(split_line[1])
		REF=split_line[3]
		ALT=split_line[4]
		INFO=split_line[7]
		split_INFO=INFO.split(";")
		for SS in split_INFO:
			temp="ExAC_AF_computed="
			if SS[0:len(temp)]==temp:
				allele_freq=float(SS[len(temp):])
			temp="MDC="
			if SS[0:len(temp)]==temp:
				MDC=SS[len(temp):]
			temp="CSQ_gene="
			if SS[0:len(temp)]==temp:
				CSQ_gene=SS[len(temp):]
		PP,ADfs,ADrs,ADfs_U,ADrs_U,rho_f_new,rho_r_new,prior_L_new,AF_unrel = PP_calc(trio_samfiles,unrelated_samfiles,chrom,pos,REF,ALT,allele_freq,MQ_thresh,BQ_thresh)
		if AF_unrel<0.01 and ALT_read_checker_in_parents(ADfs,ADrs):
			rec_single=[PP,line,MDC,ADfs,ADrs,ADfs_U,ADrs_U,allele_freq,rho_f_new,rho_r_new,prior_L_new,AF_unrel,CSQ_gene]
			record.append(rec_single)
		
		
		
		line=buff.readline()
		
	
	print
	
	record.sort(cmp_entry)
	count=1
	for rec in record:
		PP=rec[0]
		line=rec[1]
		MDC=rec[2]
		ADfs=rec[3]
		ADrs=rec[4]
		ADfs_U=rec[5]
		ADrs_U=rec[6]
		allele_freq=rec[7]
		rho_f_new=rec[8]
		rho_r_new=rec[9]
		prior_L_new=rec[10]
		AF_unrel=rec[11]
		CSQ_gene=rec[12]
		
		split_line=line.split()
		outbuff_sorted_simple.write(str(count)+")\t"+split_line[0]+"\t"+split_line[1]+"\t"+split_line[3]+"\t"+split_line[4]+"\tAF="+str(allele_freq)+"\t")
		outbuff_sorted_simple.write("rhos= "+str(rho_f_new)+","+str(rho_r_new)+"\t")
		outbuff_sorted_simple.write(("prior=%r" %np.exp(prior_L_new))+"\tPP="+str(PP)+"\t")
		outbuff_sorted_simple.write(("AF_unrel=%r" %AF_unrel)+"\t")
		outbuff_sorted_simple.write(("CSQ_gene=%r" %CSQ_gene)+"\n")
		outbuff_sorted_simple.write("trio:\n")
		for i in range(len(ADfs)):
			outbuff_sorted_simple.write("%r\t%r\n" %(ADfs[i],ADrs[i]))
		outbuff_sorted_simple.write("unrelated:\n")
		for i in range(len(ADfs_U)):
			outbuff_sorted_simple.write("%r\t%r\n" %(ADfs_U[i],ADrs_U[i]))
		
		
		count=count+1
	
	
	

	




if __name__=="__main__":
	
	argv=sys.argv
	for i in range(1,len(argv)):
		if argv[i]=="-O":
			outfilename=argv[i+1]
		if argv[i]=="-I":
			initial_filename=argv[i+1]
		if argv[i]=="-U":
			unrelated_filename=argv[i+1]
		if argv[i]=="-T":
			trio_filename=argv[i+1]
	
	print "outfilename=",outfilename
	print "initial_filename=",initial_filename
	print "unrelated_filename=",unrelated_filename
	print "trio_filename=",trio_filename
	
	runner(outfilename,initial_filename,unrelated_filename,trio_filename)
	









