#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>

#include "definitions.h"
#include "tests.h"



int main(int argc,char** argv){
	std::string arg1,arg2,arg3,arg4,arg5,arg6,arg7;
	for(int i(0);i<argc;++i){
		if(std::string(argv[i])=="-a1")
			arg1=argv[i+1];
		if(std::string(argv[i])=="-a2")
			arg2=argv[i+1];
		if(std::string(argv[i])=="-a3")
			arg3=argv[i+1];
		if(std::string(argv[i])=="-a4")
			arg4=argv[i+1];
		if(std::string(argv[i])=="-a5")
			arg5=argv[i+1];
		if(std::string(argv[i])=="-a6")
			arg6=argv[i+1];
		if(std::string(argv[i])=="-a7")
			arg7=argv[i+1];
	}
	
	
	
	{
		using namespace std;
		string infilename,trio_ID_filename,outfilename,fasta_filename,X_choice;
		double PP_thresh,ExAC_thresh;
		for(int i(0);i<argc;++i){
			if(std::string(argv[i])=="-I")
				infilename=argv[i+1];
			if(std::string(argv[i])=="-T")
				trio_ID_filename=argv[i+1];
			if(std::string(argv[i])=="-O")
				outfilename=argv[i+1];
			if(std::string(argv[i])=="-X")
				X_choice=argv[i+1];
			if(std::string(argv[i])=="-P")
				PP_thresh=atof( argv[i+1] );
			if(std::string(argv[i])=="-E")
				ExAC_thresh=atof( argv[i+1] );
		}
		cout<<"infilename="<<infilename<<"\n";
		cout<<"trio_ID_filename="<<trio_ID_filename<<"\n";
		cout<<"outfilename="<<outfilename<<"\n";
		cout<<"X_choice="<<X_choice<<"\n";
		cout<<"PP_thresh="<<PP_thresh<<"\n";
		cout<<"ExAC_thresh="<<ExAC_thresh<<"\n";
		
		ifstream fin(infilename.c_str());
		ofstream fout(outfilename.c_str());
		denovo_worker(fin,trio_ID_filename,fout,X_choice,PP_thresh,ExAC_thresh);
	}
	
	
	return 0;
}





