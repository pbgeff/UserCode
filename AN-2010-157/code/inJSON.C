//root script to make a JSON filter.
//Usage: root inJSON on afs enable machine.
//You can change the JSON input file name for other runs.
//#define INPUT "./Cert_132440-136297_7TeV_StreamExpress_Collisions10_JSON.txt"
#define INPUT "./Cert_132440-137028_7TeV_StreamExpress_Collisions10_JSON.txt"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
void inJSON(){
	ifstream orgJSON(INPUT);
	ofstream fout("inJSON.h");
	stringstream ss;
	string str;
	while(!orgJSON.eof()){
		getline(orgJSON,str,' ');//orgJSON>>str;			
		fout<<str;	
	}
	orgJSON.close();
	fout.close();	
	char inChar;
	int inInt;
	
	vector<unsigned int> VRunLumi;
	ifstream fin("inJSON.h");
	while(fin >> inChar){
		char next = fin.peek();
		if(	next == '1' || next == '2' || next == '3' ||
				next == '4' || next == '5' || next == '6' ||
				next == '7' || next == '8' || next == '9' || 
				next == '0'){				
			fin >>inInt;
			//cout<<inInt<<endl;
			VRunLumi.push_back(inInt);						
		}
	}
	fin.close();	
	ofstream fout2("inJSON.h");	
	fout2<<"//JSON Filter for nTuples"<<endl;
	fout2<<"#ifndef INJSON_H"<<endl;
	fout2<<"#define INJSON_H"<<endl;
	fout2<<"bool inJSON(unsigned int run,unsigned int lumiblock){\n";
	for(unsigned int i=0;i<VRunLumi.size();){
		//fout2<<VRunLumi[i]<<endl;
		if(VRunLumi[i] > 100000){
			if(i<3) fout2<<"\tif(";
			else fout2<<endl<<"\t|| ";
			
			fout2<<"(run=="<<VRunLumi[i];
			fout2<<" && (";
			for(unsigned int j=i;j<VRunLumi.size()-2
			&& VRunLumi[j+1] < 100000 && VRunLumi[j+2] < 100000;j=j+2){
				fout2<<"(lumiblock>="<<VRunLumi[j+1]<<" && lumiblock<="<<VRunLumi[j+2]<<")";
				//fout2<<endl<<VRunLumi[j+3]<<endl;
				if(VRunLumi[j+3] < 100000 && VRunLumi[j+3] > 0 ){
					fout2 << " || ";
				}
			
			}
			i=j+1;
		}
		fout2<<"))";
	}
	fout2<<")\n  {return true;}\n";
	fout2<<"  else return false;\n}\n";
	fout2<<"#endif //INJSON_H\n";
	fout2.close();
}
