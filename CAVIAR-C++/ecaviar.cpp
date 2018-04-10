#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <cctype>
#include <unistd.h> 

#include "Util.h"
#include "PostCal.h"
#include "CaviarModel.h"

using namespace std;

bool checkInputFiles(string zFile1, string zFile2, string ldFile1, string ldFile2, string outputFileName) {
	int size = 0;
	int snpCount1 = 0;
	int snpCount2 = 0;
	if(zFile1 == ""){
		cout << "Error: The first Z-socre file is not given" << endl;
		return(false);
	} else if (zFile2 == "") {
		cout << "Error: The second Z-socre file is not given" << endl;
                return(false);
	} else if (ldFile1 == "") {
                cout << "Error: The first LD file is not given" << endl;
                return(false);
        } else if (ldFile2 == "") {
                cout << "Error: The second LD file is not given" << endl;
                return(false);
        } else if (outputFileName=="") {
		cout << "Error: The output file is not give" << endl;
		return(false);
	} 
	fileSize(ldFile1, size);
        snpCount1 = (int)sqrt(size);
        fileSize(ldFile2, size);
        snpCount2 = (int)sqrt(size);
	if(snpCount1 != snpCount2) {
                cout << "Error: The LD files for GWAS and eQTL do not have the same number of SNPs" << endl;
		return(false);
        }
	return(true);	
}

int main( int argc, char *argv[]  ){
	int oc = 0;
	int snpCount = 0;	

	double * stat1;
	double * stat2;
	string * snpNames;

	double NCP = 5.7;
	double rho = 0;
        int totalCausalSNP = 2;
	bool histFlag;
	string tmpLDFileNames = "";
	string ldFile1 = "";
	string zFile1  = "";
	string ldFile2 = "";
        string zFile2  = "";
	string outputFileName = "";

	while ((oc = getopt(argc, argv, "vhl:o:z:r:c:f:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 2.2:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-z ZFILE, --z_file=ZFILE	the z-score and rsID files" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability " << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
				return(0);
			case 'l':
				if(ldFile1=="")
					ldFile1 = string(optarg);
				else
					ldFile2 = string(optarg);
				break;
			case 'o':
				outputFileName = string(optarg);
				break;
			case 'z':
				if(zFile1=="")
                                        zFile1 = string(optarg);
                                else
                                        zFile2 = string(optarg);
				break;
			case 'r':
				rho = atof(optarg);
				break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'f':
                                histFlag = true;
                                break;
			case ':':
			case '?':
			default:
				cout << "Strange" << endl;
				break;
		}
	}

	//program is running
	cout << "@-------------------------------------------------------------@" << endl;
	cout << "| eCAVIAR!     	| 	   v2.2         |  10/Apr/2018        | " << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| (C) 2018 Farhad Hormozdiari, GNU General Public License, v2 |" << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| For documentation, citation & bug-report instructions:      |" << endl;
	cout << "| 		http://genetics.cs.ucla.edu/caviar/           |" << endl;
	cout << "@-------------------------------------------------------------@" << endl;	

	//Check the input is right?
	if( !checkInputFiles(zFile1, zFile2, ldFile1, ldFile2, outputFileName) )
		return(0);
		
	CaviarModel gwasModel(ldFile1, zFile1, outputFileName + "_1", totalCausalSNP, NCP, rho, histFlag);
	gwasModel.run();
	gwasModel.finishUp();
	CaviarModel eqtlModel(ldFile2, zFile2, outputFileName + "_2", totalCausalSNP, NCP, rho, histFlag);	
	eqtlModel.run();
	eqtlModel.finishUp();	

	snpCount  = gwasModel.snpCount;
	snpNames  = new string [snpCount];
	stat1     = new double [snpCount];	
	stat2     = new double [snpCount];	
	importDataFirstColumn(outputFileName+"_1"+"_post", snpNames, 1);
	importDataNthColumn(outputFileName+"_1"+"_post", stat1, 3, 1);	
	importDataNthColumn(outputFileName+"_2"+"_post", stat2, 3, 1);	
	ofstream outfile( (outputFileName+"_col").c_str(), ios::out);	
	double sumCLPP = 0;
	for(int i = 0; i < snpCount; i++) {
		sumCLPP = sumCLPP + stat1[i] * stat2[i];
	}
	outfile << "SNP_ID\tProb_in_pCausalSet\tCLPP" << endl;
	for(int i = 0; i < snpCount; i++) {
		outfile << snpNames[i] << "\t" << (stat1[i] * stat2[i])/sumCLPP << "\t" << stat1[i] * stat2[i] << endl;
	}	
	outfile.close();
	//system(("rm "+ outputFileName + "_1_post").c_str());
	//system(("rm "+ outputFileName + "_2_post").c_str());	
	return 0;
}
