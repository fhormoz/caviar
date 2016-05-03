#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <cctype>
#include <unistd.h> 

#include "Util.h"

using namespace std;

bool checkInputFiles(string zFile1, string zFile2, string ldFile1, string ldFile2, string outputFileName, string rho) {
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
	} else if(rho == "") {
		cout << "Error: The rho value is not give" << endl;
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
	int tmpSize = 0;
	int snpCount1  = 0;	
	int snpCount2  = 0;

	double * sigma;
	double * stat1;
	double * stat2;
	string * snpNames;

	string rho = "";
        string totalCausalSNP = "2";
	string tmpLDFileNames = "";
	string ldFile1 = "";
	string zFile1  = "";
	string ldFile2 = "";
        string zFile2  = "";
	string outputFileName = "";

	while ((oc = getopt(argc, argv, "vhl:o:z:r:c:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 1.0:" << endl;
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
				rho = string(optarg);
				break;
			case 'c':
				totalCausalSNP = string(optarg);
				break;
			case 'f':
                                //histFlag = true;
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
	cout << "| eCAVIAR!     	| 	   v1.0         |  22/Apr/2016 | " << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| (C) 2016 Farhad Hormozdiari, GNU General Public License, v2 |" << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| For documentation, citation & bug-report instructions:      |" << endl;
	cout << "| 		http://genetics.cs.ucla.edu/caviar/            |" << endl;
	cout << "@-------------------------------------------------------------@" << endl;	

	//Check the input is right?
	if( !checkInputFiles(zFile1, zFile2, ldFile1, ldFile2, outputFileName, rho) )
		return(0);
	fileSize(ldFile1, tmpSize);
        snpCount1 = (int)sqrt(tmpSize);
	fileSize(ldFile2, tmpSize);
	snpCount2 = (int)sqrt(tmpSize);

	string command1 = "./CAVIAR -z " + zFile1 + " -l " + ldFile1 + " -c " + totalCausalSNP  + " -r " + rho  + " -o " + outputFileName+"_1";  	
	string command2 = "./CAVIAR -z " + zFile1 + " -l " + ldFile1 + " -c " + totalCausalSNP  + " -r " + rho  + " -o " + outputFileName+"_2";
	system(command1.c_str());
        system(command2.c_str());

	snpNames  = new string [snpCount1];
	stat1     = new double [snpCount1];	
	stat2     = new double [snpCount1];	

	importDataFirstColumn(outputFileName+"_1"+"_post", snpNames);
	importDataNthColumn(outputFileName+"_1"+"_post", stat1, 3);	
	importDataNthColumn(outputFileName+"_2"+"_post", stat2, 3);	

	ofstream outfile(outputFileName.c_str(), ios::out | ios::app);	
	for(int i = 0; i < snpCount1; i++) {
		outfile << snpNames[i] << "\t" << stat1[i] * stat2[i] << endl;
	}	
	outfile.close();

	return 0;
}
