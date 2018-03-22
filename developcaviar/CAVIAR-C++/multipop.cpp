#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <cctype>
#include <unistd.h> 

#include "Util.h"
#include "PostCalMultiPop.h"
#include "CaviarModelMultiPop.h"

using namespace std;

bool checkInputFiles(vector <string> zFiles, vector <string> ldFiles, string outputFileName) {
	int size = 0;
	if (zFiles.size() == 0 || ldFiles.size() ==0) {
		cout << "You did not specify any LD file of Z file" << endl;
		return (false);
	} if (zFiles.size() != ldFiles.size()) {
                cout << "The number of LD files and Z files do not match" << endl;
                return (false);
        }
	//TODO: CHECK THE PARAMETER ARE CORRECT
	return(true);	
}

int main( int argc, char *argv[]  ){
	int oc = 0;
	int snpCount = 0;	

	double NCP = 5.7;
	double rho = 0.95;
        int totalCausalSNP = 2;
	bool histFlag;
	vector <string> ldFiles;
	vector <string> zFiles;
	string outputFileName = "";

	while ((oc = getopt(argc, argv, "vhl:o:z:r:c:f:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 2.0:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-z ZFILE, --z_file=ZFILE	the z-score and rsID files" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability " << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
				cout << "Example 2 populcations command:" << endl;
				cout << "./mupCAVIAR   -l <LDFILE1> -l <LDFILE2> -z <ZFILE1> -z <ZFILE2> -r 0.95 -f <NUMCAUSAL> -o <OUTPUTFILE>" << endl;
				cout << "Example 3 populcations command:" << endl;
				cout << "./mupCAVIAR   -l <LDFILE1> -l <LDFILE2> -l <LDFILE3> ";
				cout << " -z <ZFILE1> -z <ZFILE2>  -z <ZFILE3> -r 0.95 -f <NUMCAUSAL> -o <OUTPUTFILE>" << endl;
				return(0);
			case 'l':
				ldFiles.push_back(string(optarg));
				break;
			case 'o':
				outputFileName = string(optarg);
				break;
			case 'z':
				zFiles.push_back(string(optarg));
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
	cout << "| MultiplePopAVIAR!     	| 	   v2.0         |  25/Aug/2017 | " << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| (C) 2017 Farhad Hormozdiari, GNU General Public License, v2 |" << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| For documentation, citation & bug-report instructions:      |" << endl;
	cout << "| 		http://genetics.cs.ucla.edu/caviar/            |" << endl;
	cout << "@-------------------------------------------------------------@" << endl;	

	//Check the input is right?
	if( !checkInputFiles(zFiles, ldFiles, outputFileName) )
		return(0);
		
	CaviarModelMultiPop gwasModel(ldFiles, zFiles, outputFileName, totalCausalSNP, NCP, rho, histFlag);
	gwasModel.run();
	gwasModel.finishUp();

	return 0;
}
