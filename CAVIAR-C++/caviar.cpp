#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <unistd.h> 

#include "Util.h"
#include "PostCal.h"
#include "TopKSNP.h"
#include "CaviarModel.h"

using namespace std;


int main( int argc, char *argv[]  ){
	int totalCausalSNP = 2;
	double NCP = 5.2;
	double gamma = 0.01;
	double rho = 0.95;
	bool histFlag = false;
	int oc = 0;	
	string ldFile = "";
	string zFile  = "";
	string outputFileName = "";
	string geneMapFile = "";	

	while ((oc = getopt(argc, argv, "vhl:o:z:g:r:c:f:")) != -1) {
		switch (oc) {
			case 'v':
				cout << "version 2.2:" << endl;
			case 'h':
				cout << "Options: " << endl;
  				cout << "-h, --help            		show this help message and exit " << endl;
  				cout << "-o OUTFILE, --out=OUTFILE 	specify the output file" << endl;
  				cout << "-l LDFILE, --ld_file=LDFILE  	the ld input file" << endl;
  				cout << "-z ZFILE, --z_file=ZFILE	the z-score and rsID files" << endl;
  				cout << "-r RHO, --rho-prob=RHO		set $pho$ probability (default 0.95)" << endl;
				cout << "-g GAMMA, --gamma		set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
				cout << "-c causal			set the maximum number of causal SNPs" << endl;
				cout << "-f 1				to out the probaility of different number of causal SNP" << endl;
				exit(0);
			case 'l':
				ldFile = string(optarg);
				break;
			case 'o':
				outputFileName = string(optarg);
				break;
			case 'z':
				zFile = string(optarg);
				break;
			case 'r':
				rho = atof(optarg);
				break;
			case 'c':
				totalCausalSNP = atoi(optarg);
				break;
			case 'g':
				gamma = atof(optarg);
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
	cout << "| CAVIAR!		| 	   v2.2         |  10/Apr/2018| " << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| (C) 2018 Farhad Hormozdiari, GNU General Public License, v2 |" << endl;
	cout << "|-------------------------------------------------------------|" << endl;
	cout << "| For documentation, citation & bug-report instructions:      |" << endl;
	cout << "| 		http://genetics.cs.ucla.edu/caviar/           |" << endl;
	cout << "@-------------------------------------------------------------@" << endl;	

	CaviarModel caviar(ldFile, zFile, outputFileName, totalCausalSNP, NCP, rho, histFlag, gamma);
	caviar.run();
	caviar.finishUp();		
	return 0;
}
