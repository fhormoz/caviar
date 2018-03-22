#ifndef CAVIARMODELMULTIPOP_H
#define CAVIARMODELMULTIPOP_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

#include "PostCalMultiPop.h"

using namespace std;
using namespace arma;

 
class CaviarModelMultiPop{
	
public:
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
        char * pcausalSet;
        int * rank;
        bool histFlag;
	PostCalMultiPop * post;
	string * snpNames;

	vector <double *> sigma;
        vector <double *> stat;
	vector <string> ldFile;
        vector <string> zFile;
        string outputFileName;

	CaviarModelMultiPop(vector <string> ldFile, vector <string>  zFile, string outputFileName, 
				int totalCausalSNP, double NCP, double rho, bool histFlag, double gamma=0.01) {
		int tmpSize = 0;
		this->histFlag = histFlag;
		this->NCP = NCP;
		this->rho = rho;
		this->gamma = gamma;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
		fileSize(ldFile[0], tmpSize);
	        snpCount   = (int)sqrt(tmpSize);
		pcausalSet = new char[snpCount];
                rank       = new int[snpCount];        
		snpNames   = new string [snpCount];
		importDataFirstColumn(zFile[0], snpNames);
		for (int pop = 0; pop < ldFile.size();pop++) {
			this->ldFile.push_back(ldFile[pop]);
			this->zFile.push_back(zFile[pop]);	
			sigma.push_back(new double[snpCount * snpCount]);
			stat.push_back(new double[snpCount]);
			importData(ldFile[pop], sigma[pop]);
			makeSigmaPositiveSemiDefinite(sigma[pop], snpCount);
			importDataSecondColumn(zFile[pop], stat[pop]);
		}
		post = new PostCalMultiPop(sigma, stat, snpCount, totalCausalSNP, snpNames, gamma);
	}
	void run() {
        	post->findOptimalSetGreedy(NCP, pcausalSet, rank, rho, outputFileName);
	}
	void finishUp() {
		ofstream outputFile;
                string outFileNameSet = string(outputFileName)+"_set";
                outputFile.open(outFileNameSet.c_str());
                for(int i = 0; i < snpCount; i++) {
                        if(pcausalSet[i] == '1')
                                outputFile << snpNames[i] << endl;
                }
                post->printPost2File(string(outputFileName)+"_post");
                //output the histogram data to file
                if(histFlag)
                	post->printHist2File(string(outputFileName)+"_hist");
	}

	~CaviarModelMultiPop() {
	}

};
 
#endif
