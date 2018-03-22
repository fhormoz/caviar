#ifndef CAVIARMODEL_H
#define CAVIARMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

#include "PostCal.h"

using namespace std;
using namespace arma;

 
class CaviarModel{
	
public:
	double rho;
	double NCP;
	double gamma;
	int snpCount;
	int totalCausalSNP;
	double * sigma;
        double * stat;
        char * pcausalSet;
        int * rank;
        bool histFlag;
	PostCal * post;
	string * snpNames;
	string ldFile;
        string zFile;
        string outputFileName;
        string geneMapFile;	

	CaviarModel(string ldFile, string zFile, string outputFileName, int totalCausalSNP, double NCP, double rho, bool histFlag, double gamma=0.01) {
		int tmpSize = 0;
		this->histFlag = histFlag;
		this->NCP = NCP;
		this->rho = rho;
		this->gamma = gamma;
		this->ldFile = ldFile;
		this->zFile  = zFile;
		this->outputFileName = outputFileName;
		this->totalCausalSNP = totalCausalSNP;
		fileSize(ldFile, tmpSize);
	        snpCount   = (int)sqrt(tmpSize);
         	sigma      = new double[snpCount * snpCount];
		stat       = new double[snpCount];
		pcausalSet = new char[snpCount];
		rank       = new int[snpCount];
		snpNames   = new string [snpCount];
		importData(ldFile, sigma);
		importDataFirstColumn(zFile, snpNames);
		importDataSecondColumn(zFile, stat);
		makeSigmaPositiveSemiDefinite(sigma, snpCount);
		for (int i = 0; i < snpCount; i++){
			if (abs(stat[i]) > NCP) NCP = abs(stat[i]);
		}
		post = new PostCal(sigma, stat, snpCount, totalCausalSNP, snpNames, gamma);
	}
	void run() {
        	post->findOptimalSetGreedy(stat, NCP, pcausalSet, rank, rho, outputFileName);
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

	void printLogData() {
		//print likelihood
		//print all possible configuration from the p-causal set
		post->computeALLCausalSetConfiguration(stat, NCP, pcausalSet,outputFileName+".log");
	}

	~CaviarModel() {
	}

};
 
#endif
