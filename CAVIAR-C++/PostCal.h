#ifndef POSTCAL_H
#define POSTCAL_H

#include <iostream>
#include <fstream>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

using namespace std;
using namespace arma;

void printGSLPrint(mat A, int row, int col);
 
class PostCal{

private:
	
	double gamma;		// the probability of SNP being causal
	double * postValues;	//the posterior value for each SNP being causal
	double * sigma;		//the LD matrix
	double * histValues;	//the probability of the number of causal SNPs, we make the histogram of the causal SNPs
	int snpCount;		//total number of variants (SNP) in a locus
	int maxCausalSNP;	//maximum number of causal variants to consider in a locus
	double sigmaDet;	//determinie of matrix

	double totalLikeLihoodLOG; //Compute the total log likelihood of all causal status (by likelihood we use prior)

	mat sigmaMatrix;
	mat invSigmaMatrix;
	mat statMatrix;
        mat statMatrixtTran;
	string * snpNames;

	double addlogSpace(double a, double b) {
		if (a == 0)
			return b;
		if (b == 0)
			return a;
		double base = max(a,b);
		if (base - min(a,b) > 700)
			return base;
		return(base + log(1+exp(min(a,b)-base)));
	}

public:

	PostCal(double * sigma, double * stat, int snpCount, int maxCausalSNP, string * snpNames, double gamma) {
		this->gamma = gamma;
		this->snpNames = snpNames;
		this-> snpCount = snpCount;
		this-> maxCausalSNP = maxCausalSNP;
                this->sigma = new double[snpCount * snpCount];
		this-> postValues = new double [snpCount];
		this-> histValues = new double [maxCausalSNP+1];              
	
		statMatrix                 = mat (snpCount, 1);
		statMatrixtTran            = mat (1, snpCount);
		sigmaMatrix         	   = mat (snpCount, snpCount);
	
		for(int i = 0; i < snpCount*snpCount; i++)
			this->sigma[i] = sigma[i];
		for(int i = 0; i < snpCount; i++)
                        this->postValues[i] = 0;
		for(int i= 0; i <= maxCausalSNP;i++)
			this->histValues[i] = 0;
		for(int i = 0; i < snpCount; i++) {
                	statMatrix(i,0) = stat[i];
        	        statMatrixtTran(0,i) = stat[i];
	        }
		
		for(int i = 0; i < snpCount; i++) {
                	for (int j = 0; j < snpCount; j++)
                       		sigmaMatrix(i,j) = sigma[i*snpCount+j];
       		}
		//invSigmaMatrix is depricated and the value for it is not right
		//PLASE DO NOT USE THE invSigmaMatrix;
		invSigmaMatrix = sigmaMatrix;
		sigmaDet       = det(sigmaMatrix);
	
	}
        ~PostCal() {
		delete [] histValues;
		delete [] postValues;
                delete [] sigma;
	}

	bool validConfigutation(int * configure, char * pcausalSet);
	void computeALLCausalSetConfiguration(double * stat, double NCP, char * pcausalSet, string outputFileName);


	double dmvnorm(mat Z, mat mean, mat R);
        double fracdmvnorm(mat Z, mat mean, mat R, mat diagC, double NCP);
	double fracdmvnorm2(mat Z, mat mean, mat R, mat diagC, double NCP);


        double fastLikelihood(int * configure, double * stat, double NCP);
	double likelihood(int * configure, double * stat, double NCP) ;
	int nextBinary(int * data, int size) ;
	double computeTotalLikelihood(double * stat, double NCP) ;	
	double findOptimalSetGreedy(double * stat, double NCP, char * pcausalSet, int *rank,  double inputRho, string outputFileName);
	string convertConfig2String(int * config, int size);
	void printHist2File(string fileName) {
		exportVector2File(fileName, histValues, maxCausalSNP+1);
	}

	void printPost2File(string fileName) {
		double total_post = 0;
		ofstream outfile(fileName.c_str(), ios::out );	
		for(int i = 0; i < snpCount; i++)
                	total_post = addlogSpace(total_post, postValues[i]);
		outfile << "SNP_ID\tProb_in_pCausalSet\tCausal_Post._Prob." << endl; 
		for(int i = 0; i < snpCount; i++) {
			outfile << snpNames[i] << "\t" << exp(postValues[i]-total_post) << "\t" << exp(postValues[i]-totalLikeLihoodLOG) << endl;
		}
	}

};
 
#endif
