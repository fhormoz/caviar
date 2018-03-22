#ifndef POSTCALMULTIPOP_H
#define POSTCALMULTIPOP_H

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
 
class PostCalMultiPop{

private:
	int popNum;
	double gamma;		// the probability of SNP being causal
	double * postValues;	//the posterior value for each SNP being causal
	double * histValues;	//the probability of the number of causal SNPs, we make the histogram of the causal SNPs
	int snpCount;		//total number of variants (SNP) in a locus
	int maxCausalSNP;	//maximum number of causal variants to consider in a locus
	double totalLikeLihoodLOG; //Compute the total likelihood of all causal status (by likelihood we use prior)

	vector <double *>  sigma;       //the list of all LD matrix
	vector <double *>  stat;	// the list of all stat
	vector <double> sigmaDet;	//the list of determinie of matrix

	vector <mat> sigmaMatrix;
	vector <mat> invSigmaMatrix;
	vector <mat> statMatrix;
        vector <mat> statMatrixtTran;
	vector <double> baseValue;			//base value used for calculation for overflow	
	string * snpNames;
public:

	PostCalMultiPop( vector <double *> sigma, vector <double *> stat, int snpCount, int maxCausalSNP, string * snpNames, double gamma) {
		this->gamma = gamma;
		this->snpNames = snpNames;
		this-> snpCount = snpCount;
		this-> maxCausalSNP = maxCausalSNP;
		this-> postValues = new double [snpCount];
		this-> histValues = new double [maxCausalSNP+1];              
		for(int i = 0; i < snpCount; i++)
			this->postValues[i] = 0;
		for(int i= 0; i <= maxCausalSNP;i++)
			this->histValues[i] = 0;
		for (int pop = 0; pop < sigma.size(); pop++) {
			baseValue.push_back(0);
			this->sigma.push_back(new double[snpCount * snpCount]);	
			this->stat.push_back(new double[snpCount]);
			statMatrix.push_back( mat (snpCount, 1,fill::zeros));
			statMatrixtTran.push_back(mat (1, snpCount,fill::zeros));
			sigmaMatrix.push_back(mat (snpCount, snpCount,fill::zeros));	
			for(int i = 0; i < snpCount*snpCount; i++)
				this->sigma[pop][i] = sigma[pop][i];
			for(int i = 0; i < snpCount; i++)
				this->stat[pop][i]  = stat[pop][i];
			for(int i = 0; i < snpCount; i++) {
				statMatrix[pop](i,0) = stat[pop][i];
				statMatrixtTran[pop](0,i) = stat[pop][i];
			}
			for(int i = 0; i < snpCount; i++) {
				for (int j = 0; j < snpCount; j++)
					sigmaMatrix[pop](i,j) = sigma[pop][i*snpCount+j];
			}
			invSigmaMatrix.push_back(inv(sigmaMatrix[pop]));
			sigmaDet.push_back(det(sigmaMatrix[pop]));
		}
		popNum = sigma.size();
	}
        ~PostCalMultiPop() {
		delete [] histValues;
		delete [] postValues;
		for (int pop = 0; pop < popNum; pop++)
                	delete [] sigma[pop];
	}

	double dmvnorm(mat Z, mat mean, mat R);
        double fracdmvnorm(mat Z, mat mean, mat R, mat diagC, int pop,double NCP);

        double fastLikelihood(int * configure, int pop, double NCP);
	int nextBinary(int * data, int size) ;
	double computeTotalLikelihood(double NCP) ;	
	double findOptimalSetGreedy(double NCP, char * pcausalSet, int *rank,  double inputRho, string outputFileName);
	string convertConfig2String(int * config, int size);
	void printHist2File(string fileName) {
		//exportVector2File(fileName, histValues, maxCausalSNP+1);
	}

	void printPost2File(string fileName) {
		double total_post = 0;
		ofstream outfile(fileName.c_str(), ios::out );	
		for(int i = 0; i < snpCount; i++)
                	total_post += postValues[i];
		outfile << "SNP_ID\tProb_in_pCausalSet\tCausal_Post._Prob." << endl; 
		for(int i = 0; i < snpCount; i++) {
			outfile << snpNames[i] << "\t" << postValues[i]/total_post << "\t" << postValues[i]/totalLikeLihoodLOG << endl;
		}
	}

};
 
#endif
