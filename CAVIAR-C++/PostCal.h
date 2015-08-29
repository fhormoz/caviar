#ifndef POSTCAL_H
#define POSTCAL_H

#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

void printGSLPrint(gsl_matrix *A, int row, int col);
 
class PostCal{

private:
	double * postValues;	//the posterior value for each SNP being causal
	double * sigma;		//the LD matrix
	double * histValues;	//the probability of the number of causal SNPs, we make the histogram of the causal SNPs
	int snpCount;		//total number of variants (SNP) in a locus
	int maxCausalSNP;	//maximum number of causal variants to consider in a locus

	gsl_matrix * sigmaMatrix;
	gsl_matrix * invSigmaMatrix;
	gsl_matrix ** matrixPossibleCases; 	//Compute the set of matrixes which repearsent the \sigma \sigma_c \sigma
	double baseValue;			//base value used for calculation for overflow	
public:

	PostCal(double * sigma, int snpCount, int maxCausalSNP) {
                this-> snpCount = snpCount;
		this-> maxCausalSNP = maxCausalSNP;
                this->sigma = new double[snpCount * snpCount];
		this-> postValues = new double [snpCount];
		this-> histValues = new double [maxCausalSNP+1];              
 
		for(int i = 0; i < snpCount*snpCount; i++)
			this->sigma[i] = sigma[i];
		for(int i = 0; i < snpCount; i++)
                        this->postValues[i] = 0;
		for(int i= 0; i <= maxCausalSNP;i++)
			this->histValues[i] = 0;
		baseValue = 0;
		matrixPossibleCases = new gsl_matrix*[snpCount];
		sigmaMatrix         = gsl_matrix_calloc (snpCount, snpCount);
		for(int i = 0; i < snpCount; i++) {
                	for (int j = 0; j < snpCount; j++)
                       		gsl_matrix_set(sigmaMatrix,i,j,sigma[i*snpCount+j]);
       		}
		for(int i = 0; i < snpCount; i++) {
			gsl_matrix * tmpMatrix1 = gsl_matrix_calloc (snpCount, snpCount);
			gsl_matrix * tmpMatrix2 = gsl_matrix_calloc (snpCount, snpCount);
			gsl_matrix * tmpMatrix3 = gsl_matrix_calloc (snpCount, snpCount);
			for(int j = 0; j < snpCount; j++){
				for(int k = 0; k < snpCount; k++) {	
					gsl_matrix_set(tmpMatrix1,j,k, 0);
					if(j==k && i == j)	gsl_matrix_set(tmpMatrix1,j,j, 1);
				}
			}
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, sigmaMatrix, tmpMatrix1, 0.0, tmpMatrix2);
		        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmpMatrix2, sigmaMatrix, 0.0, tmpMatrix3);	
			matrixPossibleCases[i] = tmpMatrix3;
			gsl_matrix_free(tmpMatrix1);
			gsl_matrix_free(tmpMatrix2);	
		}
		invSigmaMatrix = gsl_matrix_calloc (snpCount, snpCount);
		gsl_matrix_memcpy(invSigmaMatrix, sigmaMatrix);	
		// GOOD FOR Positive deifnite matrices
		//gsl_linalg_cholesky_decomp(invSigmaMatrix);
	        //gsl_linalg_cholesky_invert(invSigmaMatrix);
		//IF the matrix is not Positive deifnite
		int tmpS = 0;
		gsl_permutation * p = gsl_permutation_alloc (snpCount);
		gsl_linalg_LU_decomp(sigmaMatrix, p, &tmpS);
		gsl_linalg_LU_invert(sigmaMatrix, p, invSigmaMatrix);
		for(int i = 0; i < snpCount; i++) {
                        for (int j = 0; j < snpCount; j++)
                                gsl_matrix_set(sigmaMatrix,i,j,sigma[i*snpCount+j]);
                }
		gsl_permutation_free(p);	
	}
        ~PostCal() {
		delete [] histValues;
		delete [] postValues;
                delete [] sigma;
        	for(int i = 0 ; i < snpCount; i++)
			gsl_matrix_free(matrixPossibleCases[i]);
	}

	double likelihood(int * configure, double * stat, double NCP) ;
	int nextBinary(int * data, int size) ;
	double totalLikelihood(double * stat, double NCP) ;	
	double findOptimalSetGreedy(double * stat, double NCP, char * configure, int *rank,  double inputRho);
	string convertConfig2String(int * config, int size);
	void printHist2File(string fileName) {
		exportVector2File(fileName, histValues, maxCausalSNP+1);
	}

};
 
#endif
