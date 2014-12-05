#include <vector>
#include <set>
#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "Util.h"
#include "PostCal.h"

#define SMALL 0.001

void printGSLPrint(gsl_matrix *A, int row, int col) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++)
			printf("%g ", gsl_matrix_get(A, i, j));
		printf("\n");
	}	
}

string PostCal::convertConfig2String(int * config, int size) {
	string result = "0";
	for(int i = 0; i < size; i++)
		if(config[i]==1)
			result+= "_" + convertInt(i);
	return result;
}

double PostCal::likelihood(int * configure, double * stat, double NCP) {
	int coutOne = 0;
	int gsl_tmp = 0;
	double val = 0;
        double matDet = 0;
	double res    = 0;

	gsl_matrix * statMatrix                 = gsl_matrix_calloc (snpCount, 1);
	gsl_matrix * statMatrixtTran            = gsl_matrix_calloc (1, snpCount);
	gsl_matrix * configMatrix		= gsl_matrix_calloc (snpCount, 1);

	gsl_matrix * tmpResultMatrix		= gsl_matrix_calloc (snpCount, 1);
	gsl_matrix * tmpResultMatrix1 		= gsl_matrix_calloc (snpCount, snpCount);
	gsl_matrix * tmpResultMatrix2           = gsl_matrix_calloc (snpCount, snpCount);
	gsl_matrix * tmpResultMatrix1N		= gsl_matrix_calloc (1, snpCount);
	gsl_matrix * tmpResultMatrix11          = gsl_matrix_calloc (1, 1);

        for(int i = 0; i < snpCount; i++) {
                gsl_matrix_set(statMatrix,i,0,stat[i]);
                gsl_matrix_set(statMatrixtTran,0,i,stat[i]);
        	coutOne += configure[i];
	}
	gsl_matrix_memcpy(tmpResultMatrix1, sigmaMatrix);
	for(int i = 0; i < snpCount; i++) {
		if(configure[i] == 1){
			gsl_matrix_memcpy(tmpResultMatrix2, matrixPossibleCases[i]); 
			gsl_matrix_scale(tmpResultMatrix2, NCP);
			gsl_matrix_add(tmpResultMatrix1, tmpResultMatrix2);
		}
	}

	gsl_matrix_memcpy(tmpResultMatrix2, tmpResultMatrix1);
	gsl_permutation *p = gsl_permutation_alloc(snpCount);
        gsl_linalg_LU_decomp(tmpResultMatrix2, p, &gsl_tmp );
        matDet = gsl_linalg_LU_det(tmpResultMatrix2,gsl_tmp);
	gsl_linalg_cholesky_decomp(tmpResultMatrix1);
        gsl_linalg_cholesky_invert(tmpResultMatrix1);	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, statMatrixtTran, tmpResultMatrix1, 0.0, tmpResultMatrix1N);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmpResultMatrix1N, statMatrix, 0.0, tmpResultMatrix11);	

	res = gsl_matrix_get(tmpResultMatrix11,0,0);

	gsl_matrix_free(tmpResultMatrix);
	gsl_matrix_free(tmpResultMatrix1);
	gsl_matrix_free(tmpResultMatrix2);
	gsl_matrix_free(tmpResultMatrix1N);
	gsl_matrix_free(tmpResultMatrix11);
	gsl_matrix_free(statMatrix);
	gsl_matrix_free(statMatrixtTran);
	gsl_matrix_free(configMatrix);
	gsl_permutation_free(p);

	if(baseValue == 0)
		baseValue = res;
	res = res - baseValue;
	
	return( exp(-res/2)/sqrt(matDet) );	
}

int PostCal::nextBinary(int * data, int size) {
	int i = 0;
	int total_one = 0;	
	int index = size-1;
        int one_countinus_in_end = 0;

        while(index >= 0 && data[index] == 1) {
                index = index - 1;
                one_countinus_in_end = one_countinus_in_end + 1;
	}
	if(index >= 0) {
        	while(index >= 0 && data[index] == 0) {
               	 index = index - 1;	
		}
	}
        if(index == -1) {
                while(i <  one_countinus_in_end+1 && i < size) {
                        data[i] = 1;
                        i=i+1;
		}
                i = 0;
                while(i < size-one_countinus_in_end-1) {
                        data[i+one_countinus_in_end+1] = 0;
                        i=i+1;
		}
	}
        else if(one_countinus_in_end == 0) {
                data[index] = 0;
                data[index+1] = 1;
	} else {
                data[index] = 0;
                while(i < one_countinus_in_end + 1) {
                        data[i+index+1] = 1;
			if(i+index+1 >= size)
				printf("ERROR3 %d\n", i+index+1);
                        i=i+1;
		}
                i = 0;
                while(i < size - index - one_countinus_in_end - 2) {
                        data[i+index+one_countinus_in_end+2] = 0;
			if(i+index+one_countinus_in_end+2 >= size) {
				printf("ERROR4 %d\n", i+index+one_countinus_in_end+2);
			}
                        i=i+1;
		}
	}
	i = 0;
	total_one = 0;
	for(i = 0; i < size; i++)
		if(data[i] == 1)
			total_one = total_one + 1;
	
	return(total_one);		
}

double PostCal::totalLikelihood(double * stat, double NCP) {	
	int num = 0;
	double sumLikelihood = 0;
	double tmp_likelihood = 0;
	long int total_iteration = 0 ;
	int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data	

	for(long int i = 0; i <= maxCausalSNP; i++)
		total_iteration = total_iteration + nCr(snpCount, i);

	for(long int i = 0; i < snpCount; i++) 
		configure[i] = 0;
	for(long int i = 0; i < total_iteration; i++) {
                tmp_likelihood = likelihood(configure, stat, NCP) * (pow(SMALL, num))*(pow(1-SMALL, snpCount-num));
                sumLikelihood += tmp_likelihood;
		//printVector(configure, snpCount);
		//printf("\t%e\n", tmp_likelihood);
		for(int j = 0; j < snpCount; j++)
                        postValues[j] = postValues[j] + tmp_likelihood *configure[j];
		histValues[num] = histValues[num] + tmp_likelihood;
                num = nextBinary(configure, snpCount);
       		if(i % 10000 == 0)
			cout << i << endl;
	}
	for(int i = 0; i <= maxCausalSNP; i++)
		histValues[i] = histValues[i]/sumLikelihood;
        free(configure);
        return(sumLikelihood);
}

/*
	stat is the z-scpres
	sigma is the correaltion matrix
	G is the map between snp and the gene (snp, gene)
*/
double PostCal::findOptimalSetGreedy(double * stat, double NCP, char * configure, int *rank,  double inputRho) {
	int index = 0;
        double rho = 0;
        double total_likelihood = 0;

        totalLikelihood(stat, NCP);

        printf("\nRho = %lf\n", inputRho);
        printf("Total = %e\n", total_likelihood);

	for(int i = 0; i < snpCount; i++)
		total_likelihood += postValues[i];
	
        std::vector<data> items;
        std::set<int>::iterator it;
        for(int i = 0; i < snpCount; i++) {
             printf("%d==>%e ",i, postValues[i]/total_likelihood);
             items.push_back(data(postValues[i]/total_likelihood, i, 0));
        }
        printf("\n");
        std::sort(items.begin(), items.end(), by_number());
        for(int i = 0; i < snpCount; i++)
                rank[i] = items[i].index1;

        for(int i = 0; i < snpCount; i++)
                configure[i] = '0';
        do{
                rho += postValues[rank[index]]/total_likelihood;
                configure[rank[index]] = '1';
                printf("%d %e\n", rank[index], rho);
                index++;
        } while( rho < inputRho);

        printf("\n");

}
