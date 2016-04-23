#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>

#include "Util.h"
#include "PostCal.h"

#define SMALL 0.001

using namespace arma;


void printGSLPrint(mat &A, int row, int col) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++)
			printf("%g ", A(i, j));
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
	int causalCount = 0;
	int index_C = 0;
        double matDet = 0;
	double res    = 0;

	for(int i = 0; i < snpCount; i++) 
		causalCount += configure[i];
	if(causalCount == 0){
		mat tmpResultMatrix1N = statMatrixtTran * invSigmaMatrix;
		mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
		res = tmpResultMatrix11(0,0);	
		baseValue = res;
		matDet = sigmaDet;
		res = res - baseValue;
		return( exp(-res/2)/sqrt(abs(matDet)) );
	}
	mat U(snpCount, causalCount, fill::zeros);
	mat V(causalCount, snpCount, fill::zeros);
	mat VU(causalCount, causalCount, fill::zeros);
		
	for(int i = 0; i < snpCount; i++) {
                if (configure[i] == 0)	continue;
                else {
                        for(int j = 0; j < snpCount; j++) 
                                U(j, index_C) = sigmaMatrix(j,i);
			V(index_C, i) = NCP;
                        index_C++;
                }
        }
	VU = V * U;
	mat I_AA   = mat(snpCount, snpCount, fill::eye);
	mat tmp_CC = mat(causalCount, causalCount, fill::eye)+ VU;
	matDet = det(tmp_CC) * sigmaDet;
	mat tmp_AA = invSigmaMatrix - (invSigmaMatrix * U) * pinv(tmp_CC) * V ;
	//tmp_AA     = invSigmaMatrix * tmp_AA;
	mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
        mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
        res = tmpResultMatrix11(0,0);  

	res = res - baseValue;
	if(matDet==0) {
		cout << "Error the matrix is singular and we fail to fix it." << endl;
		exit(0);
	}
	/*
		We compute the log of -res/2-log(det) to see if it is too big or not. 
		In the case it is too big we just make it a MAX value.
	*/
	double tmplogDet = log(sqrt(abs(matDet)));
	double tmpFinalRes = -res/2 - tmplogDet;
	if(tmpFinalRes > 700) 
		return(exp(700));
	return( exp(-res/2)/sqrt(abs(matDet)) );	
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

double PostCal::computeTotalLikelihood(double * stat, double NCP) {	
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
		for(int j = 0; j < snpCount; j++) {
                        postValues[j] = postValues[j] + tmp_likelihood * configure[j];
		}
		histValues[num] = histValues[num] + tmp_likelihood;
                num = nextBinary(configure, snpCount);
       		if(i % 100000 == 0)
			cout << i << " "  << sumLikelihood << endl;
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
        double total_post = 0;

        totalLikeLihood = computeTotalLikelihood(stat, NCP);

	for(int i = 0; i < snpCount; i++)
		total_post += postValues[i];

	printf("Total Likelihood= %e SNP=%d \n", total_post, snpCount);
	
        std::vector<data> items;
        std::set<int>::iterator it;
	//output the poster to files
        for(int i = 0; i < snpCount; i++) {
             //printf("%d==>%e ",i, postValues[i]/total_likelihood);
             items.push_back(data(postValues[i]/total_post, i, 0));
        }
        printf("\n");
        std::sort(items.begin(), items.end(), by_number());
        for(int i = 0; i < snpCount; i++)
                rank[i] = items[i].index1;

        for(int i = 0; i < snpCount; i++)
                configure[i] = '0';
        do{
                rho += postValues[rank[index]]/total_post;
                configure[rank[index]] = '1';
                printf("%d %e\n", rank[index], rho);
                index++;
        } while( rho < inputRho);

        printf("\n");
	return(0);
}
