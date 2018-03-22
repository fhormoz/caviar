#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>

#include "Util.h"
#include "PostCalMultiPop.h"

using namespace arma;


void printGSLPrint(mat &A, int row, int col) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++)
			printf("%g ", A(i, j));
		printf("\n");
	}	
}

string PostCalMultiPop::convertConfig2String(int * config, int size) {
	string result = "0";
	for(int i = 0; i < size; i++)
		if(config[i]==1)
			result+= "_" + convertInt(i);
	return result;
}

// We compute dmvnorm(Zcc, mean=rep(0,nrow(Rcc)), Rcc + Rcc %*% Rcc) / dmvnorm(Zcc, rep(0, nrow(Rcc)), Rcc))
// // togheter to avoid numerical over flow
double PostCalMultiPop::fracdmvnorm(mat Z, mat mean, mat R, mat diagC, int pop, double NCP) {
        mat newR = R + R * diagC  * R;
        mat ZcenterMean = Z - mean;
        mat res1 = trans(ZcenterMean) * inv(R) * (ZcenterMean);
        mat res2 = trans(ZcenterMean) * inv(newR) *  (ZcenterMean);

        double v1 = res1(0,0)/2-res2(0,0)/2-baseValue[pop]/2;
        return(exp(v1)/sqrt(det(newR))* sqrt(det(R)));
}

double PostCalMultiPop::dmvnorm(mat Z, mat mean, mat R) {
        mat ZcenterMean = Z - mean;
        mat res = trans(ZcenterMean) * inv(R) * (ZcenterMean);
        double v1 = res(0,0);
        double v2 = log(sqrt(det(R)));
        return (exp(-v1/2-v2));
}

// cc=causal SNPs
// Rcc = LD of causal SNPs
// Zcc = Z-score of causal SNPs
// dmvnorm(Zcc, mean=rep(0,nrow(Rcc)), Rcc + Rcc %*% Rcc) / dmvnorm(Zcc, rep(0, nrow(Rcc)), Rcc))
//
double PostCalMultiPop::fastLikelihood(int * configure, int pop, double NCP) {
	int causalCount = 0;
	vector <int> causalIndex;

	for(int i = 0; i < snpCount; i++) {
		causalCount += configure[i];
		if(configure[i] == 1)
			causalIndex.push_back(i);
	}
	
	if (causalCount == 0) {
		int maxVal = 0;
		for(int i = 0; i < snpCount; i++) {
			if (maxVal < abs(stat[pop][i]))
				maxVal = stat[pop][i];
		}
		baseValue[pop] = maxVal * maxVal;
	}

	mat Rcc(causalCount, causalCount, fill::zeros);
	mat Zcc(causalCount, 1, fill::zeros);
	mat mean(causalCount, 1, fill::zeros);
	mat diagC(causalCount, causalCount, fill::zeros);

	for (int i = 0; i < causalCount; i++){
		for(int j = 0; j < causalCount; j++) {
			Rcc(i,j) = sigmaMatrix[pop](causalIndex[i], causalIndex[j]);
		}
		Zcc(i,0) = stat[pop][causalIndex[i]];
		diagC(i,i) = NCP;
	}
	return fracdmvnorm(Zcc, mean, Rcc, diagC,pop, NCP);
}

int PostCalMultiPop::nextBinary(int * data, int size) {
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

double PostCalMultiPop::computeTotalLikelihood(double NCP) {	
	int num = 0;
	double sumLikelihood = 0;
	double tmp_likelihood = 1;
	long int total_iteration = 0 ;
	int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data	

	for(long int i = 0; i <= maxCausalSNP; i++)
		total_iteration = total_iteration + nCr(snpCount, i);
	cout << snpCount << endl;
	cout << "Max Causal=" << maxCausalSNP << endl;
	cout << "Total="      << total_iteration << endl;
	for(long int i = 0; i < snpCount; i++) 
		configure[i] = 0;
	for(long int i = 0; i < total_iteration; i++) {
		tmp_likelihood = 1;
		for (int pop =0 ; pop < popNum; pop++) 
                	tmp_likelihood = tmp_likelihood * fastLikelihood(configure, pop, NCP) * (pow(gamma, num))*(pow(1-gamma, snpCount-num));
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
double PostCalMultiPop::findOptimalSetGreedy(double NCP, char * pcausalSet, int *rank,  double inputRho, string outputFileName) {
	int index = 0;
        double rho = 0;
        double total_post = 0;

        totalLikeLihoodLOG = computeTotalLikelihood(NCP);
	
	export2File(outputFileName+".log", totalLikeLihoodLOG); //Output the total likelihood to the log File
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
                pcausalSet[i] = '0';
        do{
                rho += postValues[rank[index]]/total_post;
                pcausalSet[rank[index]] = '1';
                printf("%d %e\n", rank[index], rho);
                index++;
        } while( rho < inputRho);

        printf("\n");
	return(0);
}
