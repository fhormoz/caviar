#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;


long int fact(int n) {
        if(n==0)
                return 1;
        return n* fact(n-1);
}

void copyConfigure(double *dest, double *src, int size) {
	for(int i = 0; i < size; i++) 
		dest[i] = src[i];
}

double min(double a, double b) {
	if(a>b)
		return b;
	else
		return a;
}

long int nCr(int n, int r) {
        long int result = 1;
        for(int i = n; i > n-r; i--)
                result *= i;
        return result/fact(r);
}

void printVector(char * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%c, ", data[i]);
}

void printVector(int * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%d, ", (int)data[i]);
}

void printVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                printf("%lf, ", data[i]);
}

void diffVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
		result[i] = data1[i] - data2[i];
}

void sumVector(double * data1, double * data2, int size, double * result) {
        for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] + data2[i];
}

double multVector(double * data1, double * data2, int size) {
	double res = 0;
	for(int i = 0; i < size; i++ ) 
                res += data1[i] * data2[i];
	return res;
}

void dotVector(double * data1, double * data2, int size, double * result) {
	for(int i = 0; i < size; i++ ) 
                result[i] = data1[i] * data2[i];
}

void multVectorMatrix(double *vector, double * matrix, int size, double * result) {
	double total_row = 0;
	for(int i = 0; i < size; i++) {
		total_row = 0;
		for(int j = 0; j < size; j++) {
			total_row += vector[j] * matrix[i + j * size];
		}
		result[i]= total_row;
	}
}

void importData(string fileName, double * vector) {
	int index = 0;
	double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        fin >> data;
        while (fin.good()) {
                vector[index] = data;
                index++;
                fin >> data;
        }
        fin.close();
}

void importData(string fileName, int * vector) {
        int index = 0;
        double data = 0;
        ifstream fin(fileName.c_str(), std::ifstream::in);
        while( fin.good()  ){
                fin >> data;
                vector[index] = (int)data;
                index++;
        }
        fin.close();
}
/*
the column index starts by 1 in this implemenation
*/
void importDataSecondColumn(string fileName, double * vector) {
	int index = 0;
	string line = "";
	string dataS = "";
	double data = 0.0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( getline(fin, line) ){
		istringstream iss(line);
		iss >> dataS;
		iss >> data;
	        vector[index] = (double)data;
                index++;
        }
	cout << "reach=" << index << endl;
        fin.close();
}

void importDataFirstColumn(string fileName, string * list) {
 	int index = 0;
        string data = "";
        string line = "";
	ifstream fin(fileName.c_str(), std::ifstream::in);
        while( getline(fin, line) ){
		istringstream iss(line);
                iss >> data;
                list[index] = data;
		index++;
        }
        fin.close();
}

void fileSize(string fileName, int & size) {
	double data = 0;	
	ifstream fin(fileName.c_str(), std::ifstream::in);
	while( fin.good()  ){
		fin >> data;
		size++;
	}
	fin.close();
}

string convertInt(int number) {
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void rmvnorm(double * mean, double * sigma, int size, double * results) {
	gsl_matrix * work = gsl_matrix_alloc(size,size);
	gsl_vector * resultVector = gsl_vector_alloc(size);
	gsl_vector * meanVector  = gsl_vector_alloc(size);
	gsl_matrix * sigmaMatrix = gsl_matrix_alloc(size,size);
	const gsl_rng_type * T;
  	gsl_rng * r;	
	gsl_rng_env_setup();
	T = gsl_rng_default;
  	r = gsl_rng_alloc (T);	
	
	for(int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			gsl_matrix_set(sigmaMatrix,i,j,sigma[i*size + j]);
		}
		gsl_vector_set(meanVector,i ,mean[i]);
	}	

	gsl_matrix_memcpy(work, sigmaMatrix);
	gsl_linalg_cholesky_decomp(work);
	for(int i=0; i<size; i++)
		gsl_vector_set( resultVector, i, gsl_ran_ugaussian(r) );

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, resultVector);
	gsl_vector_add(resultVector,meanVector);
	gsl_matrix_free(work);
	
	for(int i = 0; i < size; i++)
		results[i] = gsl_vector_get(resultVector, i);
	
	//gsl_matrix_free(work);
	//gsl_vector_free(resultVector);
	//gsl_vector_free(meanVector);
	gsl_rng_free(r);
}

void generateMean(int * causalSNP, double * sigma, int size, double * result) {
	gsl_matrix * resultMatrix = gsl_matrix_alloc(1,size);
	gsl_matrix * sigmaMatrix = gsl_matrix_alloc(size,size);
	gsl_matrix * causalSNPMatrix = gsl_matrix_alloc(1,size);
	for(int i = 0; i < size; i++) {
                for(int j = 0; j < size; j++) {
                        gsl_matrix_set(sigmaMatrix,i,j,sigma[i*size + j]);
                }
		gsl_matrix_set(causalSNPMatrix, 0, i, 5.7 * causalSNP[i]);
        }
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, causalSNPMatrix, sigmaMatrix, 0.0, resultMatrix);
	for(int i = 0; i < size; i++)                   
                result[i] = gsl_matrix_get(resultMatrix, 0, i);

	gsl_matrix_free(sigmaMatrix); 
	gsl_matrix_free(resultMatrix);
	gsl_matrix_free(causalSNPMatrix); 
}

void resetVector(char *data, int size){
	for(int i = 0; i < size; i++)
		data[i] = '0';
}

void resetVector(int * data, int size) {
	for(int i = 0; i < size; i++)
		data[i] = 0;
}

void resetVector(double * data, int size) {
        for(int i = 0; i < size; i++)
                data[i] = 0;
}

void exportVector2File(string fileName, char * data, int size) {
	ofstream outfile(fileName.c_str(), ios::out | ios::app);
	for (int i = 0; i < size; i++)
		outfile << data[i] << " ";
	outfile << endl;
	outfile.close();
}

void exportVector2File(string fileName, double * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        outfile << endl;
        outfile.close();
}

void exportVector2File(string fileName, int * data, int size) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        for (int i = 0; i < size; i++)
                outfile << data[i] << " ";
        outfile << endl;
        outfile.close();
}

void export2File(string fileName, int data) {
        ofstream outfile(fileName.c_str(), ios::out | ios::app);
        outfile << data << endl;
        outfile.close();
}

// cData = aData * bData
void matrixMul(int * aData, int * bData, int * cData, int row1, int col1, int row2, int col2) {
	gsl_matrix * aMatrix = gsl_matrix_alloc(row1,col1);
        gsl_matrix * bMatrix = gsl_matrix_alloc(row2,col2);
	gsl_matrix * cMatrix = gsl_matrix_alloc(row1,col2);
	if(col1 == row2) {
		for(int i = 0; i < row1; i++)
			for(int j = 0; j < col1; j++)
				gsl_matrix_set(aMatrix,i,j,aData[i*col1 + j]);
		for(int i = 0; i < row2; i++) 
                        for(int j = 0; j < col2; j++)
                                gsl_matrix_set(bMatrix,i,j,bData[i*col2 + j]);	
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, aMatrix, bMatrix, 0.0, cMatrix);
		for(int i = 0; i < row1; i++)
			for(int j = 0; j < col2; j++) 
				cData[i * col2 + j] =  (int)gsl_matrix_get(cMatrix, i, j);
	}	
	gsl_matrix_free(aMatrix);
	gsl_matrix_free(bMatrix);
	gsl_matrix_free(cMatrix);
}

int snp2Gene(int * G, int snpId, int snpCount, int geneCount) {
	for(int i = 0; i < geneCount; i++) {
		if(G[snpId*geneCount + i] == 1)
			return i;
	}
	return -1;
}

void setIdentitymatrix(int * G, int snpCount, int geneCount) {
	for(int i = 0; i < snpCount; i++) {
		for(int j = 0; j < geneCount; j++) {
			G[i*geneCount + j] = 0;
		}
		G[i*geneCount + (i/(snpCount/geneCount))] = 1;
	}
	  for(int i = 0; i < snpCount; i++) {
                for(int j = 0; j < geneCount; j++) {
                        printf("%d ", G[i*geneCount+j]);
                }
		printf("\n");
        }
}
