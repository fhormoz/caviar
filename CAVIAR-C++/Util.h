#ifndef UTIL_H
#define UTIL_H

#include <cmath>
#include <string>
#include <vector>

using namespace std;

	struct data {
                data(double num, int ind1, int ind2) {
                        number = num;
                        index1 = ind1;
                        index2 = ind2;
                }
                double number;
                int index1;
                int index2;
        };

        struct by_number {
            bool operator()(data const &left, data const &right) {
                return abs(left.number) > abs(right.number);
            }
        };



	string convertInt(int number);
	long int fact(int n) ;
	void copyConfigure(double *dest, double *src, int size) ;
	double min(double a, double b) ;
	long int nCr(int n, int r) ;
	void printVector(char * data, int size) ;
	void printVector(int * data, int size) ;
	void printVector(double * data, int size) ;
	void diffVector(double * data1, double * data2, int size, double * result) ;
	void sumVector(double * data1, double * data2, int size, double * result) ;
	double multVector(double * data1, double * data2, int size) ;
	void dotVector(double * data1, double * data2, int size, double * result) ;
	void multVectorMatrix(double *vector, double * matrix, int size, double * result) ;
	void fileSize(string fileName, int & size);
	void importData(string fileName, double * vector);
	void importData(string fileName, int * vector);
	void importDataSecondColumn(string fileName, double * vector);
	void importDataNthColumn(string fileName, double * vector, int colNum);
	void importDataFirstColumn(string fileName, string * list);
	void rmvnorm(double * mean, double * sigma, int size, double * results);
	void resetVector(char *data, int size);
	void resetVector(int *data, int size);
	void resetVector(double *data, int size);
	void generateMean(int * causalSNP, double * sigma, int snpCount, double * result);
	void exportVector2File(string fileName, char * data, int size);
	void exportVector2File(string fileName, int * data, int size);
	void exportVector2File(string fileName, double * data, int size);
	void export2File(string fileName, int data);
	void matrixMul(int * aData, int * bData, int * cData, int row1, int col1, int row2, int col2);
	int snp2Gene(int * G, int snpId, int snpCount, int geneCount);
	void setIdentitymatrix(int * G, int snpCount, int geneCount);
	void makeSigmaPositiveSemiDefinite(double * sigma, int size);


#endif

