#ifndef TOPKSNP_H
#define TOPKSNP_H

#include <cmath>

#include <vector>
#include <set>


using namespace std;

class TopKSNP{

private:
	double * stat;
	int snpCount;
public:
	TopKSNP(double * stat, int snpCount) {
	 	this-> snpCount = snpCount;
                this->stat = new double[snpCount];
                for(int i = 0; i < snpCount; i++)
                        this->stat[i] = stat[i];
        }
        ~TopKSNP() {
                delete [] stat;
        }
	void findCausal(int * topKConfigure);
};
 
#endif
