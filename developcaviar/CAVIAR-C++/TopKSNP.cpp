#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "Util.h"
#include "TopKSNP.h"


void TopKSNP::findCausal(int * topKConfigure) {
	std::vector<data> items;
	std::set<int> geneSet;
	std::set<int>::iterator it;
	
	for(int i = 0; i < snpCount; i++){
		items.push_back(data(stat[i], i, i));
	}
	std::sort(items.begin(), items.end(), by_number());
	printf("\n");
	for(int i = 0; i < snpCount; i++){
                printf("%f %d ", items[i].number, items[i].index1);
        	topKConfigure[i] = items[i].index1;
	}
	printf("\n");
}
