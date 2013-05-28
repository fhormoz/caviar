import scipy.stats as st
import numpy as npy, numpy.linalg
from math import *
import random 
import csv
import sys
from math import factorial


def _getAplus(A):
    eigval, eigvec = npy.linalg.eig(A)
    Q = npy.matrix(eigvec)
    xdiag = npy.matrix(npy.diag(npy.maximum(eigval, 0)))
    return Q*xdiag*Q.T

def _getPs(A, W=None):
    W05 = npy.matrix(W**.5)
    return  W05.I * _getAplus(W05 * A * W05) * W05.I

def _getPu(A, W=None):
    Aret = npy.array(A.copy())
    Aret[W > 0] = npy.array(W)[W > 0]
    return npy.matrix(Aret)

def nearPD(A, nit=10):
    n = A.shape[0]
    W = npy.identity(n) 
# W is the matrix used for the norm (assumed to be Identity matrix here)
# the algorithm should work for any diagonal W
    deltaS = 0
    Yk = A.copy()
    for k in range(nit):
        Rk = Yk - deltaS
        Xk = _getPs(Rk, W=W)
        deltaS = Xk - Rk
        Yk = _getPu(Xk, W=W)
    return Yk

def write2File(outputFile, list) :
	for data in list:
		outputFile.write(str(data));
		outputFile.write(' ');


#calculate the power for give lambda and N (number of individual), lambda* sqrt(N) is the NCP
#power ~ Norm (lambda* sqrt(N), 1)
# x the value for which we compute the power
def calculate_power( in_lambda , N, x) :
	return (1-st.norm.cdf(x, loc=in_lambda*sqrt(N), scale=1)+ st.norm.cdf(-x, loc=in_lambda*sqrt(N), scale=1));


# Given the number of individual it will find the lambda to obtain the need power
# alpha is the significant level threshold

def calculate_cutoff_power(N, power, alpha) :
	DELTA = 0.000001;

        sAlpha = -st.norm.ppf( 0.5*alpha, loc=0, scale=1);

        startLambda = 0.0;
        endLambda   = 100.0;

        middleLambda = (startLambda+endLambda) / 2;

        while (abs(middleLambda-startLambda) > DELTA ) :
                if (calculate_power(middleLambda, N , sAlpha) == power) :
                        return (middleLambda);
                elif(calculate_power(middleLambda, N , sAlpha) < power) :
                        startLambda = middleLambda;
                elif(calculate_power(middleLambda, N , sAlpha) > power) :
                        endLambda   = middleLambda;
                middleLambda = (startLambda + endLambda)/2;
        
        return (middleLambda);


def pvalues(stat_value) :
	return(2.0*st.norm.cdf(-abs(stat_value), loc=0, scale=1));


def dmvnorm(Z, mu, R, Rinv, Rdet) :
	[row, col] = Z.shape;
	return( ( 1 / npy.sqrt(pow(2*pi, col) * abs(Rdet)) )   * exp ( -0.5 * (Z-mu) *  Rdet * (Z-mu).T  ) )

def likelihood_est(configure, Z, R, Rinv, Rdet, NCP) :
	#mu is the mean for the Multivariance Normal
	mu =(NCP * configure) * R;
	return(dmvnorm(Z, mu, R, Rinv, Rdet));


def prior_cal(prob_casual, total_snp, num_casual, casual_cutoff) :
	if(num_casual <= casual_cutoff) :
		return(pow(prob_casual,num_casual)*pow(1-prob_casual,  total_snp-num_casual));
	else :
		return(0)

def next_binary_value2(curr_binary_value) :
        [row, col] = curr_binary_value.shape;

        size = col-1;

        return_value =  curr_binary_value;

        #ignore the end in the end
        index = size;
        one_countinus_in_end = 0

        while(index >= 0 and return_value[0, index] == 1) :
                index = index - 1;
                one_countinus_in_end = one_countinus_in_end + 1;


        #ignore the zeros till you reach the first 1
        while(index >= 0 and return_value[0, index] == 0) :
                index = index - 1;


        first_one_in_begin = index;

        index1 = 0;
        while(index1 <= size and return_value[0, index1] == 0) :
                index1 = index1 + 1;

        if(index == -1) :
                for i in range(one_countinus_in_end+1):
                        return_value[0,i] = 1;
                for i in range(size-one_countinus_in_end):
                        return_value[0,i+one_countinus_in_end+1] = 0;
        elif(one_countinus_in_end == 0) :
                return_value[0,index] = 0;
                return_value[0,index+1] = 1;

        else :
                return_value[0, index] = 0;
                for i in range(one_countinus_in_end+1) :
                        return_value[0,i+index+1] = 1;
                for i in range(size - index - one_countinus_in_end-1) :
                        return_value[0, i+index+one_countinus_in_end+2] = 0;
        return(return_value);

def next_binary_value1(curr_binary_value) :
	[row, col] = curr_binary_value.shape;

	size = col-1;

        return_value =  curr_binary_value;

        #ignore the end in the end
        index = size;
        one_countinus_in_end = 0
        
	while(index >= 0 and return_value[0, index] == 1) :
                index = index - 1;
                one_countinus_in_end = one_countinus_in_end + 1;
       

        #ignore the zeros till you reach the first 1
        while(index >= 0 and return_value[0, index] == 0) :
                index = index - 1;
     

        first_one_in_begin = index;

        index1 = 0;
        while(index1 <= size and return_value[0, index1] == 0) :
                index1 = index1 + 1;
      
 
        if (first_one_in_begin > 0 and index1 != index ) :
                return_value[0, first_one_in_begin] = 0;
                return_value[0, first_one_in_begin+1] = 1;
        elif(index1 == index1) :
		tmp1 = npy.matrix(npy.zeros(shape = (1, index1+1)));
		tmp2 = npy.matrix(npy.zeros(shape = (1, one_countinus_in_end+1))+1);
                tmp3 = npy.matrix(npy.zeros(shape = (1, size - (index1+one_countinus_in_end+1))));
                return_value = npy.concatenate((tmp1,tmp2,tmp3),1)
        else : 
		tmp1 = npy.matrix(npy.zeros(shape = (1, size-one_countinus_in_end)));
		tmp2 = npy.matrix(npy.zeros(shape = (1, one_countinus_in_end+1))+1);
                return_value = npy.concatenate((tmp2,tmp1),1)
	
        return(return_value);

# Generate All the possible subset of a set 
# The ouput is a vector of 1 and 0, where 1
# mean the element exist in the subset and 0
# otherwise
def subsets(set, N) :

        size = len(set);
        #intialize the zero for the container of subsets
        res_subsets = npy.matrix(npy.zeros(shape = (pow(2,size), N)));

        index = 0;
        while(index < pow(2,size)) :
                tmp = decimal2binary(index, size);
		[row, col] = tmp.shape;
                index2 = 0;
		while(index2 < col) :
                        if(tmp[0,index2] == 1) :
                        	res_subsets[index,set[index2]] = 1; 
                        index2 = index2 + 1;
               
                index = index + 1;
        
        return(res_subsets);

def nCr(n,r):
    return factorial(n) / factorial(r) / factorial(n-r);

def subsets_cutoff(set, N, cut_off) :
	size = len(set);
	
	index = 0;
	matrix_size = 0;
	
	if(size == 0) :
		return npy.matrix(npy.zeros(shape = (1, N)));

	if(size >= cut_off) :	
		while(index <= cut_off) :
			matrix_size = matrix_size + nCr(size, index);	
			index = index + 1;
	else :
		matrix_size = pow(2,size);	

        #intialize the zero for the container of subsets
        res_subsets = npy.matrix(npy.zeros(shape = (matrix_size, N)));

        index = 0;
	cur_binary_value = npy.matrix(npy.zeros(shape = (1, size)));
        while(index < matrix_size-1) :
                index2 = 0;
                while(index2 < size) :
                        if(cur_binary_value[0,index2] == 1) :
                                res_subsets[index,set[index2]] = 1;
                        index2 = index2 + 1;
		index = index  + 1;
		if (size == 1):
			cur_binary_value[0,0] = 1;
		elif(index < matrix_size) :
			cur_binary_value = next_binary_value2(cur_binary_value);

        return(res_subsets);

# Help function to generate all the possible subset 
# this function will convert the decimal to binary
# format
def decimal2binary(num, size) :
        out = npy.zeros(shape = (1, size));
        tmp = num;
        index = 0;
        while(tmp > 0):
                out[0,index] = tmp % 2;
                tmp = tmp / 2;
                index = index + 1;
        return(out[::-1]);


def all_possible_perm_cutoff(configure, cut_off) :
        set = list();
        index1 = 0;
        [row, col] =  configure.shape;

        N = col-1;

        while(index1 <= N) :
                if(configure[0, index1]==1) :
                        set.append(index1);
                index1 = index1 + 1;

        sets = subsets_cutoff(set, N+1, cut_off);
        return (sets)


def all_possible_perm(configure) :
        set = list();
        index1 = 0;
        [row, col] =  configure.shape;

        N = col-1;

	while(index1 <= N) :
                if(configure[0, index1]==1) :
                        set.append(index1);
                index1 = index1 + 1;
     
        sets = subsets(set, N+1);
        return (sets)



def find_optimal_merge_memory(potential_casual, Z, R, Rinv, Rdet, power, alpha,ro, prob_casual) :
       
	step = 0; 
	total_sum = 0;
        total_casual_cut = 5;
        N = len(potential_casual);
        NCP = calculate_cutoff_power(100, power, alpha) * sqrt(100);

	max_value_likelihood = 0;
	max_value_likelihood_config = npy.matrix(npy.zeros(shape = (1, N), dtype=npy.int8));

        print 'NCP=', NCP;

	curr_binary_value = npy.matrix(npy.zeros(shape = (1, N), dtype=npy.int8));
        while(curr_binary_value.sum() <= total_casual_cut) :
                num1 = curr_binary_value.sum();
                prior = prior_cal(prob_casual, N, num1, total_casual_cut);
		tmp_likelihood = likelihood_est(curr_binary_value, Z, R, Rinv, Rdet, NCP);
		total_sum  = total_sum + tmp_likelihood*prior;
                if(tmp_likelihood > max_value_likelihood) :
			max_value_likelihood = tmp_likelihood;
			max_value_likelihood_config = curr_binary_value.copy();
		#print curr_binary_value, '\t', likelihood_est(curr_binary_value, Z, R, NCP), '\t' ,prior;
                #print step;
		step = step + 1;
		curr_binary_value = next_binary_value2(curr_binary_value).copy();

        print 'sum=', total_sum;
	print 'max=', max_value_likelihood;
	print 'max_conf=', max_value_likelihood_config; 

        index2 = 0;
	curr_binary_value = npy.matrix(npy.zeros(shape = (1, N)));
	#while(curr_binary_value.sum() <= total_casual_cut) :
        #        all_configure = all_possible_perm_cutoff(curr_binary_value, total_casual_cut);
        #        [row, col] = all_configure.shape;
	#	index2 = 0;
        #        tmpvalue = 0;
        #        while(index2 < row) :
        #                num1  = all_configure[index2,:].sum();
        #                prior = prior_cal(prob_casual, N, num1, total_casual_cut);
        #                tmpvalue = tmpvalue + likelihood_est(all_configure[index2,:], Z, R, Rinv, Rdet, NCP) * prior;
        #                index2 = index2 + 1;
        #        if(tmpvalue/total_sum >= ro) :
        #                return(curr_binary_value);
        #        curr_binary_value = next_binary_value2(curr_binary_value);
	print "WE USE HURISTIC FOR HERE"
	#return(max_value_likelihood_config);
	index1 = 0;
	index2 = 0;
	Z2list = Z.tolist()[0];	# we convert the Z matrix to list to use
       	index_Z_sorted = sorted(range(N), key=lambda k: abs(Z2list[k]), reverse=True); 
	print index_Z_sorted;
	curr_binary_value = npy.matrix(npy.zeros(shape = (1, N)));
	tmpvalue = 0;
	while(tmpvalue / total_sum < ro and curr_binary_value.sum() <= N*0.3) :
		tmpvalue = 0;
		curr_binary_value[0, index_Z_sorted[index1]]  = 1;
		all_configure = all_possible_perm_cutoff(curr_binary_value, total_casual_cut);
                [row, col] = all_configure.shape;
                index2 = 0;
                while(index2 < row) :
                        num1  = all_configure[index2,:].sum();
                        prior = prior_cal(prob_casual, N, num1, total_casual_cut);
                        tmpvalue = tmpvalue + likelihood_est(all_configure[index2,:], Z, R, Rinv, Rdet, NCP) * prior;
                        if(tmpvalue/total_sum >= ro) :
				return curr_binary_value;
			#print "step2=", all_configure[index2,:],'\t',tmpvalue, '\t', likelihood_est(all_configure[index2,:], Z, R, NCP) , '\t', prior;
                        index2 = index2 + 1;	
		index1 = index1 + 1;			
	return(curr_binary_value);



def main(argv):
	
	
	row = 1;
	col = 1;
	ro  = 0.95;
	power = 0.5;
	alpha = 1e-08;
	prob_casual = 0.00001;	
	casual_implant = 2;
	total_snp = 0;	

	NCP = calculate_cutoff_power(100, power, alpha) * sqrt(100);

	print NCP;

	#oldR = npy.random.random_sample(size=(row,row));
	#R1 = (oldR + oldR.T)/2;
	#R = (R1 * R1.T);
	#R = nearPD(R, 10);

	inputFile = open('/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld', 'rb');
	for data in inputFile:
		total_snp = total_snp + 1;
	R = npy.matrix(npy.zeros(shape = (2*total_snp, 2*total_snp)));	
	inputFile.close();	

	row = 0;
	col = 0;

	inputFile = open('/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld', 'rb');
	for data in inputFile:
        	line = data.split(' ');
		col = 0;
		for l in line:
			if(col < total_snp):
				R[row, col] = l;
			col = col + 1;
		#row = row + 1;
		col = total_snp;
		for l in line:
			if(col < 2 * total_snp) :
				R[row, col] = abs(float(l))/2;
			col = col + 1;
		row = row + 1; 		
	inputFile.close();

	inputFile = open('/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld', 'rb');
        for data in inputFile:
                line = data.split(' ');
                col = 0;
                for l in line:
                        if(col < total_snp):
                                R[row, col] = abs(float(l))/2;
                        col = col + 1;
		col = total_snp;
                for l in line:
                        if(col < 2 * total_snp) :
                                R[row, col] = l;
                        col = col + 1;
                row = row + 1;
	
	[row,col] = R.shape;
	print row;
	print col;
	for command in argv:
		print command;

	Rinv = npy.linalg.inv(R);
	Rdet = npy.linalg.det(R);

	# Generate the casual SNP
        casual_snp = random.sample(xrange(row), casual_implant);
        true_casual = [0] * row;
        for x in casual_snp:
                true_casual[x] = 1;

	Z = npy.random.multivariate_normal((NCP * npy.matrix(true_casual) * R).tolist()[0], R);	

	pred_casual = find_optimal_merge_memory(range(row), npy.matrix(Z), R, Rinv, Rdet, power, alpha, ro, prob_casual).tolist()[0];

	error = pred_casual - npy.matrix(true_casual);
	
	[erow, ecol] = error[error < 0].shape;
	print true_casual;
	print pred_casual;
	print "Total Error", ecol; 
	print "Total SNP picked", sum(pred_casual);

	one_indexs = [i for i,val in enumerate(true_casual) if val==1];
	print one_indexs;

	outputFile = open(argv[0],'w');
	#write the total SNP picked
	outputFile.write(str(sum(pred_casual)));
	outputFile.write("\n");
	#write the set of SNPs picked
	write2File(outputFile, map(int, pred_casual));
	outputFile.write("\n");
	#write the true casual SNPs
	write2File(outputFile, true_casual);
	outputFile.write("\n");
	#write the index of sort
	sort_Z = sorted( [abs(x) for x in Z], reverse=True);
	print sort_Z;
	for i in one_indexs:
		outputFile.write(str(sort_Z.index(abs(Z[i]))+1));
		outputFile.write(' ');
	outputFile.write('\n');	
	#write the condtional restls
	outputFile.write("1 2 3 4 5 6\n");	
	#write the Z score
	write2File(outputFile, Z);
	outputFile.write("\n");
	#write if we have enough power 0 not enough power 1 enogh power
	#outputFile.write('0' if len(Z>=NCP)==0 else '1');	
	power_flag = 0;
	for val in Z:
		if(abs(val) >= NCP) :
			power_flag = 1;
	outputFile.write(str(power_flag));
	outputFile.write("\n");
	
	outputFile.close();

if __name__ == "__main__":
    main(sys.argv[1:])
