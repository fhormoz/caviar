import scipy.stats as st
import numpy as npy, numpy.linalg
from math import *
import random 
import csv
import sys
import optparse
from math import factorial
from ctypes import *

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


def p_value(stat) :
	return(2 * st.norm.pdf(-abs(stat), loc=0, scale=1));

def convertZ2Pvalue(Z) :
	result=list();
	for data in Z:
		result.append(p_value(data));
	return(result);

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
        return( ( 1 / npy.sqrt(pow(2*pi, col) * abs(Rdet)) )   * exp ( -0.5 * (Z-mu) *  Rinv * (Z-mu).T  ) )

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

def find_optimal_merge_memory(potential_casual, Z, R, Rinv, Rdet, power, alpha,ro, prob_casual, est_beta) :
      	index = 0;
	num1 = 0; 
	total_sum = 0.0;
        total_casual_cut = 5;
        N = len(potential_casual);
	NCP = est_beta;
        print 'NCP=', NCP;
	_binary = CDLL('/home/fhormoz/code/Posterior/libbin.so')
	data = '0' * N;
        c_size = c_int(N);	
	c_NCP  = c_double(NCP);
	c_Z = Z.ctypes.data_as(c_void_p);
	c_configure = c_char_p(data);
	vector_R =  npy.reshape(R, N*N);
	c_total_sum = c_double(0.0);
	c_ro = c_double(ro);
	c_R = vector_R.ctypes.data_as(c_void_p);	
	print vector_R;
	print c_Z;
	print c_NCP.value;	
	c_total_sum = _binary.findOptimalSetGreedy(c_Z, c_R,  c_size, c_NCP, c_configure, c_ro);
	curr_binary_value = npy.matrix(list(c_configure.value)[0:N], dtype=npy.uint8);
	print c_configure.value;
	print curr_binary_value;
	print c_total_sum;
	print "WE USE HURISTIC FOR HERE"
	return(curr_binary_value);

#computes the conditional causal SNP using iterative procedure
# if $i$ SNP is picked causal
# the z-score for the rest is
# Z_j^{new} = (Z_j^{old} - r_{ij} Z_j)/(sqrt(1-r_{ij}*r_{ij}))
def conditional_iterative_method(Z, R, causalCount=0) :
	P = pvalues(Z);
	Z_copy = Z;
	print Z.shape;
	print P;
	snpcount = Z.shape[0];
	outList = [];
	if(causalCount == 0) :
		while npy.any(P < 0.05/snpcount) :
			index = numpy.argmax(abs(Z_copy));
			denominator = npy.subtract(npy.ones(shape=(1,snpcount)), npy.multiply(R[index,:], R[index,:]));
			denominator[denominator==0] = 1;
			Z_copy = npy.divide( Z_copy - npy.multiply(R[index,:], Z_copy), npy.sqrt(denominator) );
			Z_copy = npy.nan_to_num(Z_copy);
			outList.append(index);
			P = pvalues(Z_copy);
			print Z_copy;
			print outList; 	
		return(outList);
	else :
		while npy.any(len(outList) != causalCount) :
			index = numpy.argmax(abs(Z_copy));
			denominator = npy.subtract(npy.ones(shape=(1,snpcount)), npy.multiply(R[index,:], R[index,:]));
			denominator[denominator==0] = 1;
			Z_copy = npy.divide( Z_copy - npy.multiply(R[index,:], Z_copy), npy.sqrt(denominator) );
			Z_copy = npy.nan_to_num(Z_copy);
			outList.append(index);
			P = pvalues(Z_copy);
			print "Z cond=", Z_copy;
			print "Causal List CM=", outList; 	
		return(outList);


#We use the approximation method same as et.al. Nature genetcs
# In R:
# p                               # Vector of p-values
# w      <- exp(.5*qchisq(1-p,1)) # Calculate AIC/BIC weights
# w      <- w/sum(w)              # Calculate approximate posterior probability
def approx_posterier_blog(Z,R, causalCount=0) :
	snps=[];
	CAP = 80;
	P = pvalues(Z);
	snpcount = Z.shape[0];
	w = st.chi2.ppf(npy.subtract(npy.ones(shape=(1,snpcount)), P),1);
	print w;
	w[npy.isinf(w)] = CAP
	tmp = [0.5]*snpcount;
	w = npy.exp(npy.multiply(npy.matrix(tmp), w));
	w = npy.divide(w, npy.sum(w));
	tmp = [-1]*snpcount;
	negw = npy.multiply(npy.matrix(tmp), w);
	print "w2=",w; 	
	index=npy.argsort(negw);
	sum = 0;
	if(causalCount==0) :
		for value in  numpy.array(index)[0]:
			sum = sum + w[0,value];
			snps.append(value);
			if(sum > 0.95) :
				print snps;
				return snps;
	else:
		for value in  numpy.array(index)[0]:
                        sum = sum + w[0,value];
                        snps.append(value);
                        if(len(snps) == causalCount) :
                                print snps;
                                return snps;
		



def main(parser):
	row = 1;
	col = 1;
	power = 0.5;
	alpha = 1e-08;
	prob_casual = 1/35;	
	casual_implant = 1;
	total_snp = 0;	
	NCP = calculate_cutoff_power(100, power, alpha) * sqrt(100);
	print NCP;
	#Parse the input
	(options, args) = parser.parse_args();
	outputFile = options.outfile;
	casual_implant =  int(options.causal_count);
	ldFile =  options.ldfile;
	est_beta = NCP * float(options.beta_fraction);
	ext_flag = int(options.ext);	#0 means use original method and 1 means use our extention model
	rho = float(options.rho)

	inputFile = open(ldFile, 'rb');
	for data in inputFile:
		total_snp = total_snp + 1;
	R = npy.matrix(npy.zeros(shape = (total_snp, total_snp)));	
	inputFile.close();	
	row = 0;
	col = 0;
	print total_snp;
	inputFile = open(ldFile, 'rb');
	for data in inputFile:
        	line = data.split(' ');
		col = 0;
		for l in line:
			if(col < total_snp):
				R[row, col] = float(l) ;
			col = col + 1;
		row = row + 1;
	R = R[1:15,1:15];	
	[row,col] = R.shape;
	print row;
	print col;
	Rinv = R;
	Rdet = npy.linalg.det(R);	
	#TODO
	# Generate the casual SNP
        casual_snp = random.sample(xrange(row), casual_implant);
	true_casual = [0] * row;
	for x in casual_snp:
                true_casual[x] = 1;
	Z = npy.random.multivariate_normal((NCP * npy.matrix(true_casual) * R).tolist()[0], R);	
	print true_casual;
	pred_casual = [];
	pred_casual = find_optimal_merge_memory(range(row), npy.matrix(Z), R, Rinv, Rdet, power, alpha, rho, prob_casual, est_beta).tolist()[0];
	pred_casual = [(-1 if num==2 else num) for num in pred_casual]
	error = 0;
	for x in range(0,row-1):
		if((true_casual[x] != 0) & (true_casual[x] != abs(pred_casual[x]))):
			error = error + 1;
	print true_casual;
	#print pred_casual;
	print "Total Error", error; 
	print "Total SNP picked", str(sum([abs(x) for x in pred_casual]));
	one_indexs = [i for i,val in enumerate(true_casual) if val!=0];
	print one_indexs;

	outputFile = open(outputFile,'w');
	#write the total SNP picked
	outputFile.write(str(sum([abs(x) for x in pred_casual])));
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
	#write the condtional results
	if(ext_flag==1) :
		rank_cond = conditional_iterative_method(Z,R, sum(pred_casual));	
	else :
		rank_cond = conditional_iterative_method(Z,R);
	cond_causal = [0]*row;  
	for value in rank_cond:
		cond_causal[value] = 1;
	write2File(outputFile, cond_causal);	
	outputFile.write("\n");
	print cond_causal;
	#wrote the posterier method
	if(ext_flag == 1):
		rank_pos = approx_posterier_blog(Z,R, sum(pred_casual));
	else:
		rank_pos = approx_posterier_blog(Z,R);
	pos_causal = [0] *row;
	for value in rank_pos:
		pos_causal[value] = 1;
	write2File(outputFile, pos_causal);
        outputFile.write("\n");
	#write the Z score
	write2File(outputFile, Z);
	outputFile.write("\n");
	write2File(outputFile, convertZ2Pvalue(Z));
	outputFile.write("\n");
	#write if we have enough power 1, not enough power 0
	#outputFile.write('1' if len(Z>=NCP)==0 else '1');	
	power_flag = 0;
	if(any(abs(Z) >= NCP)) :
		power_flag = 1;
	outputFile.write(str(power_flag));
	outputFile.write("\n");	
	outputFile.close();

if __name__ == "__main__":
	parser = optparse.OptionParser("usage: %prog [options] ")
	parser.add_option("-o", "--out", dest="outfile",
		default="out", type="string",
		help="specify the output file")
	parser.add_option("-l", "--ld_file", dest="ldfile", default="",
		type="string", help="the ld input file")
	parser.add_option("-c", "--causal_count", dest="causal_count", default=1,
                type="int", help="the number of causal SNP implanted in the region")
	parser.add_option("-b", "--beta_est", dest="beta_fraction", default=1,
                type="float", help="the fraction of true beta use to set our estimate beta")
	parser.add_option("-p", "--rho-prob", dest="rho", default=0.95,
                type="float", help="set $\rho$ probability");
	parser.add_option("-e", "--extenstion-compare", dest="ext", default=0,
                type="float", help="if 1 we will run the extention of the conditional and posterier with our number of causal")
	
	
	main(parser);
