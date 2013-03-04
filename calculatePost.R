library(MASS)
library(mvtnorm)
source('/home/fhormoz/code/Posterior/powerSim.R')

# Configuration is the set of postion which are casual is one and 
#	the rest is 0
# Z is the vector which contain the Z score for all the SNPs
# R is the correlation matrix between each SNPs
# NCP is the non-centrial 
likehood_est <- function(configure,Z, R, NCP)
{
	#mu is the mean if Multi Normal
	#sigma is variance of the Mult Normal which is R
	configure_size = length(configure);
	mu = rep(0,configure_size);
	mu = (NCP * configure) %*% R ;
	
	return (dmvnorm(Z, mu, R));	
}

# Help function to generate all the possible subset 
# this function will convert the decimal to binary
# format
decimal2binary <- function(num, size)
{
	out <- rep(0,size);
	
	tmp <- num;
	index <- 1;

	while(tmp > 0)
	{
		out[index] = tmp%%2;
		tmp = tmp %/% 2;
		index = index + 1;
	}
	return(out);
}

# Generate All the possible subset of a set 
# The ouput is a vector of 1 and 0, where 1
# mean the element exist in the subset and 0
# otherwise
subsets <- function(set, N)
{
	size <- length(set);

	#intialize the zero for the container of subsets
	zeros <- rep(0, N * 2^size);
	res_subsets <- matrix(zeros, 2^size, N);
	
	index <- 1;
	while(index < 2^size)
	{
		tmp <- decimal2binary(index, size);
		
		index2 <- 1;
		while(index2 <= length(tmp))
		{
			if(tmp[index2] == 1)
			{
				res_subsets[index,set[index2]] <- 1;
			}
			index2 <- index2 + 1;
		}  
		index <- index + 1;
	}
	return(res_subsets);	
}


m_corr <- function(X,Y)
{
	return (mean((X-mean(X))*(Y - mean(Y)))/sqrt(var(X))*sqrt(var(Y)));
}




# N is the total number of SNPs
# Z is the z-score from the EMMA
# R is the SNP correlation

find_sets <- function(potential_casual, Z, R, power, alpha)
{
	N <- length(Z);
	NCP <- calculate_cutoff_power(100, power, alpha) * sqrt(100);
	
	index1 <- 1;
	index2 <- 1;

	sum <- 0;

	sets <- subsets(potential_casual, N);
	values <- rep(0, dim(sets)[1]);

	while(index1 < dim(sets)[1])
	{
		sum <- sum + likehood_est(sets[index1,], Z, R, NCP);		
		values[index1] <-  likehood_est(sets[index1,], Z, R, NCP);
		index1 <- index1 + 1;
	} 
	
	tmpvalues <- rep(0, dim(sets)[1]);

	index1 <- 1;
	index2 <- 1;

	while(index1 < dim(sets)[1]) {
		index2 <- 1;
		while(index2 < dim(sets)[1]) {
			#if(index1 != index2) {
				if(length(sets[sets[index1,]-sets[index2,]<0]) == 0) {
					tmpvalues[index1] <- tmpvalues[index1] + values[index2];
				}
			#}
			index2 <- index2 + 1;
		}
		index1 <- index1 + 1;
	}	
	return (list(sets=sets,value=tmpvalues, sum=sum));
}

# return the optimal sets, the input results contain the return value of find_sets return value
find_optimal_sets <- function(results, ro)
{
	#results <- find_sets(potential_casual, Z, newR, 0.5, 1e-08);
	sets    <- results$sets;
	values  <- results$value;
	sumValue     <- results$sum;
	
	index <- 1;
	minIndex <- -1;
	while(index <  dim(sets)[1]) 
	{
		if(values[index]/sumValue > ro) 	
		{
			if(minIndex == -1) {
				minIndex <- index;
			}
			else if( sum(sets[index,]) <  sum(sets[minIndex,]) ) {
				minIndex <- index;
			}
		}	
		index <- index + 1;
	}
	return(minIndex);
}

# N is the total number of SNPs
# Z is the z-score from the EMMA
#

# R calculatePost outoutFILENAME seed_num ro

args <- commandArgs(trailingOnly = TRUE);

row   <- 10;
power <- 0.5;
alpha <- 1e-08;

ro <- as.double(args[3]);

true_casual <- c(0,0,1,0,0,0,0,0,0,0);
set.seed(as.integer(args[2]))
NCP <- calculate_cutoff_power(100, power, alpha) * sqrt(100);


number_snps_peak <- sqrt(length(scan("/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld")));

tmp <- matrix(scan("/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld"), number_snps_peak, number_snps_peak);
rand <- sample(number_snps_peak, row);
newR <- tmp[rand,rand];

index1 <- 1;
index2 <- 1;

while(index1 < row) {
	index2 <- index1;
	while(index2 < row) {
		newR[index2, index1] <- newR[index1, index2];
		index2 <- index2 + 1;
	}
	index1 <- index1 + 1;
}

#tmp <- rnorm(row*row, mean=0.5, sd = 0.2);
#R <- matrix(tmp, row, row);
#newR <- 0.5*R%*%t(R);	
cat(newR,'\n');

casual <- matrix(true_casual,1,row);
Z <- rmvnorm(1, NCP %*% casual %*% newR, newR);
results <- find_sets(c(1,2,3,4,5,6,7,8,9,10), Z, newR, power, alpha);
minIndex <- find_optimal_sets(results, ro);
sets <- results$sets;
output1 <- sum(sets[minIndex,]);	
output2 <- sets[minIndex,];
write(output1, file = toString(args[1]));
write(output2, file = toString(args[1]), append=TRUE, ncol=10);
