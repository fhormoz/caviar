library(MASS)
library(mvtnorm)
source('/home/fhormoz/code/Posterior/powerSim.R')


## Compute the p-value

p_value <- function (stat) 
{
	return(2 * pnorm(-abs(stat), mean=0, sd=1, lower.tail = TRUE));
}

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
	return(rev(out));
}


all_next_permutation <- function (data, size)
{
	i <- 1;
	
	while(data[i] == 0 && i <= size)
		i <- i + 1;

	last_one_index <- i;
	perm_num <- last_one_index - 1;
	
	out <- matrix(data, nrow=perm_num, ncol=size, byrow=T); 
	
	i <- 1;
	while(i <= perm_num) {
		out[i,last_one_index-i] <- 1;
		i <- i+1;
	}
	return(out)	
}


# Generate All the possible subset of a set 
# The ouput is a vector of 1 and 0, where 1
# mean the element exist in the subset and 0
# otherwise
subsets_order <- function(set, N)
{
        size <- length(set);

        #intialize the zero for the container of subsets
        zeros <- rep(0, N * 2^size);
        
	res_subsets <- matrix(zeros, 2^size, N);
	indexs      <- rep(1, 2^size);

        index <- 1;

	index1 <- 1;
	index2 <- 1;
	index3 <- 1;	

        while(index <= size)
        {
		tmpIndex <- index1;
              	while(tmpIndex <= index2) {
			res <- all_next_permutation(res_subsets[tmpIndex, ], size);
			tmpIndex2 <- 1;
			while(tmpIndex2 <= nrow(res)) {
				res_subsets[tmpIndex2+index3,] = res[tmpIndex2,];
				indexs[tmpIndex2+index3] = tmpIndex;
				#cat(tmpIndex2+index3,'\t',res_subsets[tmpIndex2+index3,], '\n');
				tmpIndex2 <- tmpIndex2 + 1;
			}
			index3 <- index3 + nrow(res);
			tmpIndex <- tmpIndex + 1;
		}
		index1 <- index2+1;
		index2 <- index3;  
		index <- index + 1;
        }
	res_subsets[2^size,] = rep(1,size)
        return(list(set=res_subsets, index=indexs));
}


fast_find_optimal_sets <- function(potential_casual, Z, R, power, alpha, ro)
{
	index1 <- 1;
	index2 <- 1;

	N <- length(Z);
	NCP <- calculate_cutoff_power(100, power, alpha) * sqrt(100);

	set_index <- subsets_order(potential_casual, N);
	
	sets <- set_index$set;
	orig_index <- set_index$index;

	stat_values <- rep(0, nrow(sets));	
	
	while(index1 <= nrow(sets)) {
		stat_values[index1] <- likehood_est(sets[index1,], Z, R, NCP);
		index1 <- index1 + 1;
	}
	
	total_sum <- sum (stat_values);
	tmpvalues <- rep(0, nrow(sets));
	index1 <- 1;

	cat(total_sum, '\n');

	while(index1 <= nrow(sets)) {
		index2 <- orig_index[index1];
		while(index2 <= index1) {
			if(length(sets[sets[index1,]-sets[index2,]<0]) == 0) {
				tmpvalues[index1] <- tmpvalues[index1] + stat_values[index2];
				if(tmpvalues[index1]/total_sum > ro) {
					return(list(sets=sets, index=index1, value=stat_values));
				}
			}
			index2 <- index2 + 1;
		}
		#cat(sets[index1,], '\n');
		index1 <- index1 + 1;
	}

	return(list(sets=sets, index=index1, value=stat_values));

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
				res_subsets[index+1,set[index2]] <- 1;
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

binary2Int<- function(configure) 
{
	sum(2^(which(rev(unlist(strsplit(as.character(configure), "")) == 1))-1));
}

all_possible <- function(configure) 
{
	N = length(configure);

	num_1 = sum(configure);
	num_0 = N - num_1;

	tmp <- matrix(configure, nrow = num_1+1, ncol = N, byrow=T);

	output <- rep(1, num_1);

	index1 <- 1;
	index2 <- 1;
	while(index1 <= N) {
		if(configure[index1]!=0) {
			tmp[index2,index1] = 0;
			output[index2] <- binary2Int(tmp[index2,])+1;
			index2 <- index2 + 1;	
		}	
		index1 <- index1 + 1;
	}
	return(output);
}

myxor <- function(X, Y) 
{
	N = length(X);
	index <- 1;
	out <- rep(0,N);
	while(index <= N) {
		if(X[index] == 1 && Y[index]==1) {
			out[index] = 1;
		}
		index <- index + 1;
	}	
	return(out);
}

remove_redundent_sets <- function(configure)
{
	N = length(configure);

        num_1 = sum(configure);
        num_0 = N - num_1;

        tmp <- matrix(configure, nrow = num_1, ncol = N, byrow=T);

        output <- c();

        index1 <- 1;
        index2 <- 1;
        
	while(index1 <= N) {
                if(configure[index1]!=0) {
                        tmp[index2,index1] = 0;
                        index2 <- index2 + 1;
                } 
                index1 <- index1 + 1;
        }
	
	index1 <- 1;
	index2 <- 1;

	while(index1 <= num_1) {
		index2 <- index1+1;
		while(index2 <= num_1) {
			set <- myxor(tmp[index2,] , tmp[index1,]);
			output <- cbind(output, binary2Int(set)+1);
			index2 <- index2 + 1;	
		}
		index1 <- index1 + 1;
	}	

        return(output);

}

one2set <- function(configure, N)
{
	out <- rep(0, N);
	index1 <- 1;
	while(index1 <= length(configure)) {
		out[configure[index1]] = 1;	
		index1 <- index1 + 1;
	}
	return(out);
}

all_possible_perm <- function(configure)
{
	set <- c();
	index1 <- 1;
	N <- length(configure);

	while(index1 <= N) {
		if(configure[index1]==1)
			set <- cbind(set,index1);
		index1 <- index1 + 1;
	}
	sets <- subsets(set, N);
	
	output <- rep(0, 2^sum(configure))
	index1 <- 1;
	total_size <- nrow(sets);
	while(index1 <= total_size) {
		output[index1] = binary2Int(sets[index1,])+1;
		index1 <- index1 + 1;
	}	

	return (output)
}

find_optimal_merge <- function(potential_casual, Z, R, power, alpha,ro) 
{
	sum <- 0;
	index1 <- 1;
	
	N = length(potential_casual);
	
	NCP <- calculate_cutoff_power(100, power, alpha) * sqrt(100);
	
	sets <- subsets(potential_casual, N);	
	results <- subsets_order(potential_casual, N);
	order_sets <- results$set;

	values <- rep(0, nrow(sets));
	tmpvalues <- rep(0, nrow(sets));
	

        while(index1 <= nrow(sets))
        {
                sum <- sum + likehood_est(sets[index1,], Z, R, NCP);
                values[index1] <-  likehood_est(sets[index1,], Z, R, NCP);
		index1 <- index1 + 1;
        }

	cat('sum=', sum, '\n');
		
#	tmpvalues <- values;

	index1 <- 2;
	index2 <- 1;
	while(index1 <= nrow(order_sets)) {
		configure <- order_sets[index1,];
		#cat(index1,'\t',configure, '\t');
		all_configure <- all_possible_perm(configure);
		index2 <- 1;
		while(index2 <= length(all_configure)) {
			tmpvalues[binary2Int(configure)+1] <- tmpvalues[binary2Int(configure)+1] + values[all_configure[index2]];
			index2 <- index2 + 1;
		}
		if(tmpvalues[binary2Int(configure)+1]/sum >= ro) {
			return(list(orderset=order_sets, index=index1, sets=sets, value=values));
		}	
		index1 <- index1 + 1;
	}
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


conditional_like_lihood <- function(Z, R, NCP)
{
	MIN_VALUE <- 0;
	result <- rep(0, length(Z));

	tmpZ <- Z;	
	index <- 1;	
	cat(tmpZ, '\n')
	while(index < length(Z)) {
		indexSortValue <- order(-abs(tmpZ));
		result[index] <- indexSortValue[1];
		tmpZ <- tmpZ - tmpZ[indexSortValue[1]] * R[indexSortValue[1],];
		cat(tmpZ, '\n');
		tmpZ[indexSortValue[1]] <- MIN_VALUE;
		index <- index + 1;
	}	
	return(result);
}


# Z is the z-score from the EMMA
# newR is the correlation between different SNPs
# ro is the cut-off of postier we like to reach 


args <- commandArgs(trailingOnly = TRUE);

row   <- 20;
power <- 0.5;
alpha <- 1e-08;

ro <- as.double(args[3]);
#ro <- 0.95;

set.seed(as.integer(args[2]))

#true_casual <- c(0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0);
#true_casual <- c(0,0,1,0,0,0,0,0,0,0);
true_casual <- rep(0,row);
true_casual[sample(1:row,2)] <- 1;	# make this SNP to be first casual

NCP <- calculate_cutoff_power(100, power, alpha) * sqrt(100);

number_snps_peak <- sqrt(length(scan("/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld")));

tmp <- matrix(scan("/home/fhormoz/code/Posterior/data/peakSNP_100kb/peakSNP_100kb.ld"), number_snps_peak, number_snps_peak);
snpIndexs <- (as.integer(number_snps_peak/2)-as.integer(row/2)):(as.integer(number_snps_peak/2)-as.integer(row/2)+row-1);

newR <- tmp[snpIndexs, snpIndexs];

index1 <- 1;
index2 <- 1;

while(index1 <= row) {
	index2 <- index1;
	while(index2 <= row) {
		newR[index2, index1] <- newR[index1, index2];
		index2 <- index2 + 1;
	}
	cat(newR[index1,],'\n');
	index1 <- index1 + 1;
}



casual <- matrix(true_casual,1,row);
Z <- rmvnorm(1, NCP %*% casual %*% newR, newR);

cat('Z=', Z, '\n');

results     <- find_optimal_merge(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), Z, newR, power, alpha, ro);
#results     <- find_optimal_merge(c(1,2,3,4,5,6,7,8,9,10), Z, newR, power, alpha, ro);
order_set   <- results$orderset;
minIndex    <- results$index;
sets        <- results$sets;
values      <- results$value;
sort_values <- values[order(-values)];

cond_sets   <- conditional_like_lihood(Z, newR, NCP)

output1 <- sum(order_set[minIndex,]);	
output2 <- order_set[minIndex,];
output3 <- which(sort_values==values[binary2Int(true_casual)+1]);

write(output1, file = toString(args[1]));
write(output2, file = toString(args[1]), append=TRUE, ncol=row);
write(true_casual, file = toString(args[1]), append=TRUE, ncol=row);
write(output3, file = toString(args[1]), append=TRUE);
write(cond_sets, file = toString(args[1]), append=TRUE, ncol=row);
