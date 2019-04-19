# CAVIAR (CAusal Variants Identication in Associated Regions):
======

## **a statistical framework that quantifies the probability of each variant to be causal while allowing with arbitrary number of causal variants.** 

original publications can be found [here](http://www.cell.com/ajhg/abstract/S0002-9297(16)30439-6) for eCAVIAR and [here](http://www.genetics.org/content/198/2/497) for CAVIAR. 

Developed at [the ZarLab](http://zarlab.cs.ucla.edu/tag/ecaviar/) at UCLA 

More information on CAVIAR and eCAVIAR can be found on the [CAVIAR website](http://genetics.cs.ucla.edu/caviar/)

### information on running CAVIAR and eCAVIAR:

to install this repository -- 

git clone https://github.com/fhormoz/caviar.git

### CAVIAR usage 

```bash 
./CAVIAR [options] 
Options:
-h, --help            		show this help message and exit 
-o OUTFILE, --out=OUTFILE 	specify the output file
-l LDFILE, --ld_file=LDFILE the ld input file
-z ZFILE, --z_file=ZFILE	the z-score and rsID files
-r RHO, --rho-prob=RHO		set $pho$ probability 
-c causal			set the maximum number of causal SNPs
-f 1				to out the probaility of different number of causal SNP
```


### eCAVIAR usage 
```bash 
Usage: ./eCAVIAR [options]

Options:
-h, --help                      show this help message and exit
-o OUTFILE, --out=OUTFILE       specify the output file
-l LDFILE, --ld_file=LDFILE 	the GWAS ld input file
-l LDFILE, --ld_file=LDFILE 	the eQTL ld input file
-z ZFILE, --z_file=ZFILE        the GWAS z-score and rsID files
-z ZFILE, --z_file=ZFILE        the eQTL z-score and rsID files
-r RHO, --rho-prob=RHO          set $pho$ probability
-c causal                       set the maximum number of causal SNPs
-f 1                            to out the probaility of different number of causal SNP
``` 

### Output

OUTFILE_1_set - causal SNP in GWAS

OUTFILE_2_set - causal SNP in eQTL 

OUTFILE_1_post - Causal posterior probability for each SNP in GWAS

OUTFILE_2_post - Causal posterior probability for each SNP in eQTL

OUTFILE_col - The Colocalization posterior probability (CLPP) for each SNP. 


### Debugging 

CAVIAR is written in C++ and must be compiled before running. If you are encountering errors in running CAVIAR or eCAVIAR try these steps:
1. check if you have the [GNU scientific library](https://www.gnu.org/software/gsl/) installed
		* for macOS this can be done using the homebrew package manager- ```brew install gsl```  
2. Next, in the caviar/CAVIAR C++ repository type ```make clean``` 
3. ```make``` 
4. ```chmod +x eCAVIAR``` may also be helpful 
CAVIAR should be able to run using these parameters 

Other helpful hints - if running eCAVIAR make sure your LD files have the same SNPs, works best for low to medium LD 
other related code developed by UCSF students [here](https://github.com/christacaggiano/GWAS-EQTL) 

CAVIAR is offered under the GNU Affero GPL (https://www.gnu.org/licenses/why-affero-gpl.html).



