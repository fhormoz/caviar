# CAVIAR (CAusal Variants Identication in Associated Regions):
======

## a statistical framework that quantifies the probability of each variant to be causal while allowing with arbitrary number of causal variants.

original publications can be found [here](http://www.cell.com/ajhg/abstract/S0002-9297(16)30439-6) for eCAVIAR and [here](http://www.genetics.org/content/198/2/497) for CAVIAR. 

Developed at [the ZarLab](http://zarlab.cs.ucla.edu/tag/ecaviar/) at UCLA 

More information on CAVIAR and eCAVIAR can be found on the [CAVIAR website](http://genetics.cs.ucla.edu/caviar/)

### information on running CAVIAR and eCAVIAR:

CAVIAR usage 

``` bash 
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










CAVIAR is offered under the GNU Affero GPL (https://www.gnu.org/licenses/why-affero-gpl.html).



