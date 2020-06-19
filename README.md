# hybridswarm.alb.nas\
This direcotry has manuscript and submission-related information. \
It will be a working progress until submission. \

/Rcodes directory contains the R code involved in analyzing the data, generating figures, and includes a sub-directory with input data. \

Analytical Pipeline: \
Step1: parental strain/species ancestry reference \
1.1 strain-specific high depth and coverage references aligned to reference (the nearest outgroup)with **bwa** \
1.2 genotype ancestry reference with [gatk] \
1.3 allelefreq calculation with vcftools > SNPs that are different between parent1 and parent2 \
Code 1.1-1.2: alb03.nas00.gatk.sh \
Code 1.3: \
 -vcftools.allelefreq.sh \
 -alb03.nas00.fixed.diff.R \

Step2: Ancestry HMM \
2.1 hybrid sequences align to the same reference (as step 1.1) \
2.2 run Ancestry_HMM on the bam files and the csv file from [Step1] \
#need ArgParse package \
Getopt::ArgParse #package which you can install by \
\
$ sudo cpan Getopt::ArgParse \
perl identify_AIMs.pl --ANGSD p1.p2.fixed.diff.csv --mpileup mpileup.subset.txt --output ahmm.in \
--ANGSD <input CSV file in the same format as before> \
--mpileup <input mpileup file> \
--output <output file for input to ahmm> \
\
There are also a few optional arguments: \
1. -m is the minimum distance in bp between two AIMs \
2. -r is the recombination rate in morgans/bp \
3. --min_p1 is the minimum number of chromosomes from population 1 \
4. --min_p2 is the minimum number of chromosomes from population 2 \
5. --freq_diff is the minimum allele frequency difference between populations for including a site \

\
Step3: processing AHMM output, filtration \
Step4: ancestry genotyp > ancestry cluster \
Step5: BGC \
