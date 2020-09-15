# hybridswarm.alb.nas
This direcotry has manuscript and submission-related information. \
It will be a working progress until submission. 

/Rcodes directory contains the R code involved in analyzing the data, generating figures, and includes a sub-directory with input data. 

Analytical Pipeline: 

**Step1: parental strain/species ancestry reference** \
1.1 strain-specific high depth and coverage references aligned to reference (the nearest outgroup)with **bwa** \
	code: **1.1.alb03.nas00.reference.alignment.sh** \
1.2 genotype ancestry reference with **gatk** \
     1.2.1 genotyping with gatk 
          **1.2.1.gatk.sh**    \
     1.2.2 filter SNPs (genotype quality>20, missing data =0, indel =0, minN =maxN alleles=2, maf=0.05)
          
     vcftools --vcf alb03.nas00.vcf --max-alleles 2 --max-missing 0.8  --minGQ 20  --recode --recode-INFO-all --out          
     alb03.nas00.filtered.vcf
1.3 allelefreq calculation with vcftools > SNPs that are different between parent1 and parent2 
      
      vcftools --vcf alb03.nas00.filtered.vcf.recode.vcf --keep nas00.list --freq --out nas00 
      vcftools --vcf alb03.nas00.filtered.vcf.recode.vcf --keep alb03.list --freq --out alb03 
      
1.4 clean up the frq out put to make:to make **alb03.nas00.diffs.csv** file - **1.4.alb03.nas00.fixed.diff.R** \
or make **alb03.nas00.diffs0.3.csv** \
& haplotyping of neo-sex chromosome, generate **alb03.male.female.csv** with **1.4.alb03.nas00.nonfixed.diff.R **
   
   --make sure that the csv has the **following columns 
     
     CHROM	POS	N_CHR.p1	ALLELE1	A1.freq.p1	ALLELE2	A2.freq.p1	N_CHR.p2	ALLELE1.p2	A1.freq.p2	     
     ALLELE2.p2	A2.freq.p2
 

**Step2: Ancestry HMM** \
2.1 hybrid sequences align to the same reference (as step 1.1) \
   2.1.1 align with bwa e.g. 
  
	while read prefix;do;
	bwa mem -M ref/kepul03FinalMaskedDrorep.fa more.fastqs/"$prefix".fastq > sam/"$prefix".sam
	samtools view -S -b  sam/"$prefix".sam > bam/"$prefix".bam 
	done < alb.nas.17add.list
   2.1.2 sort the bam files 
  
       while read prefix
          do 
       	samtools sort $prefix.bam -o $prefix.sorted.bam;
	samtools index -b $prefix.sorted.bam;
	samtools idxstats $prefix.sorted.bam > $prefix.idxstats;
         done < prefix.list 
   #this outout .idxstats file that contains read counts and length of each muller element
   2.1.3 sex each hybrid individual \
   	**2.1.3.sexing.bam.indv.R** \
2.2 run Ancestry_HMM on the bam files and the csv file from **Step1** \
   2.2.1 make *mpileup.txt* file 
    
     samtools mpileup -q20 ind1.bam ind2.bam [...] indn.bam  > mpileup.txt 
   or use a bamlist containing bam file paths in each line \
   
   	samtools mpileup -q20 --bam-list plate1.bamlistwothersp > mpileup.newlib.p1.txt

   #to input hundres of bam file here, use this code **2.2.1.print.bam.R** and copy&paste the output
    #Here, each bam would correspond to a single sample that you want to perform LAI on. 
    #need ArgParse package 
    Getopt::ArgParse #package which you can install by 
    
     sudo cpan Getopt::ArgParse 
   2.2.2 make *ahmm.in* input file \
     use the .frq file from 1.4 **alb03.nas00.diffs.csv** 
     
   	 perl identify_AIMs.pl --ANGSD alb03.nas00.diffs.csv --mpileup mpileup.txt --output ahmm.input
   #run this line for non-fixed ancestry reference, generate input 
   	
	 perl identify_AIMs.pl --ANGSD alb03.nas00.all.csv --num_pop 2 --mpileup mpileup.txt --output ahmm.input.all 
   #run this line for neo-X and neo-Y-specific sites albm, generate input 
   	 
	 perl identify_AIMs.pl --ANGSD alb03.male.female.2p.csv --mpileup mpileup.txt --output ahmm.input.albFM.2p 
	 #this line below is to haplotype in parallele
	 perl identify_AIMs.pl --ANGSD alb03.male.female.nas00.csv --num_pop 3 --mpileup mpileup.txt --output ahmm.input.albFM.nas 
   	
   2.2.3 run Ancestry HMM
    
    ancestry_hmm -i ahmm.input -s sample.list -a 2 0.5 0.5 -p 0 -3 0.5 -p 1 -3 0.5 -r 0.000005
   #run this line for non-fixed ancestry reference 
    
    ancestry_hmm -i ahmm.input0.3 -s sample.list -a 2 0.5 0.5 -p 0 -3 0.5 -p 1 -3 0.5 -r 0.000005  
    cd posterior.all.alb03.nas00
    ancestry_hmm -i ahmm.input.sp.all -s sample.list -a 2 0.5 0.5 -p 0 -3 0.5 -p 1 -3 0.5 -r 0.000005 
    
   ##run this line for alb 0=neoX (prop=0.3325), 1=neoY(prop=0.3325), and 2=nas (prop=0.0.335) ancestry
   
    cd posterior.albFM.nas
    ancestry_hmm -i ahmm.input.albneoxy.nas -s sample.list -a 3 0.3325 0.3325 0.335 -p 0 -3 0.3325 -p 1 -3 0.3325 -p 2 -3 0.335 -r 0.000005 
    
   #note: to make sample.list use **print.bam.R** and copy&paste the output
    
    There are also a few optional arguments: 
    1. -m is the minimum distance in bp between two AIMs 
    2. -r is the recombination rate in morgans/bp 
    3. --min_p1 is the minimum number of chromosomes from population 1 
    4. --min_p2 is the minimum number of chromosomes from population 2 
    5. --freq_diff is the minimum allele frequency difference between populations for including a site 

**Step3: processing AHMM output, filtration** \
	3.1 hierachical haplotyping within muller CD (5 haplotypes: nas,nas; nas,neoY; nas,neoX; neoX,neoY; neoX,neoX \
	3.2 find mode of muller CD haplotype in each individual, find recombination events, calculate het and ancestry for each muller element \
	code: **3.Genotype.blocks.neosex.species.cluster.R** \
**Step4: plotting and stats** \
	code: **4.1.introgression.recomb.R** \
	4.2. inter-chromosomal LD and network \
	steps: remove false individuals (excessive recombination)-> kmeans cluster the haplotype blocks -> LD among clusters 
	
	
	
	
	

