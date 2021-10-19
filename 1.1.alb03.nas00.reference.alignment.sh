#Alignment
ref="/scratch/silu/abo.nas/abo.nas.hybrids/ref/kepul03FinalMaskedDrorep.fa"

#align trimmed reads onto the reference for each parental individual
#eg., DBMN21-B_S2_L007
bwa mem $ref DBMN21-B_S2_L007_R1_001.fastq DBMN21-B_S2_L007_R2_001.fastq >  bam/DBMN21-B_S2_L007.bam



#If the header of the SAM file is improperly formatted from BWA (i.e. the first line does not start with @HD), you can use the reference dictionary created with Picard in Step 1 to fix it:
#java -jar picard.jar ReplaceSamHeader I=reads-mapped.sam HEADER=reference.dict O=reads-mapped-hd.sam


###STEPS to sex the individuals
#sort bam file
samtools sort DBCC035C4_S68_L008.bam -o DBCC035C4_S68_L008.sorted.bam
samtools index -b DBCC035C4_S68_L008.sorted.bam
samtools idxstats DBCC035C4_S68_L008.sorted.bam > DBCC035C4_S68_L008.idxstats

#add read group for GATK input
samtools sort bam.albomref/"$prefix".bam -o bam.albomref/"$prefix".sorted.bam
/scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java -jar /scratch/silu/x01.albom.nas/tools/picard-2.21.6/picard.jar AddOrReplaceReadGroups I=bam.albomref/$prefix.sorted.bam O=bam.albomref/$prefix.RG.bam SORT_ORDER=coordinate RGID=alb03 RGPU=CBKVCANXX RGLB=Bachtrog RGPL=ILLUMINA RGSM=$prefix CREATE_INDEX=TRUE












