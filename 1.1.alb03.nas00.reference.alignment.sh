#parental lines alb03 and nas00
/jbods/data01/READS/DBMN21_22/DBMN21-D_S4_L007_R*_001.fastq.gz	albomicans	15112-1751.03	female
/jbods/data01/READS/DBMN30/DBMN30-16_S49_L008_R*_001.fastq.gz	albomicans	15112-1751.03	Male
/jbods/data01/READS/DBMN30/DBMN30-19_S51_L008_R*_001.fastq.gz	albomicans	15112-1751.03	Male
/jbods/data01/READS/DBCC035/DBCC035C4_S68_L008_R1_001.fastq.gz	D. nasuta	00	M
/jbods/data01/READS/DBMN21_22/DBMN21-B_S2_L007_R*_001.fastq.gz	nasuta



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













