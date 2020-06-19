#!/bin/bash
sam="/scratch/silu/abo.nas/alb03.nas00.ref/sam"
bam="/scratch/silu/abo.nas/alb03.nas00.ref/bam"
log="/scratch/silu/abo.nas/alb03.nas00.ref/logs"
gvcf="/scratch/silu/abo.nas/alb03.nas00.ref/gvcf"
ref="/scratch/silu/abo.nas/abo.nas.hybrids/ref/kepul03FinalMaskedDrorep.fa"

#make dict for reference 
/scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java -jar /scratch/silu/abo.nas/picard-2.21.6/picard.jar CreateSequenceDictionary REFERENCE=/scratch/silu/abo.nas/abo.nas.hybrids/ref/kepul03FinalMaskedDrorep.fa OUTPUT=/scratch/silu/abo.nas/abo.nas.hybrids/ref/kepul03FinalMaskedDrorep.dict
#make index for reference 
samtools faidx $ref

##add readgroups to bam 
while read prefix do
/scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java -jar /scratch/silu/abo.nas/picard-2.21.6/picard.jar AddOrReplaceReadGroups \
 I=$prefix.bam O=$prefix.RG.bam\
 done < prefix.list

# tell it where the executables are
gatk='/scratch/datmai/jbods_data00_datmai/software/
GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar'

while read prefix do

# load and assign locations for executables
module load java

#Call GATK HaplotypeCaller
/scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java -jar $gatk \
-nct 8 \ 
-l INFO \ 
-R $ref 
\ -log $log/$prefix.HaplotypeCaller.log \ 
-T HaplotypeCaller \ 
-I $bam/$prefix.RG.bam \ 
--emitRefConfidence GVCF \
--max_alternate_alleles 2 \ 
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
 -o $gvcf/$prefix.GATK.gvcf.vcf \
done < prefix.list

#combine gvcf into vcf
module load java

#combina gvcf into a single file
/scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java/-jar $gatk \
-T GenotypeGVCFs \
-R $ref \
-V prefix.gatk.list \
-o alb03.nas00.vcf \
-log $log/gvcf.intoVCF.log





