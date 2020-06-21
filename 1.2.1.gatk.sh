#!/bin/bash
sam="/scratch/silu/abo.nas/alb03.nas00.ref/sam"
bam="/scratch/silu/abo.nas/alb03.nas00.ref/bam"
log="/scratch/silu/abo.nas/alb03.nas00.ref/logs"
gvcf="/scratch/silu/abo.nas/alb03.nas00.ref/gvcf"
ref="/scratch/silu/abo.nas/abo.nas.hybrids/ref/kepul03FinalMaskedDrorep.fa"

# tell it where the executables are
gatk='/scratch/datmai/jbods_data00_datmai/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar'



        # load and assign locations for executables
        module load java

        #Call GATK HaplotypeCaller
        /scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java -jar $gatk \
                -nct 8 \
                -l INFO \
                -R $ref \
                -log $log/DBCC035C4_S68_L008.HaplotypeCaller.log \
                -T HaplotypeCaller \
                -I  $bam/DBCC035C4_S68_L008.RG.bam \
                --emitRefConfidence GVCF \
                --max_alternate_alleles 2 \
                -variant_index_type LINEAR \
                -variant_index_parameter 128000 \
                -o $gvcf/DBCC035C4_S68_L008.GATK.gvcf.vcf \

##After running all the individual bam file -> gvcf files, combine gvcf files into a single vcf
#!/bin/bash
sam="/scratch/silu/abo.nas/alb03.nas00.ref/sam"
bam="/scratch/silu/abo.nas/alb03.nas00.ref/bam"
log="/scratch/silu/abo.nas/alb03.nas00.ref/logs"
gvcf="/scratch/silu/abo.nas/alb03.nas00.ref/gvcf"
ref="/scratch/silu/abo.nas/abo.nas.hybrids/ref/kepul03FinalMaskedDrorep.fa"
gatk='/scratch/datmai/jbods_data00_datmai/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar'

module load java

#combina gvcf into a single file
/scratch/datmai/bin/000program_files/jre1.8.0_25/bin/java -jar $gatk \
-T GenotypeGVCFs \
-R $ref \
-V prefix.gatk.list \
-o alb03.nas00.vcf \
-log $log/gvcf.intoVCF.log
