trimmomatic="/opt/trinityrnaseq-v2.9.0/trinity-plugins/Trimmomatic-0.36/trimmomatic-0.36.jar"
while read prefix
do
java -jar $trimmomatic PE -phred33 -threads 1 "$prefix".fq1 "$prefix".fq2 trimmed/"$prefix"_R1_001.fastq trimmed/"$prefix"_R1_001.unpaired.fastq trimmed/"$prefix"_R2_001.fastq trimmed/"$prefix"_R2_001.unpaired.fastq  TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:30 
done < prefix.last
