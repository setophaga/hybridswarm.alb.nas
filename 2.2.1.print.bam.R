library(splitstackshape)
b=read.csv("~/Library/Containers/com.microsoft.Excel/Data/Desktop/hybridswarm/1.pipeline/sorted.bam.list", header=F)
for(i in b$V1)
{cat(i, " ")}


p=cSplit(b, 'V1', sep=".", type.convert=FALSE)
for(i in p$V1_1)
{cat(i, " 2", '\n')}