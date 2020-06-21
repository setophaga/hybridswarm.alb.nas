library(splitstackshape)

setwd("~/Desktop/hybridswarm/1.pipeline/parentalstrain/")
nas=read.csv("nas00.frq", sep="\t", row.names=NULL)
alb=read.csv("alb03.frq", sep="\t", row.names=NULL)

colnames.change=function(d){colnames(d)=c("CHROM", "POS", "N_ALLELES", "N_CHR", "G1", "G2"); return(d)}
alb=colnames.change(alb)
nas=colnames.change(nas)

#because the G1 and G2 have allele base info and allele freq in the same column, we need to split each of such column into two
text.to.col=function(d)
{d=cSplit(d, 'G1', sep=":", type.convert=FALSE);d=cSplit(d, 'G2', sep=":", type.convert=FALSE)
d=na.omit(d)
d$G1_2=as.numeric(d$G1_2); d$G2_2=as.numeric(d$G2_2)
return(d)}
alb=text.to.col(alb); nas=text.to.col(nas)

#check if nas and alb have the same allele orders 
sum(alb$G1_1==nas$G1_1)==nrow(nas)

#check if each SNP has allele freqs sum to 1
sum(alb$G1_2+alb$G2_2)==nrow(nas)

#find the alleles that has difference in nas$G1_2 and alb$
#fraction of the genome with fixed difference 
length(which(abs(alb$G1_2-nas$G1_2)==1))/nrow(nas)
fixed.snps=which(abs(alb$G1_2-nas$G1_2)==1)

albf=alb[fixed.snps,]; nasf=nas[fixed.snps,]

#we need the following columns:
#CHROM	POS	

#N_CHR.p1	ALLELE1	A1.freq.p1	ALLELE2	A2.freq.p1	
#N_CHR.p2	ALLELE1.p2	A1.freq.p2  ALLELE2.p2	A2.freq.p2
d=data.frame(albf$CHROM, albf$POS, albf$N_CHR, albf$G1_1, albf$G1_2, albf$G2_1, albf$G2_2,nasf$N_CHR, nasf$G1_1, nasf$G1_2, nasf$G2_1, nasf$G2_2)
head(d)
colnames(d)=c("CHROM", "POS", "N_CHR.alb",	"ALLELE1.alb"	,"A1.freq.alb",	"ALLELE2.alb",	"A2.freq.alb",	"N_CHR.nas",	"ALLELE1.nas",	"A1.freq.nas",  "ALLELE2.nas", "A2.freq.nas")
write.csv(d, "alb03.nas00.diff.csv", row.names=F)