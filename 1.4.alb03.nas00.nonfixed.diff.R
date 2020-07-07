library(splitstackshape)

setwd("~/Desktop/hybridswarm/1.pipeline/parentalstrain/")
nas=read.csv("nas00.frq", sep="\t", row.names=NULL)
alb=read.csv("alb03.frq", sep="\t", row.names=NULL)
fstd=read.csv("nas00.alb03.weir.fst",  sep="\t", row.names=NULL)
dim(nas); dim(alb); dim(fstd)

colnames.change=function(d){colnames(d)=c("CHROM", "POS", "N_ALLELES", "N_CHR", "G1", "G2"); d=d[order(d$CHROM, d$POS), ]; return(d)}
alb=colnames.change(alb)
nas=colnames.change(nas)

#because the G1 and G2 have allele base info and allele freq in the same column, we need to split each of such column into two
text.to.col=function(d)
{d=cSplit(d, 'G1', sep=":", type.convert=FALSE);d=cSplit(d, 'G2', sep=":", type.convert=FALSE)
d=na.omit(d)
d$G1_2=as.numeric(d$G1_2); d$G2_2=as.numeric(d$G2_2)
return(d)}
alb=text.to.col(alb); nas=text.to.col(nas)
sp.rows=intersect(which(abs(alb$G1_2-nas$G1_2)>0.2), which(abs(alb$G2_2-nas$G2_2)>0.2))
albd=alb[sp.rows,]; nasd=nas[sp.rows,]
dim(albfd.mullercd=albfd[which(albfd$CHROM=="Muller_DC"),]) #find out how many informative sites are on muller CD
#we need the following columns:
#CHROM	POS	
table(albf$CHROM)
#N_CHR.p1	ALLELE1	A1.freq.p1	ALLELE2	A2.freq.p1	
#N_CHR.p2	ALLELE1.p2	A1.freq.p2  ALLELE2.p2	A2.freq.p2
d=data.frame(albd$CHROM, albd$POS, albd$N_CHR, albd$G1_1, albd$G1_2, albd$G2_1, albd$G2_2,nasd$N_CHR, nasd$G1_1, nasd$G1_2, nasd$G2_1, nasd$G2_2)
head(d)
colnames(d)=c("CHROM", "POS", "N_CHR.alb",	"ALLELE1.alb"	,"A1.freq.alb",	"ALLELE2.alb",	"A2.freq.alb",	"N_CHR.nas",	"ALLELE1.nas",	"A1.freq.nas",  "ALLELE2.nas", "A2.freq.nas")
table(d$CHROM )
#write.csv(d, "alb03.nas00.diffs0.2.csv", row.names=F)
ff=read.csv("alb03.nas00.diffs0.2.csv")
table(ff$CHROM )


#################Add female and males
#male vs females
albm=read.csv("alb.m.DBMN30-16.19.frq", sep="\t", row.names=NULL)
albf=read.csv("alb.f.DBMN21-D_S4_L007.frq", sep="\t", row.names=NULL)
albm=colnames.change(albm);albf=colnames.change(albf)
albm=text.to.col(albm);albf=text.to.col(albf);

albm.mulCD=albm[intersect(sp.rows, which(albm$CHROM=="Muller_DC")),];albf.mulCD=albf[intersect(sp.rows, which(albf$CHROM=="Muller_DC")),]; 
length(intersect(which(albm.mulCD$G1_2==0.5), which(abs(albf.mulCD$G1_2-albf.mulCD$G2_2)==1))) #species-informative && neosex chr informative
length(which(albm.mulCD$G1_2==0.5))
length(which(abs(albf.mulCD$G1_2-albf.mulCD$G2_2)==1))
rows=intersect(which(albm.mulCD$G1_2==0.5),which(abs(albf.mulCD$G1_2-albf.mulCD$G2_2)==1))

alb.mullercd=data.frame(albf.mulCD[rows,], albm.mulCD[rows,])
alb.mullercdd=alb.mullercd[,c(1:2, 4:8,12:16)]
colnames(alb.mullercdd)=c("CHROM", "POS", "N_CHR.f",	"ALLELE1.f"	,"A1.freq.f",	"ALLELE2.f",	"A2.freq.f",	"N_CHR.m",	"ALLELE1.m",	"A1.freq.m",  "ALLELE2.m", "A2.freq.m")
head(alb.mullercdd)

#write.csv(alb.mullercdd, "alb03.mullerCD.male.female.csv", row.names=F)
#Check about Muller B
albm.mulB=albm[intersect(sp.rows, which(albm$CHROM=="Muller_B")),];albf.mulB=albf[intersect(sp.rows, which(albf$CHROM=="Muller_B")),]; 
length(intersect(which(albm.mulB$G1_2==0.5), which(abs(albf.mulB$G1_2-albf.mulB$G2_2)==1))) #species-informative && neosex chr informative
length(which(albm.mulB$G1_2==0.5))
length(which(abs(albf.mulB$G1_2-albf.mulB$G2_2)==1))
rows=intersect(which(albm.mulB$G1_2==0.5),which(abs(albf.mulB$G1_2-albf.mulB$G2_2)==1))

#check about E
albm.mulE=albm[intersect(sp.rows, which(albm$CHROM=="Muller_E")),];albf.mulE=albf[intersect(sp.rows, which(albf$CHROM=="Muller_E")),]; 
length(intersect(which(albm.mulE$G1_2==0.5), which(abs(albf.mulE$G1_2-albf.mulE$G2_2)==1))) #species-informative && neosex chr informative
length(which(albm.mulE$G1_2==0.5))
length(which(abs(albf.mulE$G1_2-albf.mulE$G2_2)==1))
dim(albm.mulE)
rows=intersect(which(albm.mulE$G1_2==0.5),which(aEs(albf.mulE$G1_2-albf.mulE$G2_2)==1))
