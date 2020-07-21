library(splitstackshape)

setwd("~/Desktop/hybridswarm/1.pipeline/parentalstrain/")
nas=read.csv("nas00.frq", sep="\t", row.names=NULL)
alb=read.csv("alb03.frq", sep="\t", row.names=NULL)
#fstd=read.csv("nas00.alb03.weir.fst",  sep="\t", row.names=NULL)
dim(nas); dim(alb);

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
albd=alb; nasd=nas
#choose the informative sites with allele freq cutoff being 0.3
#sp.rows=intersect(which(abs(alb$G1_2-nas$G1_2)>0.3), which(abs(alb$G2_2-nas$G2_2)>0.3))
#sp.pos=paste(as.character(alb$CHROM[sp.rows]),as.character(alb$POS[sp.rows]), sep=".")
#albd=alb[sp.rows,]; nasd=nas[sp.rows,]
#we need the following columns:
#CHROM	POS	
#N_CHR.p1	ALLELE1	A1.freq.p1	ALLELE2	A2.freq.p1	
#N_CHR.p2	ALLELE1.p2	A1.freq.p2  ALLELE2.p2	A2.freq.p2
d=data.frame(albd$CHROM, albd$POS, albd$N_CHR, albd$G1_1, albd$G1_2, albd$G2_1, albd$G2_2,nasd$N_CHR, nasd$G1_1, nasd$G1_2, nasd$G2_1, nasd$G2_2)
head(d)
colnames(d)=c("CHROM", "POS", "N_CHR.alb",	"ALLELE1.alb"	,"A1.freq.alb",	"ALLELE2.alb",	"A2.freq.alb",	"N_CHR.nas",	"ALLELE1.nas",	"A1.freq.nas",  "ALLELE2.nas", "A2.freq.nas")
d=na.omit(d)

#plot(albd$G1_2[which(albd$CHROM=="Muller_DC")[1:1000]], nasd$G1_2[which(nasd$CHROM=="Muller_DC")[1:1000]], pch=16, col=rgb(0,1,0, 0.01))
#write.csv(d, "alb03.nas00.all.csv", row.names=F)
ff=d #ff=read.csv("alb03.nas00.diffs0.3.csv")
ff.dc=ff[which(ff$CHROM=="Muller_DC"),]
sum(ff.dc$POS %in% pos.dc)
sum(ff.dc$POS %in% fm.dc$POS)
sum(ff.dc$POS %in% dat.fm.dc$position)
sp.pos=paste(as.character(ff$CHROM),as.character(ff$POS), sep=".")

#table(ff$CHROM )
ff$alb.allele=2
ff$alb.allele[which((ff$A1.freq.alb-ff$A1.freq.nas)>0.3)]=1 #allele A1 is albomicans allele
plot(ff$A1.freq.alb, ff$alb.allele)
#################Add female and males
#male vs females
albm=read.csv("alb.m.DBMN30-16.19.frq", sep="\t", row.names=NULL)
albf=read.csv("alb.f.DBMN21-D_S4_L007.frq", sep="\t", row.names=NULL)
albm=colnames.change(albm);albf=colnames.change(albf)
albm=text.to.col(albm);albf=text.to.col(albf);

full=data.frame(albf[,-3], albm[,-c(1:3)], nas[,-c(1:3)], alb[,-c(1:3)])
colnames(full)=c("CHROM", "POS", "N_CHR.f",	"ALLELE1.f"	,"A1.freq.f",	"ALLELE2.f",	"A2.freq.f",	"N_CHR.m",	"ALLELE1.m",	"A1.freq.m",  "ALLELE2.m", "A2.freq.m","N_CHR.nas",	"ALLELE1.nas",	"A1.freq.nas",  "ALLELE2.nas", "A2.freq.nas","N_CHR.alb",	"ALLELE1.alb",	"A1.freq.alb",  "ALLELE2.alb", "A2.freq.alb" )
full=na.omit(full)
#check if A1 and A2 freq sum to 1 --NO
sum(full$A1.freq.nas+full$A2.freq.nas)==nrow(full); hist(full$A1.freq.nas+full$A2.freq.nas);full[which((full$A1.freq.nas+full$A2.freq.nas)!=1),]

#find alb =A1 locus, nas=A2
alb.a1=intersect(which((full$A1.freq.alb-full$A1.freq.nas)>0.3), (full$A2.freq.nas-full$A2.freq.alb)>0.3)
#find alb =A2 locus, nas=A1
alb.a2=intersect(which((full$A2.freq.alb-full$A2.freq.nas)>0.3), (full$A1.freq.nas-full$A1.freq.alb)>0.3)

#find neoY and neoX specific loci
hetm=intersect(which(full$A1.freq.m==full$A2.freq.m),which(full$A2.freq.m==0.5))
f.a1=intersect(which(full$A1.freq.f==1), which(full$A2.freq.f==0))
f.a2=intersect(which(full$A2.freq.f==1), which(full$A1.freq.f==0))
fm.a1f=intersect(hetm, f.a1);fm.a2f=intersect(hetm, f.a2)

dd=full[,-c(18:22)]
dd$A1.freq.m[fm.a1f]=0; dd$A2.freq.m[fm.a1f]=1 #change the rows [neoX=A1] of A1 freq for NeoY to 0, A2 freq for NeoY to 1
dd$A1.freq.m[fm.a1f]=0; dd$A2.freq.m[fm.a1f]=1 #change the rows [neoX=A1] of A1 freq for NeoY to 0, A2 freq for NeoY to 1

dd$A1.freq.m=full$A1.freq.m*full$A1.freq.alb
dd$A2.freq.m=full$A2.freq.m*full$A2.freq.alb
dd$A1.freq.f=full$A1.freq.f*full$A1.freq.alb
dd$A2.freq.f=full$A2.freq.f*full$A2.freq.alb

# write.csv(dd, "alb03.male.female.nas00.csv", row.names=F)
