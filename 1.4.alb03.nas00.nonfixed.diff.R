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
#choose the informative sites with allele freq cutoff being 0.3
sp.rows=intersect(which(abs(alb$G1_2-nas$G1_2)>0.3), which(abs(alb$G2_2-nas$G2_2)>0.3))
sp.pos=paste(as.character(alb$CHROM[sp.rows]),as.character(alb$POS[sp.rows]), sep=".")
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

#plot(albd$G1_2[which(albd$CHROM=="Muller_DC")[1:1000]], nasd$G1_2[which(nasd$CHROM=="Muller_DC")[1:1000]], pch=16, col=rgb(0,1,0, 0.01))
#write.csv(d, "alb03.nas00.diffs0.3.csv", row.names=F)
ff=d
ff.dc=ff[which(ff$CHROM=="Muller_DC"),]
sum(ff.dc$POS %in% pos.dc)
sum(ff.dc$POS %in% fm.dc$POS)
sum(ff.dc$POS %in% dat.fm.dc$position)

#table(ff$CHROM )

d$alb.allele=2
d$alb.allele[which((d$A1.freq.alb-d$A1.freq.nas)>0.3)]=1 #allele A1 is albomicans allele

#################Add female and males
#male vs females
albm=read.csv("alb.m.DBMN30-16.19.frq", sep="\t", row.names=NULL)
albf=read.csv("alb.f.DBMN21-D_S4_L007.frq", sep="\t", row.names=NULL)
albm=colnames.change(albm);albf=colnames.change(albf)
albm=text.to.col(albm);albf=text.to.col(albf);
#take out only sp-informative rows
albmm=as.matrix(albm);albfm=as.matrix(albf)
rownames(albmm)=paste(albm$CHROM, albm$POS, sep="."); rownames(albfm)=paste(albf$CHROM, albf$POS, sep=".")

albm.sp=data.frame(albmm[sp.pos,]); albf.sp=data.frame(albfm[sp.pos,]) ##take species-informative rows 

albm.sp$G1_2=as.numeric(as.character(albm.sp$G1_2))
albf.sp$G1_2=as.numeric(as.character(albf.sp$G1_2))

albm.sp$G2_2=as.numeric(as.character(albm.sp$G2_2))
albf.sp$G2_2=as.numeric(as.character(albf.sp$G2_2))
# length(intersect(which(albm.sp$G1_2==0.5), which(abs(albf.sp$G1_2-albf.sp$G2_2)==1)))  # neosex chr informative-- het in males and homo in females
# length(which(albm.sp$G1_2==0.5));length(intersect(which(albm.sp$G1_2==0.5),which(albm.sp$CHROM=="Muller_DC")))
# length(which(abs(albf.sp$G1_2-albf.sp$G2_2)==1));length(intersect(which(abs(albf.sp$G1_2-albf.sp$G2_2)==1),which(albf.sp$CHROM=="Muller_DC")) )
fm.rows=intersect(which(albm.sp$G1_2==0.5),which(abs(albf.sp$G1_2-albf.sp$G2_2)==1))
nas.sp=ff[,8:12]; alb.sp=ff[,3:7]
alb.neo=data.frame(albf.sp[fm.rows,], albm.sp[fm.rows,])
alb.fm=alb.sp[fm.rows,]
alb.neo=data.frame(alb.neo[,c(1:2, 4:8,12:16)], nas.sp[fm.rows,])
colnames(alb.neo)[1:12]=c("CHROM", "POS", "N_CHR.f",	"ALLELE1.f"	,"A1.freq.f",	"ALLELE2.f",	"A2.freq.f",	"N_CHR.m",	"ALLELE1.m",	"A1.freq.m",  "ALLELE2.m", "A2.freq.m")
neo=alb.neo
convert.num=function(d, colnames)
{for(n in colnames)
	{j=which(colnames(d)==n)
	d[,j]=as.numeric(as.character(d[,j]))}
return(d)
}
alb.neo=convert.num(alb.neo, c("A1.freq.f", "A2.freq.f","A1.freq.m","A2.freq.m"))

# neo$A1.freq.f= alb.neo$A1.freq.f*alb.fm$A1.freq.alb
# neo$A2.freq.f= alb.neo$A2.freq.f*alb.fm$A2.freq.alb
# neo$A1.freq.f= alb.neo$A1.freq.m*alb.fm$A1.freq.alb
# neo$A2.freq.f= alb.neo$A2.freq.m*alb.fm$A2.freq.alb

# table(neo$CHROM)
colSums(data.frame(neo$A1.freq.f, neo$A2.freq.f))
# #substitute the neo-Y  specific allele freq (the one different from the homo female version)
alb.neo$A1.freq.m[which(alb.neo$A1.freq.f==1)]=0; #if female homo for A1, male A1 gets assigned to 0
alb.neo$A2.freq.m[which(alb.neo$A1.freq.f==1)]=1 #if female homo for A1, male A2 gets assigned to 1
alb.neo$A1.freq.m[which(alb.neo$A1.freq.f==0)]=1;#if female homo for A2, male A1 gets assigned to 1 
alb.neo$A2.freq.m[which(alb.neo$A1.freq.f==0)]=0#if female homo for A2, male A2 gets assigned to 0 
# alb.neo$f.allele=2 #start assuming female allele is 2
# alb.neo$f.allele[which(alb.neo$A1.freq.f==1)]=1 #if female is homozygous for a1
# write.csv(alb.neo, "alb03.male.female.nas00.csv", row.names=F)


# head(alb.neo)
# sum(alb.neo$POS %in% ff$POS)
# rownames(alb.neo)


# #check that female and male allele freq should sum to 1 for the same allele
# sum(alb.neo$A1.freq.m+alb.neo$A1.freq.f)==nrow(alb.neo)
# sum(alb.neo$A2.freq.m+alb.neo$A2.freq.f)==nrow(alb.neo)
# #check that A1 and A2 allele freq should sum to 1 within female or within male
# sum(alb.neo$A2.freq.f+alb.neo$A1.freq.f)==nrow(alb.neo)
# sum(alb.neo$A2.freq.m+alb.neo$A1.freq.m)==nrow(alb.neo)
# head(alb.neo)

# alb.neo.cd=alb.neo[which(alb.neo$CHROM=="Muller_DC"),]
# d.cd=d[which(d$CHROM=="Muller_DC"),]
# d.cd$POS=as.numeric(as.character(d.cd$POS))
# #plot 
# d.cd$fm.allele=NA
# names(d.cd$fm.allele)=d.cd$POS
# d.cd$fm.allele[alb.neo.cd$POS]=alb.neo.cd$f.allele
# heatmap(t(d.cd) , Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.3)







# fm=read.csv( "alb03.male.female.csv")
# fm.dc=fm[which(fm$CHROM=="Muller_DC"),]
# sum(fm.dc$POS %in% dat.fm.dc$position)
# #Check about Muller B
# albm.mulB=albm[intersect(sp.rows, which(albm$CHROM=="Muller_B")),];albf.mulB=albf[intersect(sp.rows, which(albf$CHROM=="Muller_B")),]; 
# length(intersect(which(albm.mulB$G1_2==0.5), which(abs(albf.mulB$G1_2-albf.mulB$G2_2)==1))) #species-informative && neosex chr informative
# length(which(albm.mulB$G1_2==0.5))
# length(which(abs(albf.mulB$G1_2-albf.mulB$G2_2)==1))
# rows=intersect(which(albm.mulB$G1_2==0.5),which(abs(albf.mulB$G1_2-albf.mulB$G2_2)==1))

# #check about E
# albm.mulE=albm[intersect(sp.rows, which(albm$CHROM=="Muller_E")),];albf.mulE=albf[intersect(sp.rows, which(albf$CHROM=="Muller_E")),]; 
# length(intersect(which(albm.mulE$G1_2==0.5), which(abs(albf.mulE$G1_2-albf.mulE$G2_2)==1))) #species-informative && neosex chr informative
# length(which(albm.mulE$G1_2==0.5))
# length(which(abs(albf.mulE$G1_2-albf.mulE$G2_2)==1))
# dim(albm.mulE)
# rows=intersect(which(albm.mulE$G1_2==0.5),which(aEs(albf.mulE$G1_2-albf.mulE$G2_2)==1))

# #muller A
# albm.mulA=albm[intersect(sp.rows, which(albm$CHROM=="Muller_A")),];albf.mulA=albf[intersect(sp.rows, which(albf$CHROM=="Muller_A")),]; 
# length(intersect(which(albm.mulA$G1_2==0.5), which(abs(albf.mulA$G1_2-albf.mulA$G2_2)==1))) #species-informative && neosex chr informative
# length(which(albm.mulA$G1_2==0.5))
# length(which(abs(albf.mulA$G1_2-albf.mulA$G2_2)==1))
# dim(albm.mulA)
# rows=intersect(which(albm.mulA$G1_2==0.5),which(abs(albf.mulA$G1_2-albf.mulA$G2_2)==1))

