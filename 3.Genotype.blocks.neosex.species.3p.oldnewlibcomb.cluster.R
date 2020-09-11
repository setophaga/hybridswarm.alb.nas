

# install.packages("https://cran.r-project.org/bin/macosx/contrib/4.0/fpc_2.2-5.tgz",method="libcurl")
# install.packages("~/Downloads/mclust_5.4.6.tar")
#
#load('.RData')
library(readxl)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DescTools)
#library(fpc)
library("RColorBrewer")
library(qlcMatrix)
# refd=read.csv("../alb03.male.female.nas00.csv")#on laptop: refd=read.csv("~/Desktop/hybridswarm/1.pipeline/parentalstrain/alb03.male.female.nas00.csv")
#STEP 1 get ancestry-specific genotype from AncestryHMM output
files <- list.files(path = "/scratch/silu/abo.nas/x1.old.newlib.bam.kepmaskrep/posterior.sep.10th", pattern = "*.posterior", full.names = T)
strsplit(strsplit(files[1],"/")[[1]][7], ".posterior")[[1]][1] 
#files <- list.files(path = "~/Desktop/hybridswarm/1.pipeline/posterior.alb03.nas00.all.July22th", pattern = "*.posterior", full.names = T)
#use this above line if run on laptop
e.sp=read.csv(files[1], sep="\t");e.sp
posterior.cutoff=0.9 #only take the ancestry genotype call when posterior probability is greater than the cutoff
names={}; spd={}; h={}; 
for(ind in files)
{name=strsplit(strsplit(ind,"/")[[1]][7], ".posterior")[[1]][1] #use [8] IF run on laptop
names=c(names, name)
 d=read.csv(ind, sep="\t")
v=rep(NA,nrow(d))
an=Matrix(as.matrix(d[,3:8]))
#if the posteriors are less than cutoff, don't consider them
maxresult=rowMax(an, which=T)
anm=maxresult$which
v[which(anm[,1]==1)]=1 #albX, albX
v[which(anm[,2]==1)]=2 #albX, albY
v[which(anm[,3]==1)]=3 #albX, nas
v[which(anm[,4]==1)]=4 #albY, albY
v[which(anm[,5]==1)]=5 #albY, nas
v[which(anm[,6]==1)]=6 #nas, nas
print(maxresult$max[which(maxresult$max<posterior.cutoff)])
if(length(which(maxresult$max<posterior.cutoff))>0)
	{v[which(maxresult$max<posterior.cutoff)]=NA }
if(sum(is.na(maxresult$max))>0)
	{v[is.na(maxresult$max)]=NA }
spd=cbind(spd, v)
}
dat.sp=data.frame(e.sp[,1:2], spd)
dat.sp=dat.sp[order(dat.sp$chrom, dat.sp$position),]
table(dat.sp$chrom)
spd=dat.sp[,13:ncol(dat.sp)]

d=read.csv("alb03.nas00.X1.newold.lib.252.csv") #match generation information for each  
 bg={}; 
for( i in 11:length(names))
	{row=which(names[i]==d$prefix)
		if(length(row)==0){print(names[i]); bg=rbind(bg,NA)}
	bg=rbind(bg, d[row,])}
bgd=data.frame(names[11:length(names)], bg)
length((colMeans(spd, na.rm=T)))
#correct for sex mis-identification
nnn=paste( bgd$sex.g,bgd$Generation,bgd$Species, names[11:length(names)],sep=".")
colnames(spd)=nnn
write.csv(spd, "alb03Xnas00.hybrids.3anc.haplotype.mullercd.newoldlib.sep2020.csv")
write.csv(bgd, "alb03Xnas00.hybrids.3anc.backgroundinfo.mullercd.newoldlib.sep2020.csv")


sn=sort(nnn) #sort by generation id
sf=sn[1:154]; sm=sn[154:229]
sf.plot=sf[c(147:154,127:146, 1:126)]
sm.plot=sm[c(68:76, 50:67,1:49 )]
spd.f=spd[,rev(sf.plot)]
spd.m=spd[,rev(sm.plot)]
col.sp <- colorRampPalette(c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"))(6)

plot(x=NULL, y=NULL)
legend("topright", c( "1=(albX, albX)", "2=(albX, albY)", "3=(albX, nas)","4=(albY, albY)","5=(albY, nas)", "6=(nas, nas)"), fill=c("coral4","deepskyblue3", "gold", "chartreuse3", "lightslateblue", "aquamarine2"))

jpeg(filename="newoldlib.female.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[,rev(sf.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

jpeg(filename="newoldlib.male.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[,rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

jpeg(filename="newoldlib.male.muller.cd.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[which(dat.sp$chrom=="Muller_DC"),rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

jpeg(filename="newoldlib.female.muller.cd.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[which(dat.sp$chrom=="Muller_DC"),rev(sf.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()


