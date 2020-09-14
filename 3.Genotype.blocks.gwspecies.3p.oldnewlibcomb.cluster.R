

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
files <- list.files(path = "/scratch/silu/abo.nas/x1.old.newlib.bam.kepmaskrep/posterior.sep.13th", pattern = "*.posterior", full.names = T)
strsplit(strsplit(files[1],"/")[[1]][7], ".posterior")[[1]][1] 
#files <- list.files(path = "~/Desktop/hybridswarm/1.pipeline/posterior.alb03.nas00.all.July22th", pattern = "*.posterior", full.names = T)
#use this above line if run on laptop
e.sp=read.csv(files[1], sep="\t");head(e.sp)
posterior.cutoff=0.9 #only take the ancestry genotype call when posterior probability is greater than the cutoff
names={}; spd={}; h={}; 
for(ind in files)
{name=strsplit(strsplit(ind,"/")[[1]][7], ".posterior")[[1]][1] #use [8] IF run on laptop
names=c(names, name)
 d=read.csv(ind, sep="\t")
v=rep(NA,nrow(d))
an=Matrix(as.matrix(d[,3:5]))
#if the posteriors are less than cutoff, don't consider them
maxresult=rowMax(an, which=T)
anm=maxresult$which
v[which(anm[,1]==1)]=1 #alb, alb
v[which(anm[,2]==1)]=0.5 #alb, nas
v[which(anm[,3]==1)]=0 #nas, nas
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
#calculate genome-wide and muller element-specific heterzygosity
spd.het=spd #assigned the ancestry matrix and turn the non-het sites into 0 and het sites into1
spd.het=replace(spd.het, spd.het == 1, 0)
spd.het=replace(spd.het, spd.het == 0.5, 1)
table(spd.het[,19]); table(spd[,19]); #check
bgd$het=colMeans(spd.het, na.rm=T)
bgd$HI=colMeans(spd, na.rm=T)
#take colmeans as genome-wide het
#take different muller element and calc col means
bgd$het.a=colMeans(spd.het[which(dat.sp$chrom=="Muller_A"),], na.rm=T)
bgd$hi.a=colMeans(spd[which(dat.sp$chrom=="Muller_A"),], na.rm=T)

bgd$het.b=colMeans(spd.het[which(dat.sp$chrom=="Muller_B"),], na.rm=T)
bgd$hi.b=colMeans(spd[which(dat.sp$chrom=="Muller_B"),], na.rm=T)

bgd$het.dc=colMeans(spd.het[which(dat.sp$chrom=="Muller_DC"),], na.rm=T)
bgd$hi.dc=colMeans(spd[which(dat.sp$chrom=="Muller_DC"),], na.rm=T)

bgd$het.e=colMeans(spd.het[which(dat.sp$chrom=="Muller_E"),], na.rm=T)
bgd$hi.e=colMeans(spd[which(dat.sp$chrom=="Muller_E"),], na.rm=T)

bgd$het.f=colMeans(spd.het[which(dat.sp$chrom=="Muller_F"),], na.rm=T)
bgd$hi.f=colMeans(spd[which(dat.sp$chrom=="Muller_F"),], na.rm=T)

write.csv(spd, "alb03Xnas00.hybrids.3anc.haplotype.gw.newoldlib.sep2020.csv")


sn=sort(nnn) #sort by generation id
sf=sn[1:154]; sm=sn[154:229]
sf.plot=sf[c(147:154,127:146, 1:126)]
sm.plot=sm[c(68:76, 50:67,1:49 )]
spd.f=spd[,rev(sf.plot)]
spd.m=spd[,rev(sm.plot)]

col.sp <- colorRampPalette(c("turquoise","royalblue4"))(3)

jpeg(filename="newoldlib.gw.sp.female.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[,rev(sf.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

jpeg(filename="newoldlib.gw.sp.male.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[,rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

jpeg(filename="newoldlib.gw.sp.male.muller.cd.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[which(dat.sp$chrom=="Muller_DC"),rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

jpeg(filename="newoldlib.gw.sp.female.muller.cd.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[which(dat.sp$chrom=="Muller_DC"),rev(sf.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()


###count ancestry turnovers 
sp.ancestry.recomb=function(indv, muller)
{ indv.an=spd[which(dat.sp$chrom==muller),which(nnn==indv)]
positions=dat.sp[which(dat.sp$chrom==muller),2]
pos.notna=positions[which(is.na(indv.an)==0)]
an.notna=na.omit(indv.an)
pos.switch=pos.notna[which(abs((an.notna-c(an.notna[-1],an.notna[1])))>0)]
if(length(pos.switch)==0)
{if(length(an.notna)==0){anblock=data.frame(ancestry=NA,begin=NA,end=NA)}
else{anblock=data.frame(ancestry=an.notna[1], begin=pos.notna[1], end=tail(pos.notna, 1))}}
else{ancestry=an.notna[1]; begin=pos.notna[1]; end=pos.switch[1]; lastswitch=tail(pos.switch,1)
	for(p in pos.switch)
		{if(p==lastswitch)
			{if(length(which(pos.notna>lastswitch))>0)
				{ancestry=c(ancestry, an.notna[which(pos.notna==p)+1])#start from the second block
				begin=c(begin,pos.notna[which(pos.notna==p)+1])
				end=c(end, tail(pos.notna, 1))}
			else{break}}
		else #make a data.frame containing begin and end position of ancestry blocks
			{ancestry=c(ancestry, an.notna[which(pos.notna==p)+1])#start from the second block
			begin=c(begin,pos.notna[which(pos.notna==p)+1])
			end=c(end, (pos.switch[which(pos.switch==p)+1]+1))}
		};anblock=data.frame(ancestry, begin, end)}	
return(anblock)
}

for(i in 1:length(nnn))
{bgd$mullercd.rec[i]=nrow(sp.ancestry.recomb(nnn[i], "Muller_DC"))
bgd$mullera.rec[i]=nrow(sp.ancestry.recomb(nnn[i], "Muller_A"))
bgd$mullerb.rec[i]=nrow(sp.ancestry.recomb(nnn[i], "Muller_B"))
bgd$mullere.rec[i]=nrow(sp.ancestry.recomb(nnn[i], "Muller_E"))
bgd$mullerf.rec[i]=nrow(sp.ancestry.recomb(nnn[i], "Muller_F"))
}
write.csv(bgd, "alb03Xnas00.hybrids.3anc.backgroundinfo.gw.newoldlib.sep2020.csv")



