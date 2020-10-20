#last update: Oct.19, 2020

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
posterior.cutoff=0.7 #only take the ancestry genotype call when posterior probability is greater than the cutoff
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
# write.csv(spd, "alb03Xnas00.hybrids.3anc.haplotype.mullercd.newoldlib.sep2020.csv")
# write.csv(bgd, "alb03Xnas00.hybrids.3anc.backgroundinfo.mullercd.newoldlib.sep2020.csv")

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

#STEP 1.2 filtration: only keeping individulas with > cutoff percent of sites inferred
#take the individuals with > 80% sites
cutoff=0.5
cnsite=cutoff*nrow(spd)
list={};list.name={}
for(i in 1:ncol(spd))
{nsite=length(na.omit(spd[,i]));
		if(nsite>cnsite)
		{print(nsite/nrow(spd))
		#print(na.omit(dd[,i]))
			list=c(list, i)
			list.name=c(list.name, nnn[i])}
}

#step 2 impute chromosome
impute.chr=function(mdc, dt, bin)
{dt.mdc.27=dt
bins=round((mdc[length(mdc)]-mdc[1])/bin)+1
start=mdc[1]
#jpeg("check.imput.eg.pdf")
#plot(mdc, dt.mdc.27[67,])
for(b in 1:bins)
	{end=start+bin
	cls=intersect(which(mdc>(start-1)),which(mdc<end))
	#abline(v=end, col="red")
	for(i in 1:nrow(dt.mdc.27))
			{sub=Mode(as.numeric(dt.mdc.27[i,cls]), na.rm=T)[1]
			nas=which(is.na(dt.mdc.27[i,cls]))
			if(length(nas)>0)
				{dt.mdc.27[i,cls][nas]=sub; print(sub)}
			}
	start=start+bin
	}; #dev.off() 
return(dt.mdc.27)}
#remove DBKCL_N701-S507_S17_L002 and DBKCL_N701-S508_S18_L002 due to excessive recomb #rm=which(bgdt$mullercd.rec>30)
list=list[-c(6,7)] #becareful that you only run this once
list.name=list.name[-c(6,7)]
dt=t(spd[,list]) #genotype matrix
bgdt=bgd[list,]#background info
rownames(dt)=list.name
dt.mdc=dt[,which(dat.sp$chrom=="Muller_DC")]
mdc=dat.sp$pos[which(dat.sp$chrom=="Muller_DC")]
dt.mdc.n=impute.chr(mdc, dt.mdc, 2000000) #binsize, the smaller it is, the less imputation
#heatmap(as.matrix(dt.mdc.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
#heatmap(as.matrix(na.omit(dt.mdc.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
sln=sort(list.name)
fln=sln[c(141:148, 122:140,1:121)] ########NEED TO UPDATE IF CHANGE CUTOFFS, AS INDIVIDUAL ORDER WILL CHANGE 
mln=sln[c(203:210, 185:202,149:184)]
jpeg(filename="newoldlib.mullercd.3ance.male.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
#heatmap(as.matrix(dt.mdc.n[mln,]), Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.5)
heatmap(as.matrix(na.omit(dt.mdc.n[rev(mln),])), Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.5)
dev.off()

jpeg(filename="newoldlib.mullercd.3ance.female.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
#heatmap(as.matrix(dt.mdc.n[fln,]), Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.5)
heatmap(as.matrix(na.omit(dt.mdc.n[rev(fln),])), Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.5)
dev.off()

#Step3 count ancestry switches
###count ancestry turnovers in mullercd
sp.ancestry.recomb=function(indv, muller)
{ #ancestry at muller cd
indv.an=spd[which(dat.sp$chrom==muller),which(nnn==indv)]
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

sp.ancestry.sum=function(indv, muller)
{ #ancestry at muller cd
indv.an=spd[which(dat.sp$chrom==muller),which(nnn==indv)]
t=table(indv.an)/sum(table(indv.an))
block=rep(0, 6); names(block)=seq(1:6)
for(n in names(t))
	{block[n]=t[n]} #block[n]=t[n]
return(block)
}

anc.switch=function(sf.plot){
mullera.rec={};mullerb.rec={};mullercd.rec={};mullere.rec={};mullerf.rec={};gen={}; p={}
for(i in sf.plot)
{mullercd.rec=c(mullercd.rec,(nrow(sp.ancestry.recomb(i, "Muller_DC"))-1))
gen=c(gen, bgd$Generation[which(i==nnn)])
p=c(p, as.character(bgd$prefix[which(i==nnn)]))
}
dd.f=data.frame(sf.plot,p, gen, mullercd.rec)
colnames(dd.f)=c("names", "prefix","generation", "muller.CD.rec")
return(dd.f)}
dd.f=anc.switch(fln)
dd.m=anc.switch(mln)

#step 4: haplotype fraction
anc.sum=function(sf.plot){
d={}; gen={}; p={}
for(i in sf.plot)
{d=rbind(d,sp.ancestry.sum(i, "Muller_DC"))
gen=c(gen, bgd$Generation[which(i==nnn)])
p=c(p, as.character(bgd$prefix[which(i==nnn)]))
}
colnames(d)=c("albX.albX", "albX.albY", "albX.nas","albY.albY","albY.nas", "nas.nas")
sumd=data.frame(sf.plot, p, gen,d)
return(sumd)}
sumd.f=anc.sum(fln)
sumd.m=anc.sum(mln)

d.f=data.frame(sumd.f, dd.f)
d.m=data.frame(sumd.m, dd.m)
write.csv(d.f, "sum.switch.oldnewlib.mullercd.3anc.female.csv")
write.csv(d.m, "sum.switch.3anc.mullerced.oldnewlib.male.csv")



