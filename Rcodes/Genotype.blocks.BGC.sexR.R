

# install.packages("https://cran.r-project.org/bin/macosx/contrib/4.0/fpc_2.2-5.tgz",method="libcurl")
# install.packages("~/Downloads/mclust_5.4.6.tar")

setwd("~/Desktop/hybridswarm/AncestryHHM/output.Mar21.oldlib.alb.nas")

library(readxl)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fpc)
#STEP 1 get ancestry-specific genotype from AncestryHMM output
files <- list.files(path = "~/Desktop/hybridswarm/AncestryHHM/output.Mar21.oldlib.alb.nas", pattern = "*.posterior", full.names = T)
e=read.csv(files[1], sep="\t")
posterior.cutoff=0.6 #only take the ancestry genotype call when posterior probability is greater than the cutoff
names={}; dd={}; h={}; 
for(ind in files)
{name=strsplit(strsplit(ind,"/")[[1]][8], ".posterior")[[1]][1]
names=c(names, name)
d=read.csv(ind, sep="\t")
v=rep(NA, length.out=length(d[,1]))
an=c(d[,3], d[,4], d[,5])
v[which(d[,5]>posterior.cutoff)]=0
v[which(d[,4]>posterior.cutoff)]=0.5
v[which(d[,3]>posterior.cutoff)]=1
dd=cbind(dd, v)
}
dat=data.frame(e[,1:2], dd)
library("RColorBrewer")
col <- colorRampPalette(c("turquoise","royalblue4"))(100)
#heatmap(t(na.omit(dd)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.3)
 d=read.csv("~/Desktop/hybridswarm/nasutaXabomincans.all.simple.csv") #match generation information for each individual
 d$Gender[which(d$Gender=="F")]="female"; d$Gender=factor(d$Gender, levels=c("female", "male", "ND"));table(d$Gender)
 d28=d[which(d$Generation==28),];table(d$Gender,d$Generation)
bg={}
for( i in 1:length(names))
	{row=which(names[i]==d$prefix)
	bg=rbind(bg, d[row,])
	}
bgd=data.frame(names, bg)
nnn=paste(bgd$Generation, bgd$Gender,names,bgd$Species, sep=".")
colnames(dd)=nnn
sn=sort(nnn)
heatmap(t(dd[,nnn]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.3)

#STEP 1.2 filtration: only keeping individulas with > cutoff percent of sites inferred
#take the individuals with > 50% sites
cutoff=0.5
cnsite=cutoff*nrow(dd)
list={}
for(i in 1:ncol(dd))
{nsite=length(na.omit(dd[,i]));
		if(nsite>cnsite)
		{print(nsite/nrow(dd))
		#print(na.omit(dd[,i]))
			list=c(list, i)}
}
#heatmap(t(as.matrix(dd[,list])), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
#take the subset of individuals with >50% coverage
dt=t(dd[,list]) #genotype matrix
bgdt=bgd[list,]#background info
dt.mb=dt[,which(dat$chrom=="Muller_B")]
dt.mdc=dt[,which(dat$chrom=="Muller_DC")]
dt.me=dt[,which(dat$chrom=="Muller_E")]
dt.mf=dt[,which(dat$chrom=="Muller_F")]
dt.ma=dt[,which(dat$chrom=="Muller_A")]
mullers=c(rep("Muller_A", ncol(dt.ma)),rep("Muller_DC", ncol(dt.mdc)), rep("Muller_B", ncol(dt.mb)),rep("Muller_E", ncol(dt.me)),rep("Muller_F", ncol(dt.mf)))
###seperate each muller element
raw.female.gen27=intersect(which(bgdt$Gender=="female"), which(bgdt$Generation==27))
raw.male.gen27=intersect(which(bgdt$Gender=="male"), which(bgdt$Generation==27))
dtt=cbind(dt.ma[,rev(seq(1:ncol(dt.ma)))],dt.mdc, dt.mb[,rev(seq(1:ncol(dt.mb)))], dt.me, dt.mf )
dt.female.gen27=dtt[raw.female.gen27,]
dt.male.gen27=dtt[raw.male.gen27,]
het.female.gen27=matrix(0, nrow(dt.female.gen27), ncol(dt.female.gen27))
het.female.gen27[which(dt.female.gen27==0.5)]=1
het.male.gen27=matrix(0, nrow(dt.male.gen27), ncol(dt.male.gen27))
het.male.gen27[which(dt.male.gen27==0.5)]=1
#heterozygosity difference between female and male for all muller element
col.se=function(m) #a function return SE for each locus, either heterozygosity or ancestry
{ses={}
for(j in 1:ncol(m))
	{ses=c(ses, se(m[,j]))}
	return(ses)}
col.mean=function(m) #a function return SE for each locus, either heterozygosity or ancestry
{ms={}
for(j in 1:ncol(m))
	{ms=c(ms, mean(m[,j], na.rm=T))}
	return(ms)}


t.test(rowMeans(het.female.gen27),rowMeans(het.male.gen27))
t.test(colMeans(het.female.gen27),colMeans(het.male.gen27), paired=T)
hist(rowMeans(macd.het.female), breaks=20, xlim=c(0,1), ylim=c(0,6), col=rgb(1,0,0.2, 0.2), main="all chr", xlab="Heterozygosity")
hist(rowMeans(macd.het.male), breaks=20, add=T, col=rgb(0,0.1,1, 0.2))
#heterozygosity difference between male and female
hist(colMeans(het.female.gen27,na.rm=T), breaks=20, xlim=c(0,1),col=rgb(1,0,0.2, 0.2),ylim=c(0, 20000), main="", xlab="Heterozygosity")
hist(colMeans(het.male.gen27,na.rm=T), breaks=10, add=T, col=rgb(0,0.1,1, 0.2))
t.test(colMeans(het.female.gen27,na.rm=T),colMeans(het.male.gen27,na.rm=T),paired=T)
#a bootstrap difference between female and male
boot.diff=function(v1, v2,n,  N)
{m.diff={}
for(i in 1:N)
{length=seq(1:length(v1))
pos=sample(length,  replace=T, n)
s1=v1[pos]	; s2=v2[pos]
m.diff=c(m.diff, mean((s1-s2), na.rm=T))}
return(m.diff)}
het.mean.f.27=colMeans(het.female.gen27)
het.mean.m.27=colMeans(het.male.gen27)
nomulleradc.diff=boot.diff(het.mean.f.27[-c(which(mullers=="Muller_A"),which(mullers=="Muller_DC"))],het.mean.m.27[-c(which(mullers=="Muller_A"),which(mullers=="Muller_DC"))],n=round(length(het.mean.f.27[-c(which(mullers=="Muller_A"),which(mullers=="Muller_DC"))])*1), N=10000)
nomullera.diff=boot.diff(het.mean.f.27[-c(which(mullers=="Muller_A"))],het.mean.m.27[-c(which(mullers=="Muller_A"))],n=round(length(het.mean.f.27[-c(which(mullers=="Muller_A"))])*1), N=10000)
nomullerdc.diff=boot.diff(het.mean.f.27[-which(mullers=="Muller_DC")],het.mean.m.27[-which(mullers=="Muller_DC")],n=round(length(het.mean.f.27[-which(mullers=="Muller_DC")])*1), N=10000)
quantile(nomullerdc.diff, 0.025)
hist(nomulleradc.diff)
all.diff=boot.diff(het.mean.f.27,het.mean.m.27,n=round(length(het.mean.f.27)*1), N=10000)
quantile(nomullerdc.diff, 0.025)
hist(nomulleradc.diff, border=rgb(0,0.87,1, 0.5), xlim=c(-0.007, 0.05), breaks=200, ylim=c(0,600), xlab="Heterozygosity F-M", main="")
hist(nomullera.diff, border=rgb(0,0.73,0.1, 0.5), add=T, breaks=100)
hist(nomullerdc.diff, border=rgb(0.4,0.1,0.43, 0.5), add=T, breaks=100)
hist(all.diff, add=T, border=rgb(0.9,0.82,0, 0.5), breaks=100)
legend("topleft",  c("no muller A","no muller DC","no muller A & DC", "all"), fill=c(rgb(0,0.73,0.1, 1),rgb(0.4,0.1,0.43, 1),rgb(0,0.87,1, 01), rgb(0.9,0.82,0, 1)), cex=0.5)

##plot heterozygosity along all muller 
plot(seq(1:ncol(het.female.gen27)), colMeans(het.female.gen27), col=rgb(0.8,0,0.3, 0.01), pch=16, ylab="Heterozygosity", ylim=c(0, 0.8))
lines(seq(1:ncol(het.female.gen27)), colMeans(het.female.gen27),col=rgb(0.8,0,0.3, 0.5))
arrows(seq(1:ncol(het.female.gen27)), colMeans(het.female.gen27)-col.se(het.female.gen27), seq(1:ncol(het.male.gen27)),colMeans(het.female.gen27)+col.se(het.female.gen27), length=0.05, angle=90, code=3, col=rgb(0.8,0,0.3, 0.01))
abline(v=ncol(dt.ma))
points(seq(1:ncol(het.male.gen27)), colMeans(het.male.gen27),  col=rgb(0,0.1,0.5, 0.01), pch=16)
lines(seq(1:ncol(het.male.gen27)), colMeans(het.male.gen27),col=rgb(0,0.1,0.5, 0.5))
abline(v=ncol(dt.ma));abline(v=(ncol(dt.mdc)+ncol(dt.ma)));
abline(v=(ncol(dt.mdc)+ncol(dt.ma)+ncol(dt.mb)));abline(v=(ncol(dt.me)+ncol(dt.mdc)+ncol(dt.ma)+ncol(dt.mb)))
arrows(seq(1:ncol(het.male.gen27)), colMeans(het.male.gen27)-col.se(het.male.gen27),seq(1:ncol(het.male.gen27)), colMeans(het.male.gen27)+col.se(het.male.gen27), length=0.05, angle=90, code=3, col=rgb(0,0.1,0.5, 0.01))
abline(v=ncol(dt.ma))
##plot ancestry along all muller 
y.f=col.mean(dt.female.gen27);x=seq(1:ncol(dt.female.gen27));y.m=col.mean(dt.male.gen27);
plot(x, y.f, col=rgb(0.8,0,0.3, 0.01), pch=16, ylab="Ancestry", ylim=c(0, 1))
lines(x, y.f,col=rgb(0.8,0,0.3, 0.5))
arrows(x, y.f-col.se(dt.female.gen27), seq(1:ncol(dt.male.gen27)),y.f+col.se(dt.female.gen27), length=0.02, angle=90, code=3, col=rgb(0.8,0,0.3, 0.01))
points(x, y.m,  col=rgb(0,0.1,0.5, 0.01), pch=16)
lines(x, y.m,col=rgb(0,0.1,0.5, 0.5))
abline(v=ncol(dt.ma));abline(v=(ncol(dt.mdc)+ncol(dt.ma)));
abline(v=(ncol(dt.mdc)+ncol(dt.ma)+ncol(dt.mb)));abline(v=(ncol(dt.me)+ncol(dt.mdc)+ncol(dt.ma)+ncol(dt.mb)))
arrows(x,y.m-col.se(dt.male.gen27),x, y.m+col.se(dt.male.gen27), length=0.02, angle=90, code=3, col=rgb(0,0.1,0.5, 0.01))
abline(v=ncol(dt.ma))
t.test(y.f, y.m, paired=T)
hist(y.f, breaks=20, xlim=c(0,1),col=rgb(1,0,0.2, 0.2),ylim=c(0, 30000), main="", xlab="Ancestry")
hist(y.m, breaks=10, add=T, col=rgb(0,0.1,1, 0.2))

####focus on sex chromosomes
na.mt=matrix(NA, nrow = nrow(dt.ma), ncol = 100)
dt.macd.nagap=cbind(dt.ma[,rev(seq(1:ncol(dt.ma)))],na.mt , dt.mdc)
heatmap(as.matrix(dt.macd.nagap[raw.female.gen27,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="female gen27")
heatmap(as.matrix(dt.macd.nagap[raw.male.gen27,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="male gen27")
dt.macd=cbind(dt.ma[,rev(seq(1:ncol(dt.ma)))],dt.mdc)
macd.het.male=matrix(0, nrow(dt.macd[raw.male.gen27,]), ncol(dt.macd[raw.male.gen27,]))
macd.het.male[which(dt.macd[raw.male.gen27,]==0.5)]=1
macd.het.female=matrix(0, nrow(dt.macd[raw.female.gen27,]), ncol(dt.macd[raw.female.gen27,]))
macd.het.female[which(dt.macd[raw.female.gen27,]==0.5)]=1

#heterozygosity difference between female and male ; onnly muller A and muller CD
t.test(rowMeans(macd.het.female),rowMeans(macd.het.male))
hist(rowMeans(macd.het.female), breaks=20, xlim=c(0,1), ylim=c(0,6), col=rgb(1,0,0.2, 0.2), main="", xlab="Heterozygosity")
hist(rowMeans(macd.het.male), breaks=20, add=T, col=rgb(0,0.1,1, 0.2))
#ancestry difference between male and female
t.test(rowMeans(dt.macd[raw.male.gen27,], na.rm=T),rowMeans(dt.macd[raw.female.gen27,], na.rm=T))
hist(rowMeans(dt.macd[raw.female.gen27,], na.rm=T), breaks=20, xlim=c(0,1), ylim=c(0,6), col=rgb(1,0,0.2, 0.2), main="", xlab="Ancestry")
hist(rowMeans(dt.macd[raw.male.gen27,], na.rm=T), breaks=10, add=T, col=rgb(0,0.1,1, 0.2))
##plot heterozygosity along muller a and muller cd
plot(seq(1:ncol(macd.het.female)), colMeans(macd.het.female), col=rgb(0.8,0,0.3, 0.01), pch=16, ylab="Heterozygosity", ylim=c(0, 0.8))
lines(seq(1:ncol(macd.het.female)), colMeans(macd.het.female),col=rgb(0.8,0,0.3, 0.5))
arrows(seq(1:ncol(macd.het.female)), colMeans(macd.het.female)-col.se(macd.het.female), seq(1:ncol(macd.het.male)),colMeans(macd.het.female)+col.se(macd.het.female), length=0.05, angle=90, code=3, col=rgb(0.8,0,0.3, 0.01))
abline(v=ncol(dt.ma))
points(seq(1:ncol(macd.het.male)), colMeans(macd.het.male),  col=rgb(0,0.1,0.5, 0.01), pch=16)
lines(seq(1:ncol(macd.het.male)), colMeans(macd.het.male),col=rgb(0,0.1,0.5, 0.5))
abline(v=ncol(dt.ma))
arrows(seq(1:ncol(macd.het.male)), colMeans(macd.het.male)-col.se(macd.het.male),seq(1:ncol(macd.het.male)), colMeans(macd.het.male)+col.se(macd.het.male), length=0.05, angle=90, code=3, col=rgb(0,0.1,0.5, 0.01))
abline(v=ncol(dt.ma))
legend(30000,0.2,c("female", "male"), pch=16, col=c(rgb(0.8,0,0.3, 0.5), rgb(0,0.1,0.5, 0.5)), lwd=1, cex=1)
t.test(colMeans(macd.het.female),colMeans(macd.het.male),paired=T)

##plot ancestry along muller a and muller cd
plot(seq(1:ncol(dt.macd)), colMeans(dt.macd[raw.male.gen27,], na.rm=T), col=rgb(0,0.1,0.5, 0.01), pch=16, ylab="Ancestry", ylim=c(0.2, 0.9))
lines(seq(1:ncol(dt.macd)), colMeans(dt.macd[raw.male.gen27,], na.rm=T),col=rgb(0,0.1,0.5, 0.01))
arrows(seq(1:ncol(dt.macd)), colMeans(dt.macd[raw.male.gen27,], na.rm=T)-col.se(dt.macd[raw.male.gen27,]), seq(1:ncol(dt.macd)),  colMeans(dt.macd[raw.male.gen27,], na.rm=T)+col.se(dt.macd[raw.male.gen27,]), length=0.05, angle=90, code=3, col=rgb(0,0.1,0.5, 0.01))
abline(v=ncol(dt.ma))
points(seq(1:ncol(dt.macd)), colMeans(dt.macd[raw.female.gen27,], na.rm=T), col=rgb(0.8,0,0.3, 0.5), pch=16)
lines(seq(1:ncol(dt.macd)), colMeans(dt.macd[raw.female.gen27,], na.rm=T),col=rgb(0.8,0,0.3, 0.5))
arrows(seq(1:ncol(dt.macd)), colMeans(dt.macd[raw.female.gen27,], na.rm=T)-col.se(dt.macd[raw.female.gen27,]), seq(1:ncol(dt.macd)),  colMeans(dt.macd[raw.female.gen27,], na.rm=T)+col.se(dt.macd[raw.female.gen27,]), length=0.05, angle=90, code=3, col=rgb(0.8,0,0.3, 0.01))
abline(v=ncol(dt.ma))
legend(30000,0.5,c("female", "male"), pch=16, col=c(rgb(0.8,0,0.3, 0.5), rgb(0,0.1,0.5, 0.5)), lwd=1, cex=1)
t.test(colMeans(dt.macd[raw.female.gen27,], na.rm=T), colMeans(dt.macd[raw.male.gen27,], na.rm=T),paired=T)





# orders=names(sort(rowSums(dt.mdc.27, na.rm=T)))
# heatmap(as.matrix(dt.mdc.27[orders,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
# # d.ma=cor(t(dt.ma), na.rm=T)
#imput with 200kb windows

###STEP 1.3 chromosome blocks: generate postion m index for each chromosome; and data matrix for all generations
mdc=dat$pos[which(dat$chrom=="Muller_DC")]
ma=dat$pos[which(dat$chrom=="Muller_A")]
mb=dat$pos[which(dat$chrom=="Muller_B")]
me=dat$pos[which(dat$chrom=="Muller_E")]
mf=dat$pos[which(dat$chrom=="Muller_F")]

#step 2 impute chromosome
impute.chr=function(mdc, dt, bin)
{
dt.mdc.27=dt
bins=round((mdc[length(mdc)]-mdc[1])/bin)+1
start=mdc[1]
plot(mdc, dt.mdc.27[27,])
for(b in 1:bins)
	{
	end=start+bin
	cls=intersect(which(mdc>(start-1)),which(mdc<end))
	abline(v=end, col="red")
	for(i in 1:nrow(dt.mdc.27))
		{sub=mean(dt.mdc.27[i,cls], na.rm=T)
		nas=which(is.na(dt.mdc.27[i,cls]))
		if(length(nas)>0)
			{dt.mdc.27[i,cls][nas]=sub}
		}
	start=start+bin
	} ;return(dt.mdc.27)
}

dt.mdc.n=impute.chr(mdc, dt.mdc, 2000000) #binsize, the smaller it is, the less imputation
heatmap(as.matrix(dt.mdc.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
heatmap(as.matrix(na.omit(dt.mdc.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

dt.ma.n=impute.chr(ma, dt.ma, 2000000)
dt.macd.n=cbind(dt.ma.n[,rev(seq(1:ncol(dt.ma.n)))], na.mt, dt.mdc.n)
heatmap(as.matrix(dt.macd.n[raw.female.gen27,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="imputed, female, gen27")
heatmap(as.matrix(dt.macd.n[raw.male.gen27,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="imputed, male, gen27")

dt.mb.n=impute.chr(mb, dt.mb, 2000000)
heatmap(as.matrix(dt.mb.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

dt.me.n=impute.chr(me, dt.me, 2000000)
heatmap(as.matrix(dt.me.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
heatmap(as.matrix(na.omit(dt.me.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

dt.mf.n=impute.chr(mf, dt.mf, 200000)
heatmap(as.matrix(dt.mf.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
heatmap(as.matrix(na.omit(dt.mf.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)


#step 3 cluster each chromosome
na.rows.a=which(is.na(rowMeans(dt.ma.n))==1)
na.rows.dc=which(is.na(rowMeans(dt.mdc.n))==1)
na.rows.b=which(is.na(rowMeans(dt.mb.n))==1)
na.rows.e=which(is.na(rowMeans(dt.me.n))==1)
na.rows.f=which(is.na(rowMeans(dt.mf.n))==1)
heatmap(dt.mdc.n[na.rows,],Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5) #check if they indeed miss 
na.rows=Reduce(union, list(na.rows.a, na.rows.b,na.rows.dc, na.rows.e, na.rows.f))

dc.clean=na.omit(dt.mdc.n)
dim(dc.clean); dim(dt.mdc.n)
heatmap(as.matrix(dc.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_DC")
a.clean=na.omit(dt.ma.n)
dim(a.clean); dim(dt.ma.n)
heatmap(as.matrix(a.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_A")

b.clean=na.omit(dt.mb.n)
dim(b.clean); dim(dt.mb.n)
heatmap(as.matrix(b.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_B")

e.clean=na.omit(dt.me.n)
dim(e.clean); dim(dt.me.n)
heatmap(as.matrix(e.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_E")

f.clean=na.omit(dt.mf.n)
dim(f.clean); dim(dt.mf.n)
heatmap(as.matrix(f.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_F")

#step 4: calculate ancestry for each clusters
cluster.ancestry=function( dc.clean, clusterN)
{set.seed(26)
dck=kmeans(t(dc.clean), centers=clusterN); dck
while(dck$"betweenss"/dck$"totss"<0.95)
	{clusterN=clusterN+1; set.seed(26)
	dck=kmeans(t(dc.clean), centers=clusterN)}
print(dck)
print(clusterN)
#cluster numbers are not sorted by distance now sort by position along the chromosome
##it's okay, we can generate a factor called dc.cls with the ordering consistent with orders along the chromosome
clusters={};  n=0
cluster.new=dck$cluster
for(i in 2:(length(dck$cluster)))
{if(dck$cluster[i-1]!=dck$cluster[i])
	{print(dck$cluster[i])
	if(dck$cluster[i]%in%clusters)
		{n=n+1
		tochange=dck$cluster[i]
		all=as.numeric(as.character(which(dck$cluster==tochange)))
		want=all[c(which(all>i),which(all==i))]
		cluster.new[want]=range(dck$cluster)[2]+n
		# print("ahhh");print(dck$cluster[i])
		print(range(dck$cluster)[2]+n)
		clusters=c(clusters, range(dck$cluster)[2]+n)}
	else{clusters=c(clusters, dck$cluster[i])}
	}
}
print(length(clusters)) #the right order for clusters along the chromosome
print(table(cluster.new)) #see partitions of new cluster information
cd.cl={}
for(i in 1:nrow(dc.clean)) #loop among individuals
	{cl={} #fresh vector containing mean ancestry for each cluster 
	for(c in clusters) #loop among clusters
		{cl=c(cl,mean(dc.clean[i, which(cluster.new==c)]))}
	cd.cl=rbind(cd.cl, cl) #attach cluster ancestry row for each indivudual
	}
colnames(cd.cl)=paste("c", clusters)
rownames(cd.cl)=rownames(dc.clean)
heatmap(as.matrix(cd.cl), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5)
return(cd.cl)}
###END of the function


#run the function for each muller element
dc.clust.an=cluster.ancestry( dc.clean, 1)
a.clust.an=cluster.ancestry( a.clean, 1)
b.clust.an=cluster.ancestry( b.clean, 1)
e.clust.an=cluster.ancestry( e.clean,1)
f.clust.an=cluster.ancestry( f.clean,1)

#calculate HI for each individuals

#find intersect of all the muller.elements (individuals that are present for all)
a.list=row.names(a.clust.an)
b.list=row.names(b.clust.an)
dc.list=row.names(dc.clust.an)
e.list=row.names(e.clust.an)
f.list=row.names(f.clust.an)
all.list.notsort=Reduce(intersect, list(a.list,b.list,dc.list,e.list,f.list))
bg.info.clust=data.frame(name=all.list.notsort, sex=bgdt$Gender[-na.rows], gen=bgdt$Generation[-na.rows],)
bg.info.clust=bg.info.clust[order(bg.info.clust[,1]),]
all.list=bg.info.clust$name
all.sex=bg.info.clust$sex
all.gen=bg.info.clust$Generation
#get genotype clusters for individuals that are present for all muller elements
a.c.an=a.clust.an[all.list,]
b.c.an=b.clust.an[all.list,]
dc.c.an=dc.clust.an[all.list,]
e.c.an=e.clust.an[all.list,]
f.c.an=f.clust.an[all.list,]
all.an=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))],  dc.c.an, b.c.an[,rev(seq(1:ncol(b.c.an)))], e.c.an, f.c.an)
acd.c.an=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))], rep(NA, nrow(a.c.an)), dc.c.an)
##
####.  Figure 1.  !@@@@@
all.an.plot=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))], rep(NA, nrow(a.c.an)), dc.c.an, rep(NA, nrow(a.c.an)),b.c.an[,rev(seq(1:ncol(b.c.an)))],rep(NA, nrow(a.c.an)), e.c.an, rep(NA, nrow(a.c.an)),f.c.an)
all.an.p=rbind(all.an.plot[n28,], rep(NA, ncol(all.an.plot)), all.an.plot[n27,],rep(NA, ncol(all.an.plot)), all.an.plot[n21,], rep(NA, ncol(all.an.plot)), all.an.plot[c(n5,n4,n3, na, nn),])
heatmap(as.matrix(all.an.p), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.3,cexCol=0.3)
#####. @@@@

heatmap(as.matrix(a.clust.an[all.list,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_A")
heatmap(as.matrix(b.clust.an[all.list,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_B")
heatmap(as.matrix(e.clust.an[all.list,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_E")
heatmap(as.matrix(f.clust.an[all.list,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_F")
heatmap(as.matrix(dc.clust.an[all.list,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_DC")
nn=c(100:106,111:112);na=c(98:99,107:110);n3=77; n4=78:90; n5=91:97
n21=2:16;n27=17:47; n28=48:76

#check if the row numbers are right
all.list[n21]; all.list[n27]; #should only print out individuals from generation 21 and 27
all.list[n28]; all.list[nn]
all.sex[n27]; all.sex[n28]; 
#calculate BGC with these
##step1, calculate HI for each individual
#first HI for each muller element
dc.hi=rowMeans(dc.c.an)
a.hi=rowMeans(a.c.an)
b.hi=rowMeans(b.c.an)
e.hi=rowMeans(e.c.an)
f.hi=rowMeans(f.c.an)
dd.hi=data.frame(a.hi, b.hi, dc.hi, e.hi, f.hi)
hi=rowMeans(dd.hi)
##step2, seperate into data for each generation type: early, 21, 27
hi.early=hi[77:112] ; ne=77:112; ge=ne #check all.list[ge]
hi.21.male=hi[intersect(n21, which(all.sex=="male"))]; g21.male=intersect(n21, which(all.sex=="male"));g21.female=intersect(n21, which(all.sex=="female"))
hi.21.female=hi[intersect(n21, which(all.sex=="female"))]
hi.27.male=hi[intersect(n27, which(all.sex=="male"))]; g27.male=intersect(n27, which(all.sex=="male"));g27.female=intersect(n27, which(all.sex=="female"))
hi.27.female=hi[intersect(n27, which(all.sex=="female"))];
hi.28=hi[n28]; g28=n28
col <- colorRampPalette(c("turquoise","royalblue4"))(100)
heatmap(as.matrix(acd.c.an[g27.female,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_A, gen27, female ")
heatmap(as.matrix(acd.c.an[g27.male,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_A, gen27, male ")
heatmap(acd.c.an[g27.male,], Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5, main="Muller_A, gen27, male ")
##step3, take each muller element and calculate alpha, beta, se.alpha, se.beta
bgc=function(muller.data, gen.num) #background data, genotype data with chrom.pos, geno data
{coefi={}
for(i in 1:ncol(muller.data)) #loop through each cluster
	{h=hi[gen.num]
	p=as.numeric(as.character(muller.data[gen.num,i]))
	set.seed(23)
	m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm = "port", control=nls.control(maxiter =100, warnOnly=TRUE),lower=c(a=-3, b=-3), upper=c(a=3,b=3)) #lower=c(a=-1, b=-1), upper=c(a=1,b=1),
	a=coef(m)[1]; b=coef(m)[2]
	se.a=coef(summary(m))[1, "Std. Error"];se.b=coef(summary(m))[2, "Std. Error"];
	coefi=rbind(coefi, c(a, b, se.a, se.b))
	x<-data.frame(h=seq(0, 1, length.out=1000));
	y<-predict(m,x)
	plot(h, p, col=rgb(0,1,1, 0.5), pch=16, cex=2, main=colnames(muller.data)[i], xlim=c(0,1), ylim=c(0,1))
	abline(0,1, lty=2, lwd=1.5)
	lines(x[,1], y, col="blue", lwd=2)
	text(0.1,1, paste("alpha= ", round(a, 3)), cex=0.7)
	text(0.1,0.95, paste("beta= ", round(b, 3)), cex=0.7)
	}
colnames(coefi)=c("a", "b", "se.a", "se.b")
return(coefi)
}
#generation 28 calculate alpha, beta, and their se for each muller element
a.g28=bgc(a.c.an, n28 )
b.g28=bgc(b.c.an, n28 )
dc.g28=bgc(dc.c.an, n28 )
e.g28=bgc(e.c.an, n28 )
f.g28=bgc(f.c.an, n28 )
g28=rbind(a.g28[rev(seq(1:nrow(a.g28))),], dc.g28, b.g28[rev(seq(1:nrow(b.g28))),],  e.g28, f.g28)
#generation 27 calculate alpha, beta, and their se; seperating males and females
a.g27.male=bgc(a.c.an, g27.male )
b.g27.male=bgc(b.c.an, g27.male  )
dc.g27.male=bgc(dc.c.an, g27.male  )
e.g27.male=bgc(e.c.an, g27.male )
f.g27.male=bgc(f.c.an,g27.male  )
g27.male=rbind(a.g27.male[rev(seq(1:nrow(a.g27.male))),], dc.g27.male,b.g27.male[rev(seq(1:nrow(b.g27.male))),], e.g27.male,  f.g27.male)
a.g27.female=bgc(a.c.an, g27.female )
b.g27.female=bgc(b.c.an, g27.female  )
dc.g27.female=bgc(dc.c.an, g27.female  )
e.g27.female=bgc(e.c.an, g27.female )
f.g27.female=bgc(f.c.an,g27.female  )
g27.female=rbind(a.g27.female[rev(seq(1:nrow(a.g27.female))),], dc.g27.female,b.g27.female[rev(seq(1:nrow(b.g27.female))),], e.g27.female,  f.g27.female)

#Generation 27 bgc with all males and females
a.g27=bgc(a.c.an, n27 )
b.g27=bgc(b.c.an, n27 )
dc.g27=bgc(dc.c.an, n27 )
e.g27=bgc(e.c.an, n27 )
f.g27=bgc(f.c.an, n27 )
g27=rbind(a.g27[rev(seq(1:nrow(a.g27))),], dc.g27,b.g27[rev(seq(1:nrow(b.g27))),], e.g27,  f.g27)

#generation 21 calculate alpha, beta, and their se for all individuals in generation 21
a.g21=bgc(a.c.an, n21 )
b.g21=bgc(b.c.an, n21 )
dc.g21=bgc(dc.c.an, n21 )
e.g21=bgc(e.c.an, n21 )
f.g21=bgc(f.c.an, n21 )
g21=rbind(a.g21[rev(seq(1:nrow(a.g21))),], dc.g21, b.g21[rev(seq(1:nrow(b.g21))),],  e.g21,f.g21)
#generation early calculate alpha, beta, and their se
a.ge=bgc(a.c.an, ne )
b.ge=bgc(b.c.an, ne )
dc.ge=bgc(dc.c.an, ne )
e.ge=bgc(e.c.an, ne )
f.ge=bgc(f.c.an, ne )
ge=rbind(a.ge[rev(seq(1:nrow(a.ge))),], dc.ge,b.ge[rev(seq(1:nrow(b.ge))),], e.ge,  f.ge) #be careful, the "ge" is reassigned by a different object

#plot for all muller elements 
#plot alpha for gen 28
muller=c(nrow(a.g28),(nrow(dc.g28)+nrow(a.g28)), nrow(dc.g28)+(nrow(b.g28)+nrow(a.g28)), nrow(dc.g28)+nrow(e.g28)+(nrow(b.g28)+nrow(a.g28)))
muller.element=c(rep("MullerA", nrow(a.g28)),rep("MullerDC", nrow(dc.g28)), rep("MullerB", nrow(b.g28)),rep("MullerE", nrow(e.g28)), rep("MullerF", nrow(f.g28)))
beta.limit=3
x=seq(1:length(g28[,2]));y=g28[,2];se=g28[,4] 
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="darkolivegreen ")
abline(h=0, lty=2)
abline(v=muller)
#plot alpha for gen 28
x=seq(1:length(g28[,1]));y=g28[,1];se=g28[,3]
plot(x,y, ylim=c(-1.1, 1.2), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="dodgerblue4 ")
abline(h=0, lty=2)
abline(v=muller)
##GENEEATION 27
#plot beta for gen 27
x=seq(1:length(g27[,2]));y=g27[,2];se=g27[,4]
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="forestgreen")
abline(h=0, lty=2)
abline(v=muller)
#plot alpha for gen 27
x=seq(1:length(g27[,1]));y=g27[,1];se=g27[,3]
plot(x,y, ylim=c(-1.1, 1.2), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="royalblue1")
abline(h=0, lty=2)
abline(v=muller)
###GENERATION 27 MALES
abline(v=muller)
#plot beta for gen 27, male
x=seq(1:length(g27.male[,2]));y=g27.male[,2];se=g27.male[,4]
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions", main="G27, males")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="forestgreen")
abline(h=0, lty=2)
abline(v=muller)
#plot alpha for gen 27,male
x=seq(1:length(g27.male[,1]));y=g27.male[,1];se=g27.male[,3]
plot(x,y, ylim=c(-1.1, 1.2), pch=16, ylab=expression(alpha), xlab="Relative Positions", main="G27, males")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="royalblue1")
abline(h=0, lty=2)
abline(v=muller)
###GENERATION 27 FEMALES  
#plot beta for gen 27, female
x=seq(1:length(g27.female[,2]));y=g27.female[,2];se=g27.female[,4]
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions", main="G27, females")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="forestgreen")
abline(h=0, lty=2)
abline(v=muller)
#plot alpha for gen 27, female
x=seq(1:length(g27.female[,1]));y=g27.female[,1];se=g27.female[,3]
plot(x,y, ylim=c(-1.1, 1.2), pch=16, ylab=expression(alpha), xlab="Relative Positions", main="G27, females")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="royalblue1")
abline(h=0, lty=2)
abline(v=muller)

###plot beta difference for gen 27, between female and males
x=seq(1:length(g27.female[,2]));y=g27.female[,2]-g27.male[,2];se=sqrt(g27.female[,4]^2+g27.male[,4]^2)
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions", main="G27, females-male")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="aquamarine4")
abline(h=0, lty=2)
abline(v=muller)
rgb.convert=function(col.vect,fr)
{cols={}
for(i in col.vect)
	{cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}
muller.col=c( "cornflowerblue","gold", "turquoise","coral", "lightgreen") #color for muller a, cd, b, e
library(vioplot)
muller.element=factor(muller.element, levels=c("MullerA", "MullerDC", "MullerB","MullerE",  "MullerF" ))
vioplot(y~ muller.element,col=rgb.convert(muller.col, 0.5), border = muller.col,horizontal = F, las = 1, cex=1.2, ylab=paste("Female-male" ,expression(beta)), xlab="", lwd=2)
abline(h=0, lty=2, lwd=2)
stripchart(y~ muller.element, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=0.6, bg=muller.element, col=muller.col)  
summary(aov(y~muller.element))
pairwise.t.test(y, muller.element, p.adj = "fdr")

   
###plot alpha difference for gen 27, between female and males
x=seq(1:length(g27.female[,1]));y=g27.female[,1]-g27.male[,1];se=sqrt(g27.female[,3]^2+g27.male[,3]^2)
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(alpha), xlab="Relative Positions", main="G27, females-male")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="blue4")
abline(h=0, lty=2)
abline(v=muller)

#plot beta for gen 21
x=seq(1:length(g21[,2]));y=g21[,2];se=g21[,4]
plot(x,g21[,2], ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="limegreen")
abline(h=0, lty=2)
abline(v=muller)
#plot alpha for gen 21
x=seq(1:length(g21[,1]));y=g21[,1];se=g21[,3]
plot(x,y, ylim=c(-1.5, 1.5), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="steelblue1")
abline(h=0, lty=2)
abline(v=muller)
#plot beta for gen early
x=seq(1:length(ge[,2]));y=ge[,2];se=ge[,4]
plot(x,y, ylim=c(-beta.limit, beta.limit), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="green")
abline(h=0, lty=2)
abline(v=muller)
#plot alpha for gen early
x=seq(1:length(ge[,1]));y=ge[,1];se=ge[,3]
plot(x,y, ylim=c(-1.1, 1.2), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="lightblue1")
abline(h=0, lty=2)
abline(v=muller)

###correlation between g21 and g27
n =5
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
plot(g21[,2], g27[,2], xlim=c(-1.8,1.5), ylim=c(min(c(g27[,2],g21[,2])), max(c(g27[,2],g21[,2]))), pch=21, bg=col_vector[as.numeric(as.factor(muller.element))+9], xlab="Generation 21 beta",ylab="Generation 27 beta" )
abline(0,1); abline(h=0); abline(v=0)
m=lm(dc.g27[,2]~dc.g21[,2]); abline(m);summary(m)
legend(0.5, -0.4, names(table(muller.element)), pch=21, pt.bg=col_vector[c(10:14)], cex=0.8, pt.cex=1)

plot(g21[,1], g27[,1], xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), pch=21,bg=col_vector[as.numeric(as.factor(muller.element))+9], xlab="Generation 21 alpha",ylab="Generation 27 alpha" )
abline(0,1); abline(h=0); abline(v=0)
m=lm(g27[,1]~g21[,1]); summary(m)
legend(-0.95, 1, names(table(muller.element)), pch=21, pt.bg=col_vector[c(10:14)], cex=0.8, pt.cex=1)

###correlation between g28 and g27 beta
n =5
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#beta
plot( g27[,2],g28[,2], xlim=c(-2,2), ylim=c(-2,2), pch=21, bg=col_vector[as.numeric(as.factor(muller.element))+9], xlab="Generation 27 beta",ylab="Generation 28 beta" )
abline(0,1); abline(h=0); abline(v=0)
abline(h=0.5); abline(v=0.5)
m=lm(dc.g27[,2]~dc.g28[,2]); summary(m)
legend(1, -0.4, names(table(muller.element)), pch=21, pt.bg=col_vector[c(10:14)], cex=0.8, pt.cex=1)
#alpha
plot(g28[,1], g27[,1], xlim=c(-1,1), ylim=c(-1,1), pch=21,bg=col_vector[as.numeric(as.factor(muller.element))+9], xlab="Generation 28 alpha",ylab="Generation 27 alpha" )
abline(0,1); abline(h=0); abline(v=0)
m=lm(g27[,1]~g28[,1]); summary(m)
legend(-0.95, 1, names(table(muller.element)), pch=21, pt.bg=col_vector[c(10:14)], cex=0.8, pt.cex=1)

weight1=1/(1+7+8); weight2=7/(1+7+8); weight3=8/(1+7+8)
outliers=which((g21[,2]*weight1+g27[,2]*weight2+g28[,2]*weight3)>0.5)
# outliers=which((g27[,2]-g27[,4])>0)
plot(g27[,2], g28[,2])
points(g27[outliers,2], g28[outliers,2], pch=16, col="gold")
points(g27[others,2], g28[others,2], pch=16, col="blue")

#run bootstrap
N=100000;  r2.diff.21.27={}; r2.diff.other.21.27={};
r2.diff.21.28={}; r2.diff.other.21.28={}
r2.diff.27.28={}; r2.diff.other.27.28={}
r2.21.bo={};r2.27.bo={};r2.28.bo={}
for(per in 1:N)
	{g.c=sample(outliers, 2, replace=T)
	g.o=sample(others, 2, replace=T)
	r2.28=cor(dd.c.an.28[,g.c[1]], dd.c.an.28[,g.c[2]])^2
	r2.27=cor(dd.c.an.27[,g.c[1]], dd.c.an.27[,g.c[2]])^2
	r2.21=cor(dd.c.an.21[,g.c[1]], dd.c.an.21[,g.c[2]])^2
	r2.diff.21.27=c(r2.diff.21.27, (r2.27-r2.21))
	r2.diff.27.28=c(r2.diff.27.28, (r2.28-r2.27))
	r2.diff.21.28=c(r2.diff.21.28, (r2.28-r2.21))
	#calculate change for genome-wide clusters
	r2.27.other=cor(dd.c.an.27[,g.o[1]], dd.c.an.27[,g.o[2]])^2
	r2.21.other=cor(dd.c.an.21[,g.o[1]], dd.c.an.21[,g.o[2]])^2
	r2.28.other=cor(dd.c.an.28[,g.o[1]], dd.c.an.28[,g.o[2]])^2
	r2.diff.other.21.27=c(r2.diff.other.21.27, (r2.27.other-r2.21.other))
	r2.diff.other.21.28=c(r2.diff.other.21.28, (r2.28.other-r2.21.other))
	r2.diff.other.27.28=c(r2.diff.other.27.28, (r2.28.other-r2.27.other))
	#calculate the difference of LD in Barrier vs other cluster
	r2.21.bo=c(r2.21.bo, r2.21-r2.21.other)
	r2.27.bo=c(r2.27.bo, r2.27-r2.27.other)
	r2.28.bo=c(r2.28.bo, r2.28-r2.28.other)
	}	
hist(r2.diff.21.27, breaks=100, col=rgb(0,1,1, 0.3), xlab="Increased LD among outlier clusters (Gen27-Gen21)")
# text(0.3, 10000, "t = 7.4802, df = 99999", cex=0.7)
# text(0.3, 9000, "p-value = 7.479e-14", cex=0.7)
# text(0.3, 8000, "95% CI: 0.00199084, 0.00340455", cex=0.7)
t.test(r2.21.bo);t.test(r2.27.bo);t.test(r2.28.bo)
t.test(r2.diff.21.27);t.test(r2.diff.21.28);t.test(r2.diff.27.28)
t.test(r2.diff.21.27, r2.diff.21.28);t.test(r2.diff.21.28,r2.diff.27.28); t.test(r2.diff.21.27,r2.diff.27.28)
t.test(r2.diff.21.28, r2.diff.other.21.28, paired=T)
t.test( r2.diff.21.27,r2.diff.other.21.27, paired=T)
t.test(r2.diff.27.28,r2.diff.other.27.28, paired=T)
range(r2.diff.27.28);range(r2.diff.21.28);range(r2.diff.21.27)
range(r2.diff.other.27.28);range(r2.diff.other.21.28);range(r2.diff.other.21.27)
library(vioplot)
par(mar=c(5,6,4,1)+.1)

###Figure 4 @@@
#PART A
vioplot(r2.21.bo,r2.27.bo, r2.28.bo,names = rev(c("B-C:g21", "B-C:g27","B-C:g28")),col = rgb.convert(rev(c("limegreen","forestgreen","darkolivegreen")), 0.5), border = rev(c("limegreen","forestgreen","darkolivegreen ")),horizontal = TRUE, las = 1, cex=1.2) 
#stripchart(data.frame(r2.21.bo,r2.27.bo, r2.28.bo),vertical = F,method = "jitter", add = TRUE, pch = 16, cex=0.6, col =rgb.convert(rev(c("limegreen","forestgreen","darkolivegreen")), 0.02))  
abline(v=0, lty=2)
hist(r2.21.bo, col=col.cov("limegreen", 0.5), breaks=80, border=col.cov("limegreen", 1), main="", xlab=bquote(r^2~"Difference"))
hist(r2.27.bo, xlim=c(-1, 0.7),col=col.cov("forestgreen", 0.3), breaks=100, border=col.cov("forestgreen", 0.8), add=T,lab=bquote(Delta~r^2))
hist(r2.28.bo,  xlim=c(-1, 0.7),col=col.cov("darkolivegreen", 0.5), breaks=100, border=col.cov("darkolivegreen", 1), add=T,lab=bquote(Delta~r^2))
legend("topleft", c("g21", "g27", "g28"), fill=c(col.cov("limegreen", 0.3), col.cov("forestgreen", 0.5), col.cov("darkolivegreen", 0.5)), border=c(col.cov("limegreen", 1), col.cov("forestgreen", 1), col.cov("darkolivegreen", 1)), cex=0.7)
abline(v=0, lty=2)

###Figure 4 @@@
#PART. B
col.cov=function(col, fr){return(rgb(col2rgb(col)[1,1]/255,col2rgb(col)[2,1]/255,col2rgb(col)[3,1]/255,fr)) }
par(mar=c(5,6,4,1)+.1)
vioplot(r2.diff.27.28, r2.diff.21.28, r2.diff.21.27,names =  rev(c("B:g27-g21","B:g28-g21", "B:g28-g27" )),col = rev(c(col.cov("indianred1", 1),col.cov("deepskyblue", 1),col.cov("gold", 0.5))), border = rev(c(col.cov("indianred1", 1),col.cov("deepskyblue", 1),col.cov("gold", 1))),horizontal = TRUE, las = 1, cex=1.2)
abline(v=0, lty=2)
#colMed=rev(c(rgb(0.2,0.8,1,1),rgb(0.5,0.7,0.9,1), rgb(0.5,0.8,0.1,1),rgb(0.5,0.7,0.5,1),rgb(1,0.843137,0,1),rgb(1,0.3,0.1,1))),
hist(r2.diff.27.28, xlim=c(-0.5, 0.5),col=col.cov("gold", 0.5), breaks=80, border=col.cov("gold", 1), main="", xlab=bquote(Delta~r^2))
hist(r2.diff.21.27, xlim=c(-1, 0.7),col=col.cov("indianred1", 0.3), breaks=100, border=col.cov("indianred1", 0.8), add=T,lab=bquote(Delta~r^2))
hist(r2.diff.21.28, xlim=c(-1, 0.7),col=col.cov("deepskyblue", 0.5), breaks=100, border=col.cov("deepskyblue", 1), add=T,lab=bquote(Delta~r^2))
legend(-0.45, 13000, title="Barrier LD Change", c("g27-g21", "g28-g21", "g28-g27"), fill=c(col.cov("indianred1", 0.3), col.cov("deepskyblue", 0.5), col.cov("gold", 0.5)), border=c(col.cov("indianred1", 1), col.cov("deepskyblue", 1), col.cov("gold", 1)), cex=0.7)
abline(v=0, lty=2)


hist(r2.diff.27.28, xlim=c(-1, 0.7),col=rgb(1,0.843137,0,0.5), breaks=80, border=rgb(1,0.843137,0,0.6), main="", xlab=bquote(Delta~r^2))
hist(r2.diff.other.21.28, xlim=c(-1, 0.8),add=T,col=rgb(0.5,0.7,0.5,0.5), breaks=200, border=rgb(0.5,0.7,0.5,0.1))
hist(r2.diff.21.28, xlim=c(-1, 0.7),col=rgb(0.5,0.8,0.1,0.5), breaks=100, border=rgb(0.5,0.7,0.1,0.5), add=T,lab=bquote(Delta~r^2))
hist(r2.diff.other.21.27, xlim=c(-1, 0.8),add=T,col=rgb(0.5,0.7,0.9,0.5), breaks=200, border=rgb(0.5,0.7,0.9,0.1))
hist(r2.diff.21.27, xlim=c(-1, 0.7),col=rgb(0.2,0.8,1,0.5), breaks=100, border=rgb(0.2,0.7,1,0.5), add=T,lab=bquote(Delta~r^2))
hist(r2.diff.other.27.28, xlim=c(-1, 0.8),add=T,col=rgb(1,0.3,0.1,0.5), breaks=100, border=rgb(1,0.3,0,0.1))
hist(r2.diff.27.28, xlim=c(-1, 0.7),col=rgb(1,0.843137,0,0.1), breaks=80, border=rgb(1,0.843137,0,0.2), add=T,lab=bquote(Delta~r^2))
legend(-1, 12000, title="Clusters",  c("Barrier: g27-g21", "Control: g27-g21","Barrier: g28-g21", "Control: g28-g21", "Barrier: g28-g27", "Control: g28-g27"), fill=c(rgb(0.2,0.8,1,0.5),rgb(0.5,0.7,0.9,0.5), rgb(0.5,0.8,0.1,0.5),rgb(0.5,0.7,0.5,0.5),rgb(1,0.843137,0,0.5),rgb(1,0.3,0.1,0.5)), border=c(rgb(0.2,0.8,1,0.5),rgb(0.5,0.7,0.9,0.5), rgb(0.5,0.8,0.1,0.5),rgb(0.5,0.7,0.5,0.5),rgb(1,0.843137,0,0.5),rgb(1,0.3,0.1,0.5)), cex=0.7)

par(mar=c(5,6,4,1)+.1)
vioplot(r2.diff.other.27.28,r2.diff.27.28, r2.diff.other.21.28,r2.diff.21.28,r2.diff.other.21.27, r2.diff.21.27,names =  rev(c("B:g27-g21", "C:g27-g21","B:g28-g21", "C:g28-g21", "B:g28-g27", "C:g28-g27")),col = rev(c(rgb(0.2,0.8,1,0.5),rgb(0.5,0.7,0.9,0.5), rgb(0.5,0.8,0.1,0.5),rgb(0.5,0.7,0.5,0.5),rgb(1,0.843137,0,0.5),rgb(1,0.3,0.1,0.5))), border = "brown",horizontal = TRUE, las = 1, cex=1.2)
#colMed=rev(c(rgb(0.2,0.8,1,1),rgb(0.5,0.7,0.9,1), rgb(0.5,0.8,0.1,1),rgb(0.5,0.7,0.5,1),rgb(1,0.843137,0,1),rgb(1,0.3,0.1,1))),
abline(v=0, lty=2)


#network analysis
###plot cluster in muller B and muller E
dd.c.an.be=data.frame( b.c.an[,rev(seq(1:ncol(b.c.an)))], e.c.an);
dd.c.an.be.28=dd.c.an.be[n28,]; dd.c.an.be.27=dd.c.an.be[n27,]; dd.c.an.be.21=dd.c.an.be[n21,]; dd.c.an.be.e=dd.c.an.be[ne,]; 
cutoff=0.6
m28.be=cor(dd.c.an.be.28); m28.be[lower.tri(m28.be, diag = T)]=0;m28.be[which(m28.be<cutoff)]=0
m27.be=cor(dd.c.an.be.27); m27.be[lower.tri(m27.be, diag = T)]=0;m27.be[which(m27.be<cutoff)]=0; 
m21.be=cor(dd.c.an.be.21); m21.be[lower.tri(m21.be, diag = T)]=0;m21.be[which(m21.be<cutoff)]=0
me.be=cor(dd.c.an.be.e); me.be[lower.tri(me.be, diag = T)]=0;me.be[which(me.be<cutoff)]=0
####
###plot cluster in muller A and CD
dd.c.an.acd=data.frame( a.c.an[,rev(seq(1:ncol(a.c.an)))], dc.c.an);
dd.c.an.acd.28=dd.c.an.acd[n28,]; dd.c.an.acd.27=dd.c.an.acd[n27,]; dd.c.an.acd.21=dd.c.an.acd[n21,]; dd.c.an.acd.e=dd.c.an.acd[ne,]; 
m28.acd=cor(dd.c.an.acd.28); m28.acd[lower.tri(m28.acd, diag = T)]=0;m28.acd[which(m28.acd<cutoff)]=0
m27.acd=cor(dd.c.an.acd.27); m27.acd[lower.tri(m27.acd, diag = T)]=0;m27.acd[which(m27.acd<cutoff)]=0; 
m21.acd=cor(dd.c.an.acd.21);m21.acd[lower.tri(m21.acd, diag = T)]=0; m21.acd[which(m21.acd<cutoff)]=0
me.acd=cor(dd.c.an.acd.e); me.acd[lower.tri(me.acd, diag = T)]=0;me.acd[which(me.acd<cutoff)]=0
###network function
network.plot.muller=function(m28, n1,n2,main, col1, col2)
{g28=graph_from_adjacency_matrix(as.matrix(m28), mode = "undirected", weighted = T, diag = F)
 colors.muller <- c("black",col1, col2)
 ecolor = colorRampPalette(c("black","darkolivegreen1" ))(100-cutoff*10)
plot.igraph(g28, vertex.cex=0.1, edge.color=ecolor[E(g28)$weight*100-cutoff*10],layout=layout_in_circle,edge.width=E(g28)$weight*2,vertex.size=3,vertex.color=colors.muller[c(1,rep(2, n1),rep(3, n2))],vertex.label=NA, margin=0,main=main)#(20*authority.score(g28)$vector)
print(summary(g28))
} 
###plot all clusters in muller E and B
par(mfrow=c(2,2),mar=c(1, 1, 1, 1))
network.plot.muller(me.be,n1=ncol(b.c.an), n2=ncol(e.c.an),  main="Generation 0-5", "turquoise", "coral")
network.plot.muller(m21.be,n1=ncol(b.c.an), n2=ncol(e.c.an),  main="Generation 21", "turquoise", "coral")
network.plot.muller(m27.be, n1=ncol(b.c.an), n2=ncol(e.c.an), main="Generation 27", "turquoise", "coral")
network.plot.muller(m28.be, n1=ncol(b.c.an), n2=ncol(e.c.an), main="Generation 28", "turquoise", "coral")
dev.off()

###plot all clusters in muller A and CD
par(mfrow=c(2,2),mar=c(1, 1, 1, 1))
network.plot.muller(me.acd,n1=ncol(a.c.an), n2=ncol(dc.c.an),  main="Generation 0-5", "cornflowerblue", "gold")
network.plot.muller(m21.acd,n1=ncol(a.c.an), n2=ncol(dc.c.an),  main="Generation 21", "cornflowerblue", "gold")
network.plot.muller(m27.acd, n1=ncol(a.c.an), n2=ncol(dc.c.an), main="Generation 27", "cornflowerblue", "gold")
network.plot.muller(m28.acd, n1=ncol(a.c.an), n2=ncol(dc.c.an), main="Generation 28", "cornflowerblue", "gold")
dev.off()

###Figure 4 @@@@
###.        @@@@
#Barrier clusters
#run LD test
#genetic data 
dd.c.an=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))],  dc.c.an, b.c.an[,rev(seq(1:ncol(b.c.an)))], e.c.an, f.c.an)
dd.c.an.27=dd.c.an[n27, ]
dd.c.an.28=dd.c.an[n28, ]
dd.c.an.21=dd.c.an[n21, ]
dd.c.an.e=dd.c.an[ne, ]
weight1=1/(1+7+8); weight2=7/(1+7+8); weight3=8/(1+7+8)
outliers=which((g21[,2]*weight1+g27[,2]*weight2+g28[,2]*weight3)>0.5)
others=seq(1:nrow(g21))[-outliers]
table(muller.element[outliers])
dd.c.an.27.b=dd.c.an.27[,outliers] #
dd.c.an.21.b=dd.c.an.21[,outliers] #
dd.c.an.28.b=dd.c.an.28[,outliers] #
dd.c.an.e.b=dd.c.an.e[,outliers] #

ld.cutoff=0.7
m28=cor(dd.c.an.28.b);m28[lower.tri(m28, diag = T)]=0; m28[which(m28<ld.cutoff)]=0; 
m27=cor(dd.c.an.27.b); m27[lower.tri(m27, diag = T)]=0;m27[which(m27<ld.cutoff)]=0
m21=cor(dd.c.an.21.b); m21[lower.tri(m21, diag = T)]=0; m21[which(m21<ld.cutoff)]=0
me=cor(dd.c.an.e.b); me[lower.tri(me, diag = T)]=0; me[which(me<ld.cutoff)]=0

#label clusters from different muller element with 
muller.barriers=muller.element[outliers]
colnames.muller=function(m){colnames(m)=paste(muller.element[outliers],outliers, sep="."); rownames(m)=paste(muller.element[outliers],outliers, sep="."); return(m)}
m28=colnames.muller(m28);m21=colnames.muller(m21);m27=colnames.muller(m27);me=colnames.muller(me)

library(igraph); library(network)
#m=data.frame(read.csv("ld.Dprime.May7.2018.csv", header=T, sep=","))
network.plot.muller.barrier=function(m28, ecut, main)
{g28=graph_from_adjacency_matrix(as.matrix(m28), mode = "undirected", weighted = T, diag = F)
 colors.muller <- c("turquoise","gold", "coral")
 ecolor = colorRampPalette(c("black", "darkolivegreen1"))(100-ld.cutoff*10)
plot.igraph(g28, vertex.cex=0.1, edge.color=ecolor[E(g28)$weight*100-ld.cutoff*10],layout=layout_in_circle,edge.width=E(g28)$weight*2,vertex.size=6,vertex.color=colors.muller[as.numeric(as.factor(muller.barriers))],vertex.label=NA, margin=0,main=main)#(20*authority.score(g28)$vector)
print(summary(g28))
degree_distribution(g28, cumulative=T) #edge.width=E(g28)$weight*2
} # layout=layout_in_circle,  #layout_with_fr #layout_on_sphere #add_layout_, component_wise, layout_as_bipartite, layout_as_star, layout_as_tree, layout_in_circle, layout_nicely, layout_on_grid, layout_on_sphere, layout_randomly, layout_with_dh, layout_with_fr, layout_with_gem, layout_with_graphopt, layout_with_kk, layout_with_lgl, layout_with_mds, layout_with_sugiyama, merge_coords, norm_coords, normalize
#edge.color=ecolor[(E(g28)$weight>ecut)+1]
#pdf("Figure.LDnetwork.barrier.pdf")
par(mfrow=c(2,2),mar=c(1, 1, 1, 1) )
network.plot.muller.barrier(me, 0.9, main="Generation 0-5")
network.plot.muller.barrier(m21, 0.9, main="Generation 21")
network.plot.muller.barrier(m27, ecut=0.9, main="Generation 27")
network.plot.muller.barrier(m28, ecut=0.9, main="Generation 28")
#dev.off()
#find modules and then plot it
g28=graph_from_adjacency_matrix(as.matrix(m28), mode = "undirected", weighted = T, diag = F)
  fc.28 <- fastgreedy.community(g28)
  colors <- rainbow(max(membership(fc.28)))
  plot(g28,vertex.color=colors[membership(fc.28)], layout=layout.fruchterman.reingold)

   #here's how to call a group, and individual vertex within a group 
   #g6=unlist(strsplit(fc[[6]], split=" "))
   #v=g6[which(g6=="5060967")]
   #write_graph(g.copy, "graphLDlessthan0.5.txt", "edgelist")


#some zoom-in plots
#for Muller DC only
plot(dc.g21[,2], dc.g27[,2], xlim=c(-1,1), ylim=c(-1,1))
abline(0,1); abline(h=0); abline(v=0)
m=lm(dc.g27[,2]~dc.g21[,2]); summary(m)
plot(dc.g21[,1], dc.g27[,1], xlim=c(-1,1), ylim=c(-1,1))
abline(0,1); abline(h=0); abline(v=0)
m=lm(dc.g27[,1]~dc.g21[,1]); summary(m)

#plot beta for gen 27
x=seq(1:length(dc.g27[,2]));y=dc.g27[,2];se=dc.g27[,4]
plot(x,y, ylim=c(-2, 1.5), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="forestgreen")
abline(h=0, lty=2)
#plot alpha for gen 27
x=seq(1:length(dc.g27[,1]));y=dc.g27[,1];se=dc.g27[,3]
plot(x,y, ylim=c(-1, 0.8), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="royalblue1")
abline(h=0, lty=2)

#plot beta for gen 21
x=seq(1:length(dc.g21[,2]));y=dc.g21[,2];se=dc.g21[,4]
plot(x,dc.g21[,2], ylim=c(-2, 2.3), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="limegreen")
abline(h=0, lty=2)
#plot alpha for gen 21
x=seq(1:length(dc.g21[,1]));y=dc.g21[,1];se=dc.g21[,3]
plot(x,y, ylim=c(-1.5, 0.8), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="steelblue1")
abline(h=0, lty=2)

#plot beta for gen early
x=seq(1:length(dc.ge[,2]));y=dc.ge[,2];se=dc.ge[,4]
plot(x,y, ylim=c(-2.6, 0.6), pch=16, ylab=expression(beta), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="green")
abline(h=0, lty=2)
#plot alpha for gen early
x=seq(1:length(dc.ge[,1]));y=dc.ge[,1];se=dc.ge[,3]
plot(x,y, ylim=c(-1, 0.8), pch=16, ylab=expression(alpha), xlab="Relative Positions")
arrows(x, y-se, x, y+se, length=0.05, angle=90, code=3, col="lightblue1")
abline(h=0, lty=2)



