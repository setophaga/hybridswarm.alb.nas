

# install.packages("https://cran.r-project.org/bin/macosx/contrib/4.0/fpc_2.2-5.tgz",method="libcurl")
# install.packages("~/Downloads/mclust_5.4.6.tar")

#setwd("~/Desktop/hybridswarm/1.pipeline/posterior.FM.3p.July20/")
library(readxl)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
#library(fpc)
library("RColorBrewer")
library(qlcMatrix)
refd=read.csv("../alb03.male.female.nas00.csv")#on laptop: ~/Desktop/hybridswarm/1.pipeline/parentalstrain
#STEP 1 get ancestry-specific genotype from AncestryHMM output
files <- list.files(path = "/scratch/silu/abo.nas/alb03.nas00.ref/posterior.all.alb03.nas00", pattern = "*.posterior", full.names = T)
#files <- list.files(path = "~/Desktop/hybridswarm/1.pipeline/posterior.alb03.nas00.all.July21th", pattern = "*.posterior", full.names = T)
#use this above line if run on laptop
e.sp=read.csv(files[1], sep="\t")
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
v[which(anm[,1]==1)]=0 #albomicans
v[which(anm[,2]==1)]=0.5 #heterozygotes
v[which(anm[,3]==1)]=1 #nasulat
print(maxresult$max[which(maxresult$max<posterior.cutoff)])
if(length(which(maxresult$max<posterior.cutoff))>0)
	{v[which(maxresult$max<posterior.cutoff)]=NA }
if(sum(is.na(maxresult$max))>0)
	{v[is.na(maxresult$max)]=NA }
# for(i in 1:nrow(an))
	# {if(max(an[i,],na.rm=T)>posterior.cutoff)
		# {v.g=(3-which.max(an[i,]))/2}
	# else{v.g=NA}
	# v=c(v, v.g)
	# }
spd=cbind(spd, v)
}
dat.sp=data.frame(e.sp[,1:2], spd)
dat.sp=dat.sp[order(dat.sp$chrom, dat.sp$position),]
table(dat.sp$chrom)
spd=dat.sp[,3:148]
spd.het=spd
for(j in 1:ncol(spd))
{spd.het[which(spd[,j]==1),j]=0;spd.het[which(spd[,j]==0.5),j]=1}
 spd.dc=spd[which(dat.sp$chrom=="Muller_DC"),];spd.het.dc=spd.het[which(dat.sp$chrom=="Muller_DC"),]
 spd.a=spd[which(dat.sp$chrom=="Muller_A"),];spd.het.a=spd.het[which(dat.sp$chrom=="Muller_A"),]
 spd.b=spd[which(dat.sp$chrom=="Muller_B"),];spd.het.b=spd.het[which(dat.sp$chrom=="Muller_B"),]
 spd.e=spd[which(dat.sp$chrom=="Muller_E"),];spd.het.e=spd.het[which(dat.sp$chrom=="Muller_E"),]
 spd.f=spd[which(dat.sp$chrom=="Muller_F"),];spd.het.f=spd.het[which(dat.sp$chrom=="Muller_F"),]
 
 d=read.csv("nasutaXabomincans.all.simple.csv") #match generation information for each individual
 #d=read.csv("~/Desktop/hybridswarm/nasutaXabomincans.all.simple.csv") #run this if on laptop
 mul.map=read.csv("old.lib.mullers.mapreads.csv")
 #mul.map=read.csv("~/Desktop/hybridswarm/old.lib.mullers.mapreads.csv") #tun this if on laptop
 d$Gender[which(d$Gender=="F")]="female"; d$Gender=factor(d$Gender, levels=c("female", "male", "ND"));table(d$Gender)
 d28=d[which(d$Generation==28),];table(d$Gender,d$Generation)
bg={}; mm={}
for( i in 1:length(names))
	{row=which(names[i]==d$prefix)
	mrow=which(names[i]==mul.map$id)
	bg=rbind(bg, d[row,])
	mm=rbind(mm, mul.map[mrow,])
	}
bgd=data.frame(names, bg, mm)
length((colMeans(spd, na.rm=T)))
bgd$HI=(colMeans(spd, na.rm=T))
bgd$het=(colMeans(spd.het, na.rm=T))
d$Gendersim=as.character(d$Gender); d$Gendersim[which(d$Gendersim=="female")]="F"; d$Gendersim[which(d$Gendersim=="male")]="M"; 
bgd$sex.g[which(is.na(as.character(bgd$sex.g)))]=d$Gendersim[which(is.na(as.character(bgd$sex.g)))]
bgd$hi.a=(colMeans(spd.a, na.rm=T))
bgd$hi.dc=(colMeans(spd.dc, na.rm=T))
bgd$hi.e=(colMeans(spd.e, na.rm=T))
bgd$hi.b=(colMeans(spd.b, na.rm=T))
bgd$hi.f=(colMeans(spd.f, na.rm=T))
bgd$het.a=(colMeans(spd.het.a, na.rm=T))
bgd$het.dc=(colMeans(spd.het.dc, na.rm=T))
bgd$het.e=(colMeans(spd.het.e, na.rm=T))
bgd$het.b=(colMeans(spd.het.b, na.rm=T))
bgd$het.f=(colMeans(spd.het.f, na.rm=T))
#correct for sex mis-identification
prob=intersect(which(bgd$sex.auto.fr>0.7),which(bgd$sex.auto.fr<0.9)) 
ok=c(intersect(which(bgd$sex.g[prob]=="F"), which(bgd$Gender[prob]=="male")), intersect(which(bgd$sex.g[prob]=="M"), which(bgd$Gender[prob]=="male")))
prob=prob[-ok]
bgd$sex.g[prob]=NA
nnn=paste( bgd$sex.g,bgd$Generation,bgd$Species, names,sep=".")
colnames(spd)=nnn
write.csv(spd, "alb03Xnas00.hybrids.haplotype.csv")
write.csv(bgd, "alb03Xnas00.hybrids.backgroundinfo.csv")

sn=sort(nnn) #sort by generation id
#make string with row numbers ordered by chromosomes in the right order with respect to the fusions
sf=sn[1:65]; sm=sn[66:140]
sf.plot=sf[c(59:65,56:58, 1:55)]
sm.plot=sm[c(67:69,70:75, 49:66,1:48 )]

spd.f=spd[,rev(sf.plot)]
spd.m=spd[,rev(sm.plot)]
col.sp <- colorRampPalette(c("royalblue4","turquoise"))(100)
jpeg(filename="female.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[,rev(sf.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()
jpeg(filename="female.muller.cd.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd.f[which(dat.sp$chrom=="Muller_DC"),]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()


jpeg(filename="male.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd[,rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()
jpeg(filename="male.muller.cd.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd.m[which(dat.sp$chrom=="Muller_DC"),]) , Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3)
dev.off()

spd.dc.f=spd.f[which(dat.sp$chrom=="Muller_DC"),]
spd.dc.m=spd.m[which(dat.sp$chrom=="Muller_DC"),]
spd.dc=spd[which(dat.sp$chrom=="Muller_DC"),]

#STEP 2 get neoSexChrom-specific genotype from AncestryHMM output
files.neosex <- list.files(path = "/scratch/silu/abo.nas/alb03.nas00.ref/posterior.albFM.nas.2p", pattern = "*.posterior", full.names = T)
#files <- list.files(path = "~/Desktop/hybridswarm/1.pipeline/posterior.alb03.nas00.all.July21th", pattern = "*.posterior", full.names = T)
#use this above line if run on laptop
e.neosex=read.csv(files.neosex[1], sep="\t")
posterior.cutoff=0.9 #only take the ancestry genotype call when posterior probability is greater than the cutoff
names={}; spd.neosex={}; h={}; 
for(ind in files.neosex)
{name=strsplit(strsplit(ind,"/")[[1]][7], ".posterior")[[1]][1] #use [8] IF run on laptop
names=c(names, name)
 d=read.csv(ind, sep="\t")
v=rep(NA,nrow(d))
an=Matrix(as.matrix(d[,3:5]))
#if the posteriors are less than cutoff, don't consider them
maxresult=rowMax(an, which=T)
anm=maxresult$which
v[which(anm[,1]==1)]=0
v[which(anm[,2]==1)]=0.5
v[which(anm[,3]==1)]=1
print(maxresult$max[which(maxresult$max<posterior.cutoff)])
if(length(which(maxresult$max<posterior.cutoff))>0)
	{v[which(maxresult$max<posterior.cutoff)]=NA }
if(sum(is.na(maxresult$max))>0)
	{v[is.na(maxresult$max)]=NA }
# for(i in 1:nrow(an))
	# {if(max(an[i,],na.rm=T)>posterior.cutoff)
		# {v.g=(3-which.max(an[i,]))/2}
	# else{v.g=NA}
	# v=c(v, v.g)
	# }
spd.neosex=cbind(spd.neosex, v)
}
dat.neosex=data.frame(e.neosex[,1:2], spd.neosex)
dat.neosex=dat.neosex[order(dat.neosex$chrom, dat.neosex$position),]
table(dat.neosex$chrom)
spd.neosex=dat.neosex[,3:148]
 
nnn=paste( bgd$sex.g,bgd$Generation,bgd$Species,names, sep=".")
colnames(spd.neosex)=nnn

#write.csv(spd.neosex, "alb03Xnas00.hybrids.neosex.haplotype.csv")

spd.f.neosex=spd.neosex[,rev(sf.plot)]
spd.m.neosex=spd.neosex[,rev(sm.plot)]
col.neosex <- colorRampPalette(c("gold","forestgreen"))(100)
jpeg(filename="female.neosex.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd.f.neosex) , Colv = NA, Rowv = NA, scale="none", col=col.neosex, cexRow=0.3)
dev.off()
jpeg(filename="female.muller.cd.neosex.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd.f.neosex[which(dat.neosex$chrom=="Muller_DC"),]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex, cexRow=0.3)
dev.off()
#males
jpeg(filename="male.neosex.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd.m.neosex) , Colv = NA, Rowv = NA, scale="none", col=col.neosex, cexRow=0.3)
dev.off()
jpeg(filename="male.muller.cd.neosex.heatmap.jpeg",width=8,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(spd.m.neosex[which(dat.neosex$chrom=="Muller_DC"),]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex, cexRow=0.3)
dev.off()

spd.neosex.dc.f=spd.neosex.f[which(dat.neosex$chrom=="Muller_DC"),]
spd.neosex.dc.m=spd.neosex.m[which(dat.neosex$chrom=="Muller_DC"),]
spd.neosex.dc=spd.neosex[which(dat.neosex$chrom=="Muller_DC"),]
#check if the pposition on chrom MUllder CD are included in the neosex haplotyping vs species haplotyping
sum(dat.neosex$position[dat.neosex$chrom=="Muller_DC"] %in% dat.sp$position[dat.sp$chrom=="Muller_DC"])
sum(dat.sp$position[dat.sp$chrom=="Muller_DC"]  %in% dat.neosex$position[dat.neosex$chrom=="Muller_DC"])

ref=read.csv("../alb03.male.female.2p.csv")
ref.cd=ref[which(ref$CHROM=="Muller_DC"),]
dat.neosex.cd=dat.neosex[which(dat.neosex$chrom=="Muller_DC"),]
dat.sp.cd=dat.sp[which(dat.sp$chrom=="Muller_DC"),]
length(dat.neosex.cd$position ); length(ref.cd$POS)
sum(dat.neosex.cd$position %in% ref.cd$POS)==nrow(dat.neosex.cd)
ref.cd.hap=ref.cd[(which(ref.cd$POS %in% dat.neosex.cd$position)),]
#save.image("species.neosexchrom.haplotype.RData")
spd.neosex.dc=spd.neosex[which(dat.neosex$chrom=="Muller_DC"),]
spd.sp.dc=spd[which(dat.sp$chrom=="Muller_DC"),]

##$$$$$$$plot ancestry for each individual
plot.indv=function(indv)
{indv.an=spd.sp.dc[,indv] ; indv.neosex=spd.neosex.dc[,indv]
jpeg(filename=paste("ancestry.neosex.",indv,".jpeg",sep=""),width=10,height=4,units="in",res=500)
plot(dat.sp.cd$position, indv.an, col="blue",ylim=c(0,1))
lines(dat.sp.cd$position, indv.an, col="blue")
points(dat.neosex.cd$position, indv.neosex, col=rgb(0, 1, 0, 0.2),pch=16,cex=0.5)
lines(dat.neosex.cd$position, indv.neosex, col=rgb(0, 1, 0, 0.6))
dev.off()}
#indv="F.28.A1.DBZW33_163_N711_S503_S59_L008_R1_001" #column name from nnn
#"F.28.A1.DBZW33_163_N711_S503_S59_L008_R1_001"                  
# [36] "F.28.A1.DBZW33_164_N711_S504_S83_L008_R1_001" 
sp.ancestry.block.indv=function(indv.an)
{ #ancestry at muller cd
#names(indv.an)=dat.sp.cd$position
indv.an.switches=abs(sign(indv.an-c(indv.an[-1],indv.an[1])))
pos.switch=dat.sp.cd$position[which(indv.an.switches==1)];lastswitch=tail(pos.switch,1)
if(length(pos.switch)==0)
{anblock=data.frame(ancestry=indv.an[1], begin=dat.sp.cd$position[1], end=tail(dat.sp.cd$position, 1))}else{ancestry=indv.an[1]; begin=dat.sp.cd$position[1]; end=pos.switch[1];
	for(p in pos.switch)
		{if(p==lastswitch)
			{if(length(which(dat.sp.cd$position>lastswitch))>0)
				{ancestry=c(ancestry, indv.an[which(dat.sp.cd$position==p)+1])#start from the second block
				begin=c(begin,dat.sp.cd$position[which(dat.sp.cd$position==p)+1])
				end=c(end, tail(dat.sp.cd$position, 1))}
			else{break}}
		else #make a data.frame containing begin and end position of ancestry blocks
			{ancestry=c(ancestry, indv.an[which(dat.sp.cd$position==p)+1])#start from the second block
			begin=c(begin,dat.sp.cd$position[which(dat.sp.cd$position==p)+1])
			end=c(end, (pos.switch[which(pos.switch==p)+1]+1))}
		};anblock=data.frame(ancestry, begin, end)}	
return(anblock)
}
#&&&&&&&&&&&&&&& Combine species haplotyping with ancestry haplotyping
#Re-assign ancestry scores accounting for both species and neosex chrom
#&&&&&&&&&&&&&&& #&&&&&&&&&&&&&&& 
update.ancestry=function(indv)
{#&&&&&&&&&&&&&&& Homozygotes (nas, nas) 
anblock=sp.ancestry.block.indv(indv.an=spd.sp.dc[,indv])
#change indv.neosex within (nas, nas) ranges into 0
indv.neosex=spd.neosex.dc[,indv] 
nas.neosex.rows={};indv.neosex.new=rep(NA, length(indv.neosex))
nas.block=anblock[which(anblock$ancestry==1),] #the ancestry blocks with nasuta
if(nrow(nas.block)>0)
{for(i in 1:nrow(nas.block))
	{if(length(intersect(which(dat.neosex.cd$position > nas.block$begin[i]), which(dat.neosex.cd$position < nas.block$end[i])))>0){nas.neosex.rows=c(nas.neosex.rows,intersect(which(dat.neosex.cd$position > nas.block$begin[i]), which(dat.neosex.cd$position < nas.block$end[i])))}}}
if(length(nas.neosex.rows)>0)
{indv.neosex.new[nas.neosex.rows]=0} #assign zero for all sites within nasuta block
#&&&&&&&&&&&&&&& Homozygous (alb.alb)
#change indv.neosex within (alb, alb) ranges into 1 (if neosex=0), into 0.75 (if neosex =0.5), into -1 (if neosex=1, impossible!!)
alb.block=anblock[which(anblock$ancestry==0),]
alb.neosex.rows={}
if(nrow(alb.block)>0)
{for(i in 1:nrow(alb.block))
	{if(length(intersect(which(dat.neosex.cd$position > alb.block$begin[i]), which(dat.neosex.cd$position < alb.block$end[i])))>0){alb.neosex.rows=c(alb.neosex.rows,intersect(which(dat.neosex.cd$position > alb.block$begin[i]), which(dat.neosex.cd$position < alb.block$end[i])))}}}
if(length(intersect(alb.neosex.rows, which(indv.neosex==0)))>0)
	{indv.neosex.new[intersect(alb.neosex.rows, which(indv.neosex==0))]=1} #assign zero for sites within block with neoXs
if(length(intersect(alb.neosex.rows, which(indv.neosex==0.5)))>0)
	{indv.neosex.new[intersect(alb.neosex.rows, which(indv.neosex==0.5))]=0.75} #assign 0.5 sites within block with alb.neoX, alb.neoYs
if(length(intersect(alb.neosex.rows, which(indv.neosex==1)))>0)
	{indv.neosex.new[intersect(alb.neosex.rows, which(indv.neosex==1))]=-1} #assign -1 sites within block with both neoYs
#all.ancestry.correctedsites=Reduce(union, list(nas.neosex.rows, het.neosex.rows, alb.neosex.rows))
#indv.neosex.new[-all.ancestry.correctedsites]=NA #the sites that are not ancestry matched
#&&&&&&&&&&&&&&& #&&&&&&&&&&&&&&& 
#&&&&&&&&&&&&&&& Heterozygotes: within (nas, alb)--more complicated
#change indv.neosex within (nas, alb) ranges into 0.25 (if neosex=1), into 0.5 (if neosex =0)
het.block=anblock[which(anblock$ancestry==0.5),]
het.neosex.rows={}
for(i in 1:nrow(het.block))
	{if(length(intersect(which(dat.neosex.cd$position > het.block$begin[i]), which(dat.neosex.cd$position<het.block$end[i])))>0){het.neosex.rows=c(het.neosex.rows,intersect(which(dat.neosex.cd$position > het.block$begin[i]), which(dat.neosex.cd$position<het.block$end[i])))}}
if(length(intersect(het.neosex.rows, which(indv.neosex==1)))>0)
{indv.neosex.new[intersect(het.neosex.rows, which(indv.neosex==1))]=0.25} #assign zero for sites within heterozygous block with neoXs
if(length(intersect(het.neosex.rows, which(indv.neosex==0)))>0)
{indv.neosex.new[intersect(het.neosex.rows, which(indv.neosex==0))]=0.5} #assign zero for sites within heterozygous block 
#idnv.neosex==0.5 sites are the (neox, neoy)-- should be either (neoX.alb, mullercd.nas)=0.5 or (neoY.alb, mullercd.nas)=0.25
### use the source ref data to figure out
het.het=intersect(het.neosex.rows, which(indv.neosex==0.5)) #sites that are called (alb.nas) and (neox, neoy)
refd.cd=refd[which(refd$CHROM=="Muller_DC"),]
het.het.pos=dat.neosex.cd$position[het.het]
rows=which(refd.cd$POS %in% het.het.pos); sum(refd.cd$POS[rows] %in% het.het.pos)==length(rows)
nas.anc=(refd.cd$A1.freq.nas-refd.cd$A2.freq.nas); neoy.anc=(refd.cd$A1.freq.m-refd.cd$A2.freq.m); neox.anc=(refd.cd$A1.freq.f-refd.cd$A2.freq.f); 
nas.a1.rows=which(nas.anc>0);nas.a2.rows=which(nas.anc<0); nas.na.rows=which(nas.anc==0)
neox.a1.rows=intersect(which(neox.anc>0),which(neox.anc>neoy.anc));neox.a2.rows=intersect(which(neox.anc<0), which(neox.anc<neoy.anc));neox.na.rows=which(neox.anc==0);
neoy.a1.rows=intersect(which(neoy.anc>0),which(neoy.anc>neox.anc));neoy.a2.rows=intersect(which(neoy.anc<0), which(neoy.anc<neox.anc));neoy.na.rows=which(neoy.anc==0);
length(intersect(neox.a1.rows,neoy.a1.rows)) #check neox.a1 be mutually exclusive to neoy.a1
neox.nas.rows=union(intersect(intersect(neox.a1.rows, nas.a1.rows), rows), intersect(intersect(neox.a2.rows, nas.a2.rows), rows)) 
neoy.nas.rows=union(intersect(intersect(neoy.a1.rows, nas.a1.rows), rows), intersect(intersect(neoy.a2.rows, nas.a2.rows), rows))
length(intersect(neoy.nas.rows, neox.nas.rows))#check this should be zero
neoy.nas.pos=refd.cd$POS[neoy.nas.rows]#nas=neoY, so (neoX.alb, mullercd.nas)=0.5
neox.nas.pos=refd.cd$POS[neox.nas.rows]#nas=neoX, so (neoY.alb, mullercd.nas)=0.25
names(indv.neosex)=dat.neosex.cd$position; 
# sum(neoy.nas.pos %in% names(indv.neosex))==length(neoy.nas.pos)
if(length(neoy.nas.pos)>0)
{indv.neosex.new[as.character(neoy.nas.pos)]=0.5} #nas=neoy, so (neox, nasCD), assign 0.5 for all sites within het block with neo
if(length(neox.nas.pos)>0)
{indv.neosex.new[as.character(neox.nas.pos)]=0.25}#nas=neox, so (neoY, nasCD), assign 0.25 for all sites within het block 
return(indv.neosex.new)
}
###########&&&&&&&&&&&&&&& END OF HAPLOTYPING UPDATES *update.ancestry* function#&&&&&&&&&&&&&&& 
updated.sp.neosex={}
for(indv in nnn)
{updated.anc=update.ancestry(indv)
print(indv)
print(levels(as.factor(updated.anc)))
updated.sp.neosex=cbind(updated.sp.neosex, updated.anc)}
colnames(updated.sp.neosex)=nnn

plot.5hap=updated.sp.neosex
head(plot.5hap)
plot.5hap=updated.sp.neosex
 for(j in 1:ncol(plot.5hap))
 {if(length(which(plot.5hap[,j]==(-1)))>0)
	 {print(colnames(plot.5hap)[j])
	 print(length(plot.5hap[which(plot.5hap[,j]==(-1)),j]))
	  plot.5hap[which(plot.5hap[,j]==(-1)),j]=NA}
 }
 colnames(plot.5hap)=nnn
#count number of recombinations
rec={}; info.sites={};mode={};mode.freq={}; s0={};s0.25={};s0.5={};s0.75={};s1={};l0={};l0.25={};l0.5={};l0.75={};l1={};
for(j in 1:length(nnn))
{recd=sp.ancestry.block.indv(plot.5hap[,nnn[j]])
md=Mode(plot.5hap[,nnn[j]], na.rm=T)
mode=c(mode, md[1])
mode.freq=c(mode.freq, attributes(md)$freq)
rec=c(rec, nrow(recd))
s0=c(s0, length(which(recd$ancestry==0)))
l0=c(l0,sum(recd$end[which(recd$ancestry==0)]-recd$begin[which(recd$ancestry==0)]))
s0.25=c(s0.25, length(which(recd$ancestry==0.25)))
l0.25=c(l0.25,sum(recd$end[which(recd$ancestry==0.25)]-recd$begin[which(recd$ancestry==0.25)]))
s0.5=c(s0.5, length(which(recd$ancestry==0.5)))
l0.5=c(l0.5,sum(recd$end[which(recd$ancestry==0.5)]-recd$begin[which(recd$ancestry==0.5)]))
s0.75=c(s0.75, length(which(recd$ancestry==0.75)))
l0.75=c(l0.75,sum(recd$end[which(recd$ancestry==0.75)]-recd$begin[which(recd$ancestry==0.75)]))
s1=c(s1, length(which(recd$ancestry==1)))
l1=c(l1,sum(recd$end[which(recd$ancestry==1)]-recd$begin[which(recd$ancestry==1)]))
info.sites=c(info.sites,length(na.omit(plot.5hap[,j])))
}
recomb=data.frame(nnn, bgd$name,bgd$HI, bgd$het,rec, info.sites, bgd$Generation, bgd$sex.g ,mode,mode.freq, s0, s0.25, s0.5, s0.75, s1,l0, l0.25, l0.5, l0.75, l1)
colnames(recomb)=c("nnn", "name", "HI", "het","rec", "info.sites", "generation", "sex", "mode","mode.freq","s0", "s0.25", "s0.5", "s0.75", "s1","l0", "l0.25", "l0.5", "l0.75", "l1")
recomb.d=data.frame(recomb, bgd)
write.csv(recomb.d, "hi.recomb.csv")
jpeg(filename="recomb.hi.jpeg",width=12,height=8,units="in",res=500)
plot(recomb$HI, -log(recomb$rec/recomb$info.sites, base=10))
dev.off()


#plot Muller CD haplotype with individuals that have >30% of the haplotypes
all.plot=nnn[-which(recomb.d$info.sites<max(recomb.d$info.sites, na.rm=T)*0.5)]
splot=sort(all.plot)
splot.f=splot[c(41:47,39:40,1:38)]; splot.m=splot[c(83:90,70:82,48:69)]
col.neosex.sp <- colorRampPalette(c("aquamarine2","darkolivegreen3", "gold", "deepskyblue3", "darkslateblue"))(5)
jpeg(filename="mullercd.5hap.female.0.5moresites.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(plot.5hap[,rev(splot.f)]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()

jpeg(filename="mullercd.5hap.male.0.5moresites.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(plot.5hap[,rev(splot.m)]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()

sn=sort(nnn) #sort by generation id
for(j in sf.plot)
{print(j)
print(levels(as.factor(plot.5hap[,j])))}

#make string with row numbers ordered by chromosomes in the right order with respect to the fusions
#0=(CD_nas, CD_nas), 0.25=(NeoY_alb, CD_nas),0.5=(NeoX_alb, CD_nas), 0.75=(NeoX_alb, NeoY_alb), 1=(NeoX_alb, NeoX_alb)
col.neosex.sp <- colorRampPalette(c("aquamarine2","darkolivegreen3", "gold", "deepskyblue3", "darkslateblue"))(5)
jpeg(filename="mullercd.5hap.female.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(plot.5hap[,rev(sf.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()

jpeg(filename="mullercd.5hap.male.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(plot.5hap[,rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()

#plot(x=NULL, y=NULL)
#legend("topright", c("0=(CD_nas, CD_nas)", "0.25=(NeoY_alb, CD_nas)", "0.5=(NeoX_alb, CD_nas)", "0.75=(NeoX_alb, NeoY_alb)", "1=(NeoX_alb, NeoX_alb)"), fill=c("aquamarine2","darkolivegreen3", "gold", "deepskyblue3", "darkslateblue"))

jpeg(filename="mullercd.5hap.male.heatmap.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(plot.5hap[,rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()












############### For future, (if imputation is needed)
##Impute
impute.chr=function(pdc, dt, bin) #pdc=positions for each muller element; dt=genotype matrix; updated.sp.neosex bib=2000000 
{bins=round((pdc[length(pdc)]-pdc[1])/bin)+1
start=pdc[1]
#jpeg(filename="imputation.jpeg",width=12,height=8,units="in",res=500)
#plot(pdc, dt[,27])
for(b in 1:bins)
	{end=start+bin
	cls=intersect(which(pdc>(start-1)),which(pdc<end))
	#abline(v=end, col="red")
	for(j in 1:ncol(dt))
		{sub=mode(dt[cls,j]) #if use mean(), need to add "na.rm=T"
		nas=which(is.na(dt[cls,j]))
		if(length(nas)>0)
			{dt[cls,j][nas]=sub}
		}
	start=start+bin
	#dev.off()
	} ;return(dt)
}
imputed.sp.neosex.block=impute.chr(dat.neosex.cd$position, updated.sp.neosex, bin=2000000)
imputed.sp.neosex=apply(imputed.sp.neosex.block, c(1,2), as.numeric)
#levels(as.factor(imputed.sp.neosex.block))
col.neosex.sp <- colorRampPalette(c("grey","aquamarine2","darkolivegreen3", "gold", "deepskyblue3", "darkslateblue"))(6)
jpeg(filename="mullercd.5hap.female.heatmap.imputed.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t((imputed.sp.neosex[,rev(sf.plot)])), Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()

jpeg(filename="mullercd.5hap.male.heatmap.imputed.jpeg",width=12,height=8,units="in",res=500)
mar=c(1,1,1,1)
heatmap(t(imputed.sp.neosex[,rev(sm.plot)]) , Colv = NA, Rowv = NA, scale="none", col=col.neosex.sp, cexRow=0.3)
dev.off()




