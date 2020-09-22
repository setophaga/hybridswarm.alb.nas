load("x1.haplotypes.sp.new.old.combo.RData")
install.packages("sommer")
library(sommer)
#STEP 1.2 filtration: only keeping individulas with > cutoff percent of sites inferred
#take the individuals with > 80% sites
cutoff=0.8
cnsite=cutoff*nrow(spd)
list={}
for(i in 1:ncol(spd))
{nsite=length(na.omit(spd[,i]));
		if(nsite>cnsite)
		{print(nsite/nrow(spd))
		#print(na.omit(dd[,i]))
			list=c(list, i)}
}

#heatmap(t(as.matrix(dd[,list])), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
#take the subset of individuals with >50% coverage
#remove DBKCL_N701-S507_S17_L002 and DBKCL_N701-S508_S18_L002 due to excessive recomb #rm=which(bgdt$mullercd.rec>30)
list=list[-c(6,7)] #becareful that you only run this once
dt=t(spd[,list]) #genotype matrix
bgdt=bgd[list,]#background info
dt.mb=dt[,which(dat.sp$chrom=="Muller_B")]
dt.mdc=dt[,which(dat.sp$chrom=="Muller_DC")]
dt.me=dt[,which(dat.sp$chrom=="Muller_E")]
dt.mf=dt[,which(dat.sp$chrom=="Muller_F")]
dt.ma=dt[,which(dat.sp$chrom=="Muller_A")]
mullers=c(rep("Muller_A", ncol(dt.ma)),rep("Muller_DC", ncol(dt.mdc)), rep("Muller_B", ncol(dt.mb)),rep("Muller_E", ncol(dt.me))) #rep("Muller_F", ncol(dt.mf)) 

dt=t(spd[,list]) #genotype matrix
bgdt=bgd[list,]#background info
dt.mb=dt[,which(dat$chrom=="Muller_B")]
dt.mdc=dt[,which(dat$chrom=="Muller_DC")]
dt.me=dt[,which(dat$chrom=="Muller_E")]
dt.mf=dt[,which(dat$chrom=="Muller_F")]
dt.ma=dt[,which(dat$chrom=="Muller_A")]

mdc=dat.sp$pos[which(dat.sp$chrom=="Muller_DC")]
ma=dat.sp$pos[which(dat.sp$chrom=="Muller_A")]
mb=dat.sp$pos[which(dat.sp$chrom=="Muller_B")]
me=dat.sp$pos[which(dat.sp$chrom=="Muller_E")]
mf=dat.sp$pos[which(dat.sp$chrom=="Muller_F")]

#step 2 impute chromosome
impute.chr=function(mdc, dt, bin)
{
dt.mdc.27=dt
bins=round((mdc[length(mdc)]-mdc[1])/bin)+1
start=mdc[1]
#jpeg("check.imput.eg.pdf")
#plot(mdc, dt.mdc.27[67,])
for(b in 1:bins)
	{
	end=start+bin
	cls=intersect(which(mdc>(start-1)),which(mdc<end))
	#abline(v=end, col="red")
	for(i in 1:nrow(dt.mdc.27))
		{sub=mean(dt.mdc.27[i,cls], na.rm=T)
		nas=which(is.na(dt.mdc.27[i,cls]))
		if(length(nas)>0)
			{dt.mdc.27[i,cls][nas]=sub}
		}
	start=start+bin
	#dev.off()
	} ;return(dt.mdc.27)
}

dt.mf.n=impute.chr(mf, dt.mf, 2000000)
#heatmap(as.matrix(dt.mf.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
#heatmap(as.matrix(na.omit(dt.mf.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

dt.ma.n=impute.chr(ma, dt.ma, 2000000)
#heatmap(as.matrix(dt.macd.n[raw.female.gen27,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="imputed, female, gen27")
#heatmap(as.matrix(dt.macd.n[raw.male.gen27,]), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="imputed, male, gen27")

dt.mb.n=impute.chr(mb, dt.mb, 2000000)
#heatmap(as.matrix(dt.mb.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

dt.mdc.n=impute.chr(mdc, dt.mdc, 2000000) #binsize, the smaller it is, the less imputation
#heatmap(as.matrix(dt.mdc.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
#heatmap(as.matrix(na.omit(dt.mdc.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

dt.me.n=impute.chr(me, dt.me, 2000000)
#heatmap(as.matrix(dt.me.n), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)
#heatmap(as.matrix(na.omit(dt.me.n)), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5)

#remove weird ones
imputed.plot=data.frame(dt.ma.n[,rev(seq(1:ncol(dt.ma.n)))], matrix(NA, nrow(dt.ma.n), 100),dt.mdc.n, matrix(NA, nrow(dt.ma.n), 100),dt.mb.n[,rev(seq(1:ncol(dt.mb.n)))], matrix(NA, nrow(dt.ma.n), 200),dt.me.n, matrix(NA, nrow(dt.ma.n), 100),dt.mf.n)
rownames(imputed.plot)=paste(bgdt$Generation,bgdt$Species,bgdt$prefix, sep=".")
sorted.imp=sort(rownames(imputed.plot))#exclude the two weird samples
rownames(imputed.plot)
ordered.imp=sorted.imp[c(199:204, 206:214, 160:198, 1:159)]
plotm=as.matrix(imputed.plot[rev(ordered.imp),])
jpeg("imputated.all216.jpeg",width=15,height=8,units="in",res=500)
heatmap(plotm, Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.1)
dev.off()
#@@@@@

#Step3: kmeans clusters
na.rows.a=which(is.na(rowMeans(dt.ma.n))==1)
na.rows.dc=which(is.na(rowMeans(dt.mdc.n))==1)
na.rows.b=which(is.na(rowMeans(dt.mb.n))==1)
na.rows.e=which(is.na(rowMeans(dt.me.n))==1)
na.rows.f=which(is.na(rowMeans(dt.mf.n))==1)
na.rows=Reduce(union, list(na.rows.a, na.rows.b,na.rows.dc, na.rows.e, na.rows.f)) #, na.rows.f

a.clean=na.omit(dt.ma.n)
dim(a.clean); dim(dt.ma.n)

b.clean=na.omit(dt.mb.n)
dim(b.clean); dim(dt.mb.n)
#heatmap(as.matrix(b.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_B")

dc.clean=na.omit(dt.mdc.n)
dim(dc.clean); dim(dt.mdc.n)
#heatmap(as.matrix(dc.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_DC")

e.clean=na.omit(dt.me.n)
dim(e.clean); dim(dt.me.n)
#heatmap(as.matrix(e.clean), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5, main="Muller_E")

f.clean=na.omit(dt.mf.n)
dim(f.clean); dim(dt.mf.n)


#step 4: calculate ancestry for each clusters
cluster.ancestry=function( dc.clean, clusterN, cutoff) #input: clean hapdata muller element; number of clusters to start with
{set.seed(26)
dck=kmeans(t(dc.clean), centers=clusterN); dck
while(dck$"betweenss"/dck$"totss"<cutoff) #change this percent of variation to adjust the extent of clustering !!!!!
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
#heatmap(as.matrix(cd.cl), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5)
return(cd.cl)}
###END of the function
#run the function for each muller element
a.clust.an=cluster.ancestry( a.clean, 1, 0.7);dim(a.clean); dim(a.clust.an)
b.clust.an=cluster.ancestry( b.clean, 1, 0.7);dim(b.clean); dim(b.clust.an)
dc.clust.an=cluster.ancestry( dc.clean, 1, 0.7);dim(dc.clean); dim(dc.clust.an)
e.clust.an=cluster.ancestry( e.clean,1, 0.7); dim(e.clean); dim(e.clust.an)
f.clust.an=cluster.ancestry( f.clean,1, 0.5); dim(f.clean); dim(f.clust.an)

#calculate HI for each individuals
#find intersect of all the muller.elements (individuals that are present for all)
a.list=row.names(a.clust.an)
b.list=row.names(b.clust.an)
dc.list=row.names(dc.clust.an)
e.list=row.names(e.clust.an)
f.list=row.names(f.clust.an)
all.list.notsort=Reduce(intersect, list(a.list,b.list,dc.list,e.list,f.list)) #, f.list
check=data.frame(name=all.list.notsort, p=bgdt$prefix[-na.rows])
bg.info.clust=data.frame(name=all.list.notsort, sex=bgdt$sex.g[-na.rows], gen=bgdt$Generation[-na.rows], sp=bgdt$Species[-na.rows])
bg.info.clust$gen[is.na(bg.info.clust$gen)]=0
bg.info.clust$gen[which(bg.info.clust$gen==0)[1]]=NA
table(bg.info.clust$gen)
n0=which(bg.info.clust$gen==0)
n0.alb=intersect(n0, which(bg.info.clust$sp=="albomicans"))
n0.nas=intersect(n0, which(bg.info.clust$sp=="nasuta"))
n3.4.5=c(which(bg.info.clust$gen==3), which(bg.info.clust$gen==4), which(bg.info.clust$gen==5))
n9.10=c(which(bg.info.clust$gen==10),which(bg.info.clust$gen==9))
n11.12=c(which(bg.info.clust$gen==11),which(bg.info.clust$gen==12))
n16.17.18=c(which(bg.info.clust$gen==16),which(bg.info.clust$gen==17),which(bg.info.clust$gen==18))
n21=which(bg.info.clust$gen==21)
n27.28=c(which(bg.info.clust$gen==27),which(bg.info.clust$gen==28))
#get genotype clusters for individuals that are present for all muller elements
a.c.an=a.clust.an[all.list.notsort,]
b.c.an=b.clust.an[all.list.notsort,]
dc.c.an=dc.clust.an[all.list.notsort,]
e.c.an=e.clust.an[all.list.notsort,]
f.c.an=f.clust.an[all.list.notsort,]
all.an=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))],  dc.c.an, b.c.an[,rev(seq(1:ncol(b.c.an)))], e.c.an, f.c.an) #, f.c.an
acd.c.an=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))], rep(NA, nrow(a.c.an)), dc.c.an)
##
####.  Figure 1.  !@@@@@
all.an.plot=data.frame(a.c.an[,rev(seq(1:ncol(a.c.an)))], rep(NA, nrow(a.c.an)), dc.c.an, rep(NA, nrow(a.c.an)),b.c.an[,rev(seq(1:ncol(b.c.an)))],rep(NA, nrow(a.c.an)), e.c.an, rep(NA, nrow(a.c.an)),f.c.an) #,f.c.an
all.an.p=rbind(all.an.plot[n27.28,], rep(NA, ncol(all.an.plot)), all.an.plot[n21,],rep(NA, ncol(all.an.plot)), all.an.plot[n16.17.18,], rep(NA, ncol(all.an.plot)), all.an.plot[n11.12,],rep(NA, ncol(all.an.plot)),all.an.plot[n11.12,],rep(NA, ncol(all.an.plot)),all.an.plot[n9.10,],rep(NA, ncol(all.an.plot)),all.an.plot[n3.4.5,],rep(NA, ncol(all.an.plot)),all.an.plot[c(n0.alb,n0.nas),])
jpeg("kmeans.clust.jpeg",width=15,height=8,units="in",res=500)
heatmap(as.matrix(all.an.p), Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.3,cexCol=0.3)
dev.off()
#####. @@@@
#step 5 run LD inter-chrom LD 
col.cor <- colorRampPalette(c("aquamarine","coral1"))(500)
cor.mullers.withDC=function(mclusters, indvs, mname) #input should be a ancestry matrix of a muller element, indviduals, muller element name
{cm.f.dc=cor(cbind(mclusters[indvs,],dc.c.an[indvs,]))
heatmap(t(cm.f.dc), Colv = NA, Rowv = NA, scale="none", col=col.cor, cexRow=0.3,cexCol=0.3, main=mname)
# f.dc.mean.r={}
# for(i in 1:ncol(mclusters[indvs,]))
	# {j=c(ncol(mclusters[indvs,]):(ncol(mclusters[indvs,])+ncol(dc.c.an[indvs,])))
	# f.dc.mean.r=c(f.dc.mean.r,mean(cm.f.dc[i,j], na.rm=T))}
interm=cm.f.dc[1:(ncol(mclusters[indvs,])),(ncol(mclusters[indvs,])+1):(ncol(mclusters[indvs,])+ncol(dc.c.an[indvs,]))]
heatmap(t(interm), Colv = NA, Rowv = NA, scale="none", col=col.cor, cexRow=0.3,cexCol=0.3, main=mname)
#inter-chromosome clusters correlation test corrected by fdr
c1={}; c2={}; pval={}; r={}
for(i in 1:ncol(mclusters[indvs,])) #loop through clusters in muller X
	{for(j in 1:ncol(dc.c.an[indvs,])) #loop through clusters in muller cd
		{ pval=c(pval,cor.test(mclusters[indvs,i],dc.c.an[indvs,j])[[3]]); c1=c(c1, i); c2=c(c2, j); r=c(r, cor.test(mclusters[indvs,i],dc.c.an[indvs,j])[[4]])}}
pdat=data.frame(c1, c2, pval,r)
pdat$rank=rank(pdat$pval)
pdat$c1.c2=paste(pdat$c1, pdat$c2, sep=".")
psig=pdat[which((pdat$pval*length(pdat$rank)/(pdat$rank))<0.1),]
return(psig)}

#take out females
fn=which(bg.info.clust$sex=="F")
fn.noGen0=intersect(fn, c(which(bg.info.clust$gen==9),which(bg.info.clust$gen!=0)))
table(bg.info.clust$gen[fn])
fn9.10=intersect(fn, c(which(bg.info.clust$gen==9),which(bg.info.clust$gen==10)))
fn11.12=intersect(fn, c(which(bg.info.clust$gen==11),which(bg.info.clust$gen==12)))
fn16.17=intersect(fn, c(which(bg.info.clust$gen==16),which(bg.info.clust$gen==17)))
fn4=intersect(fn, c(which(bg.info.clust$gen==4)))
fn.bf17=intersect(fn, c(which(bg.info.clust$gen<17)))
fn.af16=intersect(fn, c(which(bg.info.clust$gen>16)))
fn18=intersect(fn, c(which(bg.info.clust$gen==18)))
fn21=intersect(fn, c(which(bg.info.clust$gen==21)))
fn27=intersect(fn, c(which(bg.info.clust$gen==27)))
fn28=intersect(fn, c(which(bg.info.clust$gen==28)))
mn=which(bg.info.clust$sex=="M")
table(bg.info.clust$gen[mn])
mn27=intersect(mn, c(which(bg.info.clust$gen==27)))
mn21=intersect(mn, c(which(bg.info.clust$gen==21)))

ld.gen=function(fn)
{pdf(paste(deparse(substitute(fn)),"mullerABEF.withCD.pdf",sep="."),width=40,height=40)
a.dc.sig.r=cor.mullers.withDC(a.c.an[,rev(seq(1:ncol(a.c.an)))], indvs=fn, mname="Muller A")
b.dc.sig.r=cor.mullers.withDC(b.c.an, indvs=fn, mname="Muller B")
e.dc.sig.r=cor.mullers.withDC(e.c.an, indvs=fn, mname="Muller E")
f.dc.sig.r=cor.mullers.withDC(f.c.an,indvs=fn, mname="Muller F")
dev.off()
return(list(a.dc.sig.r, b.dc.sig.r, e.dc.sig.r, f.dc.sig.r))}
ld.full=ld.gen(fn.noGen0)
# ld.fgen9.10=ld.gen(fn9.10);ld.fgen11.12=ld.gen(fn11.12);ld.fgen16.17=ld.gen(fn16.17)
# ld.fgen18=ld.gen(fn18);ld.fgen21=ld.gen(fn21);ld.fgen27=ld.gen(fn27);ld.fgen28=ld.gen(fn28);
# f.bf16=c(fn4, fn9.10, fn11.12, fn)

ld.f.bf17=ld.gen(fn.bf17)
f.bf17.a=ld.f.bf17[[1]]; f.bf17.b=ld.f.bf17[[2]];f.bf17.e=ld.f.bf17[[3]];f.bf17.f=ld.f.bf17[[4]]
ld.f.af16=ld.gen(fn.af16)
f.af16.a=ld.f.af16[[1]]; f.af16.b=ld.f.af16[[2]]; f.af16.e=ld.f.af16[[3]]; f.af16.f=ld.f.af16[[4]]; 

intersect(f.af16.f$c1.c2, f.bf17.f$c1.c2)
ld.mgen21=ld.gen(mn21);ld.mgen27=ld.gen(mn27);
#save.image("4.3.kmeans.clust.LD.0.7kmeans.RData")

