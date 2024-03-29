load("x1.haplotypes.sp.new.old.combo.RData")
#install.packages("sommer")
library(sommer)
library(DescTools)#Mode function for imputation
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

mdc=dat.sp$pos[which(dat.sp$chrom=="Muller_DC")]
ma=dat.sp$pos[which(dat.sp$chrom=="Muller_A")]
mb=dat.sp$pos[which(dat.sp$chrom=="Muller_B")]
me=dat.sp$pos[which(dat.sp$chrom=="Muller_E")]
mf=dat.sp$pos[which(dat.sp$chrom=="Muller_F")]

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
heatmap(as.matrix(plotm), Colv = NA, Rowv = NA, scale="none", col=col.sp, cexRow=0.1)
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


# load("4.3.kmeans.clust.LD.0.7kmeans.RData")
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
	{#print(dck$cluster[i])
	if(dck$cluster[i]%in%clusters)
		{n=n+1
		tochange=dck$cluster[i]
		all=as.numeric(as.character(which(dck$cluster==tochange)))
		want=all[c(which(all>i),which(all==i))]
		cluster.new[want]=range(dck$cluster)[2]+n
		# print("ahhh");print(dck$cluster[i])
		#print(range(dck$cluster)[2]+n)
		clusters=c(clusters, range(dck$cluster)[2]+n)}
	else{clusters=c(clusters, dck$cluster[i])}
	}
}
#print(length(clusters)) #the right order for clusters along the chromosome
#print(table(cluster.new)) #see partitions of new cluster information
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
return(list(cd.cl, table(cluster.new)))}
###END of the function
#run the function for each muller element
a.clust.results=cluster.ancestry( a.clean, 1, 0.5); a.clust.results[[2]]; a.clust.an=a.clust.results[[1]];dim(a.clean); dim(a.clust.an)
b.clust.results=cluster.ancestry( b.clean, 1, 0.5); b.clust.results[[2]]; b.clust.an=b.clust.results[[1]];dim(b.clean); dim(b.clust.an)
dc.clust.results=cluster.ancestry( dc.clean, 1, 0.5); dc.clust.results[[2]]; dc.clust.an=dc.clust.results[[1]];dim(dc.clean); dim(dc.clust.an)
e.clust.results=cluster.ancestry( e.clean, 1, 0.5); e.clust.results[[2]]; e.clust.an=e.clust.results[[1]]; dim(e.clean); dim(e.clust.an)
f.clust.results=cluster.ancestry( f.clean, 1, 0.5); f.clust.results[[2]]; f.clust.an=f.clust.results[[1]]; dim(f.clean); dim(f.clust.an)

#calculate HI for each individuals
#find intersect of all the muller.elements (individuals that are present for all)
a.list=row.names(a.clust.an)
b.list=row.names(b.clust.an)
dc.list=row.names(dc.clust.an)
e.list=row.names(e.clust.an)
f.list=row.names(f.clust.an)
all.list.notsort=Reduce(intersect, list(a.list,b.list,dc.list,e.list,f.list)) #, f.list
check=data.frame(name=all.list.notsort, p=bgdt$prefix[-na.rows])
bg.info.clust=data.frame(name=all.list.notsort, sex=bgdt$sex.g[-na.rows], gen=bgdt$Generation[-na.rows], sp=bgdt$Species[-na.rows], hi=bgdt$HI[-na.rows])
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

#take out females vs males and different generations
fn=which(bg.info.clust$sex=="F")
fn.bg=bg.info.clust[fn,]
fn.noGen0=intersect(fn, c(which(bg.info.clust$gen==9),which(bg.info.clust$gen!=0)))
table(bg.info.clust$gen[fn])
fn9.10=intersect(fn, c(which(bg.info.clust$gen==9),which(bg.info.clust$gen==10)))
fn11.12=intersect(fn, c(which(bg.info.clust$gen==11),which(bg.info.clust$gen==12)))
fn16.17=intersect(fn, c(which(bg.info.clust$gen==16),which(bg.info.clust$gen==17)))
fn4=intersect(fn, c(which(bg.info.clust$gen==4)))
fn5=intersect(fn, c(which(bg.info.clust$gen==5)))
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


#step 5, BGC analysis of the clusters.
#5.1 first calculate genome-wide hybrid index hi
dc.hi=rowMeans(dc.c.an)
a.hi=rowMeans(a.c.an)
b.hi=rowMeans(b.c.an)
e.hi=rowMeans(e.c.an)
f.hi=rowMeans(f.c.an)
dd.hi=data.frame(a.hi, b.hi, dc.hi, e.hi, f.hi)
hi=rowMeans(dd.hi)

bgc=function(muller.data, gen.num) #background data, genotype data with chrom.pos, geno data
{coefi={}
for(i in 1:ncol(muller.data)) #loop through each cluster
	{h=bg.info.clust$hi[gen.num]
	p=as.numeric(as.character(muller.data[gen.num,i]))
	set.seed(23)
	m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm = "port", control=nls.control(maxiter =100, warnOnly=TRUE),lower=c(a=-3, b=-3), upper=c(a=3,b=3)) #lower=c(a=-1, b=-1), upper=c(a=1,b=1),
	a=coef(m)[1]; b=coef(m)[2]
	se.a=coef(summary(m))[1, "Std. Error"];se.b=coef(summary(m))[2, "Std. Error"];
	a.lw95CI=a-1.96*se.a; a.up95CI=a+1.96*se.a; b.lw95CI=b-1.96*se.b; b.up95CI=b+1.96*se.b; 
	coefi=rbind(coefi, c(a, b, se.a, se.b, a.lw95CI,a.up95CI, b.lw95CI, b.up95CI))
	#x<-data.frame(h=seq(0, 1, length.out=1000));
	#y<-predict(m,x)
	#plot(h, p, col=rgb(0,1,1, 0.5), pch=16, cex=2, main=colnames(muller.data)[i], xlim=c(0,1), ylim=c(0,1))
	#abline(0,1, lty=2, lwd=1.5)
	#lines(x[,1], y, col="blue", lwd=2)
	#text(0.1,1, paste("alpha= ", round(a, 3)), cex=0.7)
	#text(0.1,0.95, paste("beta= ", round(b, 3)), cex=0.7)
	}
colnames(coefi)=c("a", "b", "se.a", "se.b", "a.95CI.lw", "a.95CI.up", "b.95CI.lw", "b.95CI.up")
return(data.frame(coefi))
}
#generation 28 calculate alpha, beta, and their se for each muller element
bgc.generation=function(fn28, gen)
{a.g28=bgc(a.c.an, fn28 ); #a.g28[which((a.g28$b.95CI.lw*a.g28$b.95CI.up)>0),]
b.g28=bgc(b.c.an, fn28 )
dc.g28=bgc(dc.c.an, fn28 );#dc.g28[which((dc.g28$b.95CI.lw*dc.g28$b.95CI.up)>0),]
e.g28=bgc(e.c.an, fn28 )
f.g28=bgc(f.c.an, fn28 )
g28=rbind(a.g28[rev(seq(1:nrow(a.g28))),], dc.g28, b.g28[rev(seq(1:nrow(b.g28))),],  e.g28, f.g28)
muller=c(nrow(a.g28),(nrow(dc.g28)+nrow(a.g28)), nrow(dc.g28)+(nrow(b.g28)+nrow(a.g28)), nrow(dc.g28)+nrow(e.g28)+(nrow(b.g28)+nrow(a.g28)))
g28$muller.element=c(rep("MullerA", nrow(a.g28)),rep("MullerDC", nrow(dc.g28)), rep("MullerB", nrow(b.g28)),rep("MullerE", nrow(e.g28)), rep("MullerF", nrow(f.g28)))
#plot alpha for gen 28
g28$cluster=c(seq(1:nrow(a.g28)), seq(1:nrow(dc.g28)), seq(1:nrow(b.g28)),seq(1:nrow(e.g28)),seq(1:nrow(f.g28)))
x=seq(1:length(g28[,1]));y2=g28[,2];se2=g28[,4];y1=g28[,1];se1=g28[,3]
pdf(paste("gen", gen,"alpha.beta.pdf",sep="."),width=7,height=6)
par(mfrow = c(2, 1))
#plot alpha
plot(x,y1, ylim=c(-2, 2), pch=16, ylab=expression(alpha), xlab="Relative Positions", main=paste("Generation", gen))
arrows(x, y1-1.96*se1, x, y1+1.96*se1, length=0.05, angle=90, code=3, col="forestgreen")
sig.pos1=intersect(which((g28$a.95CI.up* g28$a.95CI.lw)>0),which(g28$a.95CI.up>0))
sig.neg1=intersect(which((g28$a.95CI.up* g28$a.95CI.lw)>0),which(g28$a.95CI.up<0))
points(x[sig.pos1], y1[sig.pos1], col="gold", pch=16, cex=0.8)
points(x[sig.neg1], y1[sig.neg1], col="turquoise1", pch=16, cex=0.8)
abline(h=0, lty=2)
abline(v=muller+0.5)
#plot beta
plot(x,y2, ylim=c(-4, 4), pch=16, ylab=expression(beta), xlab="Relative Positions", main=paste("Generation", gen))
arrows(x, y2-1.96*se2, x, y2+1.96*se2, length=0.05, angle=90, code=3, col="dodgerblue4 ")
sig.pos2=intersect(which((g28$b.95CI.up* g28$b.95CI.lw)>0),which(g28$b.95CI.up>0))
sig.neg2=intersect(which((g28$b.95CI.up* g28$b.95CI.lw)>0),which(g28$b.95CI.up<0))
points(x[sig.pos2], y2[sig.pos2], col="gold", pch=16, cex=0.8)
points(x[sig.neg2], y2[sig.neg2], col="turquoise2", pch=16, cex=0.8)
abline(h=0, lty=2)
abline(v=muller+0.5)
dev.off()
return(g28)}

#table(fn.bg$gen)
fn4to10.bgc=bgc.generation(c(fn4, fn5, fn9.10), "4-10"); 
fn9to10.bgc=bgc.generation( fn9.10, "9-10")
fn11to12.bgc=bgc.generation(fn11.12, "11-12")
fn16to17.bgc=bgc.generation(c(fn16.17), "16-17")
fn16to18.bgc=bgc.generation(c(fn16.17, fn18), "16-18")
fn16to21.bgc=bgc.generation(c(fn16.17, fn18, fn21), "16-21")
fn18.bgc=bgc.generation(fn18, "18")
fn21.bgc=bgc.generation(fn21, "21"); 
fn27.bgc=bgc.generation(fn27, "27")
fn28.bgc=bgc.generation(fn28, "28")
write.csv(fn18.bgc, "fn18.bgc.csv")
write.csv(fn21.bgc, "fn21.bgc.csv")
write.csv(fn27.bgc, "fn27.bgc.csv")
write.csv(fn28.bgc, "fn28.bgc.csv")

fn27to28.bgc=bgc.generation(c(fn27, fn28), "27-28")
fn.bgc=bgc.generation(fn.noGen0, "all.female")
fn21to28.bgc=bgc.generation(c(fn21,fn27, fn28), "21-28")

#testing whether there is less rec per sites in clusters with high beta
####count ancestry turnovers within clusters
ancestry.recomb=function(indv, muller, cluster) #cluster is the a.clust.results[[2]], with number of sites for each cluster
{ #ancestry at muller cd
cluster.size={}; switches={}
for(c in names(cluster))
{nc=which(names(cluster)==c)
start.seq=sum(cluster[0:(nc-1)]);end.seq=sum(cluster[0:nc])
indv.an=spd[which(dat.sp$chrom==muller)[start.seq:end.seq],which(nnn==indv)]
positions=dat.sp[which(dat.sp$chrom==muller)[start.seq:end.seq],2]
pos.notna=positions[which(is.na(indv.an)==0)]
an.notna=na.omit(indv.an)
pos.switch=pos.notna[which(abs((an.notna-c(an.notna[-1],an.notna[1])))>0)]
if(length(pos.switch)==0)
{switch.cluster=0}
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
	switch.cluster=length(ancestry)
		};#anblock=data.frame(ancestry, begin, end)
	}
	cluster.size=c(cluster.size, cluster[nc])
	switches=c(switches, switch.cluster)
}
switchd=data.frame(switches, cluster.size)	
return(switchd)
}
anc.recomb.avg=function(indv.list) #indv.list is a list of names, e.g. bg.info.clust$name[fn21]
{rec.cd.rated={};rec.a.rated={}; rec.b.rated={};  rec.e.rated={}; rec.f.rated={}; 
for(indv in indv.list)
	{rec.a.clust=ancestry.recomb(indv, "Muller_A", a.clust.results[[2]])
	rec.cd.clust=ancestry.recomb(indv, "Muller_DC", dc.clust.results[[2]])
	rec.b.clust=ancestry.recomb(indv, "Muller_B", b.clust.results[[2]])
	rec.e.clust=ancestry.recomb(indv, "Muller_E", e.clust.results[[2]])
	rec.f.clust=ancestry.recomb(indv, "Muller_F", f.clust.results[[2]])
	rec.a.rate=rec.a.clust[,1]/rec.a.clust[,2]; rec.a.rated=cbind(rec.a.rated, rev(rec.a.rate)) 
	rec.cd.rate=rec.cd.clust[,1]/rec.cd.clust[,2]; rec.cd.rated=cbind(rec.cd.rated, rec.cd.rate) #each row is an individual, each col is a cluster
	rec.b.rate=rec.b.clust[,1]/rec.b.clust[,2]; rec.b.rated=cbind(rec.b.rated, rec.b.rate) 
	rec.e.rate=rec.e.clust[,1]/rec.e.clust[,2]; rec.e.rated=cbind(rec.e.rated, rec.e.rate) 
	rec.f.rate=rec.f.clust[,1]/rec.f.clust[,2]; rec.f.rated=cbind(rec.f.rated, rec.f.rate) }
	rec.rated=rbind(rec.a.rated, rec.cd.rated, rec.b.rated, rec.e.rated, rec.f.rated)
	colnames(rec.rated)=paste("indv", 1:length(indv.list), sep=".")
	muller.elements=c(rep("Muller_A", length(a.clust.results[[2]])),rep("Muller_DC", length(dc.clust.results[[2]])),rep("Muller_B", length(b.clust.results[[2]])), rep("Muller_E", length(e.clust.results[[2]])),rep("Muller_F", length(f.clust.results[[2]])))
	rec.rated=data.frame(rec.rated, muller.elements)
return(rec.rated)}

recd.f21=anc.recomb.avg(c(bg.info.clust$name[fn21]))
x.f21=rowMeans(recd.f21[,-ncol(recd.f21)])
ya.f21=fn21.bgc[,1];yb.f21=fn21.bgc[,2]
ma.21=lm(ya.f21~x.f21*fn21.bgc$muller.element);summary(ma.21); mb.21=lm(yb.f21~x.f21); summary(mb.21); 

recd.f27=anc.recomb.avg(bg.info.clust$name[fn27])
x.f27=rowMeans(recd.f27[,-ncol(recd.f27)])
ya.f27=fn27.bgc[,1];yb.f27=fn27.bgc[,2]
ma.27=lm(ya.f27~x.f27*fn27.bgc$muller.element);summary(ma.27); mb.27=lm(yb.f27~x.f27); summary(mb.27); 

recd.f28=anc.recomb.avg(bg.info.clust$name[fn28])
x.f28=rowMeans(recd.f28[,-ncol(recd.f28)])
ya.f28=fn28.bgc[,1];yb.f28=fn28.bgc[,2]
ma.28=lm(ya.f28~x.f28*fn28.bgc$muller.element);summary(ma.28); mb.28=lm(yb.f28~x.f28); summary(mb.28); 

recd.f21to28=anc.recomb.avg(bg.info.clust$name[c(fn21, fn27, fn28)])
x.f21to28=rowMeans(recd.f21to28[,-ncol(recd.f21to28)])
ya.f21to28=fn21to28.bgc[,1];yb.f21to28=fn21to28.bgc[,2]
ma.21to28=lm(ya.f21to28~x.f21to28);summary(ma.21to28); mb.21to28=lm(yb.f21to28~x.f21to28); summary(mb.21to28); 
summary(ma.21to28); mb.21to28=lm(yb.f21to28~x.f21to28); summary(mb.21to28); 


recd.21to28=data.frame(fn21.bgc$muller.element, x.f21, ya.f21, yb.f21, x.f27, ya.f27, yb.f27, x.f28, ya.f28, yb.f28, x.f21to28, ya.f21to28, yb.f21to28)
colnames(recd.21to28)=c("muller.element", "rec.rate21", "alpha21", "beta21", "rec.rate27", "alpha27", "beta27", "rec.rate28", "alpha28", "beta28", "rec.rate21to28", "alpha21to28","beta21to28")
write.csv(recd.21to28, "recomb.bgc.muller.csv")

### Now associate this with fn 21 muller cd beta and alpha!!!!!!!!
##############################################################################################################################
##########################################


barrier.mullercd=function(fn21.bgc)
{return(f21.b=fn21.bgc$cluster[intersect(intersect(which((fn21.bgc$b.95CI.up* fn21.bgc$b.95CI.lw)>0),which(fn21.bgc$b.95CI.up>0)), which(fn21.bgc$muller.element=="MullerDC"))])}

sig.cd=union(barrier.mullercd(fn21.bgc),union(barrier.mullercd(fn27.bgc),barrier.mullercd(fn28.bgc)))


#####. @@@@
#step 6 run LD inter-chrom LD 
col.cor <- colorRampPalette(c("aquamarine","coral1"))(500)
cor.mullers.withDC=function(mclusters, indvs, mname, sig.cd) #input should be a ancestry matrix of a muller element, indviduals, muller element name
{cm.f.dc=cor(cbind(mclusters[indvs,],dc.c.an[indvs,sig.cd])) #sig.cd is the cluster that is significantly positive for beta
heatmap(t(cm.f.dc), Colv = NA, Rowv = NA, scale="none", col=col.cor, cexRow=0.3,cexCol=0.3, main=mname)
# f.dc.mean.r={}
# for(i in 1:ncol(mclusters[indvs,]))
	# {j=c(ncol(mclusters[indvs,]):(ncol(mclusters[indvs,])+ncol(dc.c.an[indvs,])))
	# f.dc.mean.r=c(f.dc.mean.r,mean(cm.f.dc[i,j], na.rm=T))}
#interm=cm.f.dc[1:(ncol(mclusters[indvs,])),(ncol(mclusters[indvs,])+1):(ncol(mclusters[indvs,])+ncol(dc.c.an[indvs,sig.cd]))]
heatmap(t(cm.f.dc), Colv = NA, Rowv = NA, scale="none", col=col.cor, cexRow=0.3,cexCol=0.3, main=mname)
#inter-chromosome clusters correlation test corrected by fdr
c1={}; c2={}; pval={}; r={}
for(i in 1:ncol(mclusters[indvs,])) #loop through clusters in muller X
	{for(j in sig.cd) #loop through clusters in muller cd
		{ pval=c(pval,cor.test(mclusters[indvs,i],dc.c.an[indvs,j])[[3]]); c1=c(c1, i); c2=c(c2, j); r=c(r, cor.test(mclusters[indvs,i],dc.c.an[indvs,j])[[4]])}}
pdat=data.frame(c1, c2, pval,r)
pdat$rankwithinchr=rank(pdat$pval)
pdat$c1.c2=paste(pdat$c1, pdat$c2, sep=".");pdat$chrom=mname
psig=pdat[which((pdat$pval*length(pdat$rank)/(pdat$rank))<0.1),]
return(list(psig, pdat))}

ld.gen=function(fn)
{pdf(paste(deparse(substitute(fn)),"mullerABEF.withCD.pdf",sep="."),width=40,height=40)
a.dc.sig.r=cor.mullers.withDC(a.c.an[,rev(seq(1:ncol(a.c.an)))], indvs=fn, mname="Muller A", sig.cd=sig.cd)
b.dc.sig.r=cor.mullers.withDC(b.c.an, indvs=fn, mname="Muller B",sig.cd=sig.cd)
e.dc.sig.r=cor.mullers.withDC(e.c.an, indvs=fn, mname="Muller E",sig.cd=sig.cd)
f.dc.sig.r=cor.mullers.withDC(f.c.an,indvs=fn, mname="Muller F",sig.cd=sig.cd)
dev.off()
return(list(rbind(a.dc.sig.r[[2]], b.dc.sig.r[[2]], e.dc.sig.r[[2]], f.dc.sig.r[[2]]), a.dc.sig.r[[1]], b.dc.sig.r[[1]], e.dc.sig.r[[1]], f.dc.sig.r[[1]]))}

sig.ld=function(fn.noGen0){ld.full=ld.gen(fn.noGen0)
ld.fulld=ld.full[[1]]
ld.fulld$rank=rank(ld.fulld$pval)
print(ld.fulld[which((ld.fulld$pval*length(ld.fulld$rank)/(ld.fulld$rank))<0.05),])
}
# ld.fgen9.10=ld.gen(fn9.10);ld.fgen11.12=ld.gen(fn11.12);ld.fgen16.17=ld.gen(fn16.17)
# ld.fgen18=ld.gen(fn18);ld.fgen21=ld.gen(fn21);ld.fgen27=ld.gen(fn27);ld.fgen28=ld.gen(fn28);
# f.bf16=c(fn4, fn9.10, fn11.12, fn)
sig.ld(fn.noGen0)
sig.ld(fn.bf17)
sig.ld(fn.af16)
sig.ld(c(fn21, fn27,fn28) )
sig.ld(fn21);sig.ld(fn27); sig.ld(fn28)
sig.ld(c(fn21,fn27,fn28))


ld.f.bf17=ld.gen(fn.bf17)
# sig.ld(fn9.10); sig.ld(fn11.12); sig.ld(fn16.17); sig.ld(fn18); 
f.bf17.a=ld.f.bf17[[1]]; f.bf17.b=ld.f.bf17[[2]];f.bf17.e=ld.f.bf17[[3]];f.bf17.f=ld.f.bf17[[4]]
ld.f.af16=ld.gen(fn.af16)
f.af16.a=ld.f.af16[[1]]; f.af16.b=ld.f.af16[[2]]; f.af16.e=ld.f.af16[[3]]; f.af16.f=ld.f.af16[[4]]; 

intersect(f.af16.f$c1.c2, f.bf17.f$c1.c2)
ld.mgen21=ld.gen(mn21);ld.mgen27=ld.gen(mn27);
#save.image("4.3.kmeans.clust.LD.0.7kmeans.RData")

#plot networks
library(igraph)
graph_from_edgelist(, directed = TRUE)
network.plot.muller=function(m28, n1,n2,main, col1, col2)
{g28=graph_from_adjacency_matrix(as.matrix(m28), mode = "undirected", weighted = T, diag = F)
 colors.muller <- c("black",col1, col2)
 ecolor = colorRampPalette(c("black","darkolivegreen1" ))(100-cutoff*10)
plot.igraph(g28, vertex.cex=0.1, edge.color=ecolor[E(g28)$weight*100-cutoff*10],layout=layout_in_circle,edge.width=E(g28)$weight*2,vertex.size=3,vertex.color=colors.muller[c(1,rep(2, n1),rep(3, n2))],vertex.label=NA, margin=0,main=main)#(20*authority.score(g28)$vector)
print(summary(g28))
} 


