
setwd("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/kepluana.ref")
#kepluana.ref
library(RColorBrewer)
library(vioplot)
 #1 Reference sequence identifier
  #   2 Reference sequence length
  #   3 Number of mapped reads
   #  4 Number of placed but unmapped reads # (typically unmapped partners of mapped reads)
p1alb.files <- list.files(path = "~/Desktop/hybridswarm/newlibrary/p1vsp3.redemulti/p1.albref", pattern = "*.idxstats", full.names = T)
p3alb.files <- list.files(path = "~/Desktop/hybridswarm/newlibrary/p1vsp3.redemulti/p3.albref", pattern = "*.idxstats", full.names = T)
p1sul.files <- list.files(path = "~/Desktop/hybridswarm/newlibrary/p1vsp3.redemulti/p1.sulref", pattern = "*.idxstats", full.names = T)
p3sul.files <- list.files(path = "~/Desktop/hybridswarm/newlibrary/p1vsp3.redemulti/p3.sulref", pattern = "*.idxstats", full.names = T)

e=read.csv(p1alb.files[3], sep="\t")
#column names of e = c("chr", "length", "mapped", "unmapped")
mapd=function(files)
{names={}; dd={}; Column.Barcode={}; Row.Barcode={}
for(ind in files)
{name=strsplit(strsplit(ind,"/")[[1]][9], ".idxstats")[[1]][1]
names=c(names, name); barcode=strsplit(name, "_")[[1]][2]
row=as.numeric(as.character(strsplit(strsplit(barcode, "-")[[1]][2],"S")[[1]][2]))
col=as.numeric(as.character(strsplit(strsplit(barcode, "-")[[1]][1],"N")[[1]][2]))
d=read.csv(ind, sep="\t", header=F)
ma.row=which(d[,1]=="Muller_A"); ma.fr=d[ma.row, 3]/d[ma.row, 2]
mb.row=which(d[,1]=="Muller_B"); mb.fr=d[mb.row, 3]/d[mb.row, 2]
mdc.row=which(d[,1]=="Muller_DC"); mdc.fr=d[mdc.row, 3]/d[mdc.row, 2]
me.row=which(d[,1]=="Muller_E"); me.fr=d[me.row, 3]/d[me.row, 2]
mf.row=which(d[,1]=="Muller_F"); mf.fr=d[mf.row, 3]/d[mf.row, 2]
total.mreads=sum(na.omit(d[, 3]))
total.unmreads=sum(na.omit(d[, 4]))
totalm.fr=total.mreads/(total.unmreads+total.mreads)
unmapped=sum(na.omit(d[, 4]))
add=c(ma.fr,mb.fr, mdc.fr, me.fr, mf.fr,total.mreads, row,col,total.unmreads, totalm.fr); #print(add)
dd=rbind(dd, add)
}
dat=data.frame(names, dd)
colnames(dat)=c("id", "ma.fr","mb.fr", "mdc.fr","me.fr","mf.fr", "total.mreads","Row.Barcode","Column.Barcode", "unmapped", "totalmapped.fr")
return(dat)}
# write.csv(dat, "~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/kepluana.ref/mapping.rate.newlib.kepref.csv")
p1alb.d=mapd(p1alb.files)
p3alb.d=mapd(p3alb.files)
p1sul.d=mapd(p1sul.files)
p3sul.d=mapd(p3sul.files)
p1sul.d=p1sul.d[-which(p1sul.d$Row.Barcode==503),]
p1alb.d=p1alb.d[-which(p1alb.d$Row.Barcode==503),]
plot(p1alb.d$totalmapped.fr,p1sul.d$totalmapped.fr, xlim=c(0.8,.98),ylim=c(0.8,.98),pch=16,col=rgb(1,0,1,0.3),cex=1.2, xlab="mapping rate to albomicans ref", ylab="mapping rate to sulfurigaster ref")
abline(0,1)
points(p3alb.d$totalmapped.fr,p3sul.d$totalmapped.fr,col=rgb(0,1,0,0.5),add=T,cex=1.2,pch=16)
legend("topleft",pch=16, c("plate1", "plate3"),col=c(rgb(1,0,1,0.3),rgb(0,1,0,0.5)),cex=1.2)
p1.weird=intersect(which(p1alb.d$totalmapped.fr<0.6), which(p1sul.d$totalmapped.fr>0.8))
text(p1alb.d$totalmapped.fr[p1.weird],p1sul.d$totalmapped.fr[p1.weird], paste(p1alb.d$Row.Barcode[p1.weird], p1alb.d$Column.Barcode[p1.weird]), sep=".")


rows=rev(c(501, 502, 504, 505, 506, 507, 508, 509))
cols=c(701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712)
fill.plate=function(vect)
{m=matrix(nrow=8, ncol=12)
rownames(m)=rows; colnames(m)=cols
for(i in 1:length(rows))
{row=which(dat$Row.Barcode==rows[i])
for(j in 1:length(cols))
	{col=which(dat$Column.Barcode[row]==cols[j])
		if(length(col)==1)
			{m[i,j]=vect[intersect(row, which(dat$Column.Barcode==cols[j]))]}
		else{m[i,j]=NA}
	}
}
return(m)}
m.treads=fill.plate(dat$total.mreads)
m.fr=fill.plate(dat$total.fr)
m.placed=fill.plate(dat$placed)
m.placed.fr=fill.plate(dat$placed.fr)
m.mu=fill.plate(dat$totalmapped.vsunmapped)
col.well <- colorRampPalette(c("limegreen","coral"))(100)
heatmap(m.treads , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)
heatmap(m.fr , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)

heatmap(m.placed , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)
heatmap(m.placed.fr , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)

heatmap(m.mu , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)


hist(dat$totalmapped.vsunmapped,breaks=30)
hist(dat$totalmapped.vsunmapped,breaks=30, add=T)

ath=read.csv("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/athabasca.ref/mapping.rate.newlib.athref.csv") #import ath map info
dat$id ==ath$id #check each indiv matches
ath.kep.fr=ath$totalmapped.vsunmapped/dat$totalmapped.vsunmapped #fr of mapped in ath over kep ref 
ath.kep=ath$total.mreads/dat$total.mreads
m.ath.kep=fill.plate(ath.kep)
m.ath.kep.fr=fill.plate(ath.kep.fr)
col.well <- colorRampPalette(c("turquoise","hotpink"))(100)
heatmap(m.ath.kep , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)

heatmap(m.ath.kep , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)
heatmap(m.ath.kep.fr , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)

dat$prediction=min(m.ath.kep.fr, na.rm=T)
dat$prediction[intersect(which(dat$Column.Barcode==703),which(dat$Row.Barcode==504))]= max(m.ath.kep.fr, na.rm=T)
dat$prediction[intersect(which(dat$Column.Barcode==702),which(dat$Row.Barcode==505))]= max(m.ath.kep.fr, na.rm=T)
dat$prediction[intersect(which(dat$Column.Barcode==705),which(dat$Row.Barcode==505))]= max(m.ath.kep.fr, na.rm=T)
dat$prediction[intersect(which(dat$Column.Barcode==703),which(dat$Row.Barcode==508))]= max(m.ath.kep.fr, na.rm=T)
dat$prediction[intersect(which(dat$Column.Barcode==701),which(dat$Row.Barcode==509))]= max(m.ath.kep.fr, na.rm=T)
m.p=fill.plate(dat$prediction)
heatmap(m.p , Colv = NA, Rowv = NA, scale="none", col=col.well, cexRow=1)

hist(ath$total.mreads, breaks=30, col=rgb(0,1, 0, 0.6),xlab="mapped/unmapped reads to ref", xlim=c(0, max(max(ath$total.mreads), max(dat$total.mreads))))
hist(dat$total.mreads, breaks=60, col=rgb(0,0,1, 0.5),  add=T)
legend("topright",c("athabascans.ref", "kepluana.ref"), fill=c(rgb(0,1, 0, 0.6), rgb(0,0,1, 0.5)))

hist(dat$totalmapped.fr)
m=lm(ath$totalmapped.fr~dat$totalmapped.fr); summary(m)
plot( dat$total.mreads,ath$total.mreads, xlab="kepluana ref. #mapped reads", ylab="athabascan ref. #mapped reads", pch=16, col="lightgrey")
out=which((ath$total.mreads/dat$total.mreads)>0.3)
abline(m)
text(dat$total.mreads[out],ath$total.mreads[out], paste(dat$Row.Barcode[out], dat$Column.Barcode[out], sep="-"), cex=0.8)
points(dat$total.mreads[intersect(which(dat$Row.Barcode==504),which(dat$Column.Barcode==703))],ath$total.mreads[intersect(which(dat$Row.Barcode==504),which(dat$Column.Barcode==703))],col="green", pch=16)
points(dat$total.mreads[intersect(which(dat$Row.Barcode==505),which(dat$Column.Barcode==702))],ath$total.mreads[intersect(which(dat$Row.Barcode==505),which(dat$Column.Barcode==702))],col="green", pch=16)
points(dat$total.mreads[intersect(which(dat$Row.Barcode==505),which(dat$Column.Barcode==705))],ath$total.mreads[intersect(which(dat$Row.Barcode==505),which(dat$Column.Barcode==705))],col="green", pch=16)
points(dat$total.mreads[intersect(which(dat$Row.Barcode==503),which(dat$Column.Barcode==708))],ath$total.mreads[intersect(which(dat$Row.Barcode==503),which(dat$Column.Barcode==708))],col="green", pch=16)
points(dat$total.mreads[intersect(which(dat$Row.Barcode==509),which(dat$Column.Barcode==701))],ath$total.mreads[intersect(which(dat$Row.Barcode==509),which(dat$Column.Barcode==701))],col="green", pch=16)
