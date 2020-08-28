#kepluana.ref
library(RColorBrewer)
library(vioplot)
 #1 Reference sequence identifier
  #   2 Reference sequence length
  #   3 Number of mapped reads
   #  4 Number of placed but unmapped reads # (typically unmapped partners of mapped reads)
# files <- list.files(path = "~/Desktop/hybridswarm/1.alb.nas.backgroundinfo/testing.athref.align", pattern = "*.idxstats", full.names = T)
# files <- list.files(path = "~/Desktop/hybridswarm/1.alb.nas.backgroundinfo/testing.kepref.align", pattern = "*.idxstats", full.names = T)
e=read.csv(files[4], sep="\t")
#column names of e = c("chr", "length", "mapped", "unmapped")
names={}; dd={}; Column.Barcode={}; Row.Barcode={}
for(ind in files)
{name=strsplit(strsplit(ind,"/")[[1]][8], ".idxstats")[[1]][1]
names=c(names, name); barcode=strsplit(name, "_")[[1]][2]
row=as.numeric(as.character(strsplit(strsplit(barcode, "-")[[1]][2],"S")[[1]][2]))
col=as.numeric(as.character(strsplit(strsplit(barcode, "-")[[1]][1],"N")[[1]][2]))
d=read.csv(ind, sep="\t", header=F)
total.mreads=sum(na.omit(d[, 3]))
total.unmapped=sum(na.omit(d[, 4]))
add=c(total.mreads, total.unmapped); #print(add)
dd=rbind(dd, add)
}
dat=data.frame(names, dd)
colnames(dat)=c("id", "totalmapped", "totalunmapped")
# write.csv(dat, "~/Desktop/hybridswarm/1.alb.nas.backgroundinfo/testing.athref.align/mapping.rate.testing.athref.csv")
#write.csv(dat, "~/Desktop/hybridswarm/1.alb.nas.backgroundinfo/testing.kepref.align/mapping.rate.testing.kepref.csv")

ath=read.csv("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/athabasca.ref/mapping.rate.newlib.athref.csv") #import ath map info
kep=read.csv("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/kepluana.ref/mapping.rate.newlib.kepref.csv")
testing.ath=read.csv("~/Desktop/hybridswarm/1.alb.nas.backgroundinfo/testing.athref.align/mapping.rate.testing.athref.csv")
testing.kep=read.csv("~/Desktop/hybridswarm/1.alb.nas.backgroundinfo/testing.kepref.align/mapping.rate.testing.kepref.csv")
testing.ath=testing.ath[1:4,];testing.kep=testing.kep[testing.kep$id %in% testing.ath$id,]
testing.kep$totalmapped.fr=testing.kep$totalmapped/(testing.kep$totalunmapped+testing.kep$totalmapped)
testing.ath$totalmapped.fr=testing.ath$totalmapped/(testing.ath$totalunmapped+testing.ath$totalmapped)

data.frame(testing.ath$id, testing.ath$totalmapped, testing.kep$totalmapped)
#mapping ratio ath vs 
hist(ath$totalmapped.fr, breaks=30)
abline(v=testing.ath$totalmapped/(testing.ath$totalmapped+testing.ath$totalunmapped), col="purple")

hist(kep$totalmapped.fr, breaks=30)
abline(v=testing.kep$totalmapped/(testing.kep$totalmapped+testing.kep$totalunmapped), col="purple")


#m=lm(ath$totalmapped.fr~kep$totalmapped.fr); summary(m)
plot( kep$total.mreads,ath$total.mreads, xlab="kepluana ref. #mapped reads", ylab="athabascan ref. #mapped reads", pch=16, col="lightgrey")#, xlim=c(0, 61273314), ylim=c(0,12898658)
out=which((ath$total.mreads/kep$total.mreads)>0.3)
#abline(m)
text(kep$total.mreads[out],ath$total.mreads[out], paste(kep$Row.Barcode[out], kep$Column.Barcode[out], sep="-"), cex=0.7)
points(kep$total.mreads[intersect(which(kep$Row.Barcode==504),which(kep$Column.Barcode==703))],ath$total.mreads[intersect(which(kep$Row.Barcode==504),which(kep$Column.Barcode==703))],col="green", pch=16)
points(kep$total.mreads[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==702))],ath$total.mreads[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==702))],col="green", pch=16)
points(kep$total.mreads[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==705))],ath$total.mreads[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==705))],col="green", pch=16)
points(kep$total.mreads[intersect(which(kep$Row.Barcode==503),which(kep$Column.Barcode==708))],ath$total.mreads[intersect(which(kep$Row.Barcode==503),which(kep$Column.Barcode==708))],col="green", pch=16)
points(kep$total.mreads[intersect(which(kep$Row.Barcode==509),which(kep$Column.Barcode==701))],ath$total.mreads[intersect(which(kep$Row.Barcode==509),which(kep$Column.Barcode==701))],col="green", pch=16)
points( testing.kep$totalmapped, testing.ath$totalmapped, col=rgb(160/255, 32/255, 240/255, 0.8), pch=16)
range(testing.kep$totalmapped)
range(testing.ath$totalmapped)

###fraction
plot( kep$totalmapped.fr,ath$totalmapped.fr, xlab="kepluana ref. fraction mapped", ylab="athabascan ref. fraction mapped", pch=16, col="lightgrey")
#out=which((ath$totalmapped.fr/kep$totalmapped.fr)>0.3)
#text(kep$totalmapped.fr[out],ath$totalmapped.fr[out], paste(kep$Row.Barcode[out], kep$Column.Barcode[out], sep="-"), cex=0.6)
points(kep$totalmapped.fr[intersect(which(kep$Row.Barcode==504),which(kep$Column.Barcode==703))],ath$totalmapped.fr[intersect(which(kep$Row.Barcode==504),which(kep$Column.Barcode==703))],col="green", pch=16)
points(kep$totalmapped.fr[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==702))],ath$totalmapped.fr[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==702))],col="green", pch=16)
points(kep$totalmapped.fr[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==705))],ath$totalmapped.fr[intersect(which(kep$Row.Barcode==505),which(kep$Column.Barcode==705))],col="green", pch=16)
points(kep$totalmapped.fr[intersect(which(kep$Row.Barcode==503),which(kep$Column.Barcode==708))],ath$totalmapped.fr[intersect(which(kep$Row.Barcode==503),which(kep$Column.Barcode==708))],col="green", pch=16)
points(kep$totalmapped.fr[intersect(which(kep$Row.Barcode==509),which(kep$Column.Barcode==701))],ath$totalmapped.fr[intersect(which(kep$Row.Barcode==509),which(kep$Column.Barcode==701))],col="green", pch=16)
points( testing.kep$totalmapped.fr, testing.ath$totalmapped.fr, col="purple", pch=16)
legend("bottomleft", c("ath. wells", "positive control alb"), pch=16, col=c("green", "purple"), cex=.8)
outfr=intersect(which(ath$totalmapped.fr>0.33), which(kep$totalmapped.fr<0.8))
text(kep$totalmapped.fr[outfr],ath$totalmapped.fr[outfr], paste(kep$Row.Barcode[outfr], kep$Column.Barcode[outfr]), cex=0.5)

