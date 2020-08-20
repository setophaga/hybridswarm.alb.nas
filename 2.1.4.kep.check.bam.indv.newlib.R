
setwd("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/kepluana.ref")
#kepluana.ref
library(RColorBrewer)
library(vioplot)
 #1 Reference sequence identifier
  #   2 Reference sequence length
  #   3 Number of mapped reads
   #  4 Number of placed but unmapped reads # (typically unmapped partners of mapped reads)
files <- list.files(path = "~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/kepluana.ref", pattern = "*.idxstats", full.names = T)
e=read.csv(files[34], sep="\t")
#column names of e = c("chr", "length", "mapped", "unmapped")
names={}; dd={}; Column.Barcode={}; Row.Barcode={}
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
total.mreads=sum(na.omit(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 3]))
total.fr=total.mreads/sum(na.omit(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 2]))
totalm.u=total.mreads/sum(na.omit(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 4]))
placed=sum(na.omit(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 4]))
placed.fr=sum(na.omit(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 4]))/sum(na.omit(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 2]))
add=c(mad.fr,mb.fr, mdc.fr, me.fr, mf.fr,total.fr,total.mreads, row,col, placed, placed.fr, totalm.u); #print(add)
dd=rbind(dd, add)
}
dat=data.frame(names, dd)
colnames(dat)=c("id", "ma.fr","mb.fr", "mdc.fr","me.fr","mf.fr", "total.fr","total.mreads","Row.Barcode","Column.Barcode", "placed", "placed.fr", "totalmapped.vsunmapped")
# write.csv(dat, "mapping.rate.newlib.kepref.csv")

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





