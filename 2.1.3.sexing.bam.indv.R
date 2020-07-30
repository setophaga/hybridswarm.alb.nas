
setwd("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.oldlib")
library(RColorBrewer)
library(vioplot)

files <- list.files(path = "~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.oldlib/", pattern = "*.idxstats", full.names = T)
e=read.csv(files[1], sep="\t")
#column names of e = c("chr", "length", "mapped", "unmapped")
names={}; dd={}; 
for(ind in files)
{name=strsplit(strsplit(ind,"/")[[1]][9], ".idxstats")[[1]][1]
names=c(names, name)
d=read.csv(ind, sep="\t", header=F)
ma.row=which(d[,1]=="Muller_A"); ma.fr=d[ma.row, 3]/d[ma.row, 2]
mb.row=which(d[,1]=="Muller_B"); mb.fr=d[mb.row, 3]/d[mb.row, 2]
mdc.row=which(d[,1]=="Muller_DC"); mdc.fr=d[mdc.row, 3]/d[mdc.row, 2]
me.row=which(d[,1]=="Muller_E"); me.fr=d[me.row, 3]/d[me.row, 2]
mf.row=which(d[,1]=="Muller_F"); mf.fr=d[mf.row, 3]/d[mf.row, 2]
total.mreads=sum(d[c(ma.row,mb.row,mdc.row,me.row,mf.row), 3])
dd=rbind(dd, c(ma.fr,mb.fr, mdc.fr, me.fr, mf.fr,total.mreads ))
}
dat=data.frame(names, dd)
colnames(dat)=c("id", "ma.fr","mb.fr", "mdc.fr","me.fr","mf.fr", "total.mreads")
hist(dat$total.mreads)
hist(dat$mb.fr/dat$ma.fr)
plot(dat$mb.fr/dat$ma.fr,dat$total.mreads)
plot(dat$mdc.fr/dat$ma.fr,dat$total.mreads)
plot(dat$me.fr/dat$ma.fr,dat$total.mreads)
autofr=(dat$mb.fr+dat$mdc.fr+dat$me.fr)/3 #note: if include F, it is weird +dat$mf.fr
sexfr=dat$ma.fr
sex.auto.fr=sexfr/autofr

muller.col=c( "cornflowerblue","gold", "turquoise","coral", "lightgreen") #color for muller a, cd, b, e
rgb.convert=function(col.vect,fr)
{cols={}
for(i in col.vect)
	{cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}
muller.element=c("MullerA", "MullerDC", "MullerB","MullerE",  "MullerF" )
vioplot(dat$ma.fr, dat$mdc.fr,dat$mb.fr,dat$me.fr,dat$mf.fr,col=rgb.convert(muller.col, 0.5), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Mapped reads/Chr length", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(dat$ma.fr, dat$mdc.fr,dat$mb.fr,dat$me.fr,dat$mf.fr), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 0.5))  



color.gradient <- function(x, colors=c("forestgreen","gold"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )}
#plot A/(B, DC, E)
plot(sex.auto.fr,dat$total.mreads, xlim=c(0.4, 1.2), col=color.gradient(sex.auto.fr), pch=16, xlab="Muller_A/(Muller_B.DC.E)", ylab="Mapped Reads")
abline(v=0.5, lty=2);abline(v=1, lty=2)
abline(v=0.7, lty=2, col="forestgreen", lwd=2);abline(v=.9, lty=2, col="gold", lwd=2)
arrows(x=0.7, y0=quantile(dat$total.mreads, 0.99), x1 = 0.63, y1 = quantile(dat$total.mreads,.99), col="forestgreen", lwd=2, cex=0.7) #, angle = 30
arrows(x=0.9, y0=quantile(dat$total.mreads, 0.99), x1 = 0.96, y1 = quantile(dat$total.mreads,.99), col="gold", lwd=2, cex=0.7) #, angle = 30

##Muller A/B
plot(dat$ma.fr/dat$mb.fr,dat$total.mreads, xlim=c(0.4, 1.2), col=color.gradient(dat$ma.fr/dat$mb.fr), pch=16, xlab="Muller_A/Muller_B", ylab="Mapped Reads")
abline(v=0.5, lty=2);abline(v=1, lty=2)

##Muller A/DC
plot(dat$ma.fr/dat$mdc.fr,dat$total.mreads, xlim=c(0.4, 1.2), col=color.gradient(dat$mdc.fr/dat$mb.fr), pch=16, xlab="Muller_A/Muller_DC", ylab="Mapped Reads")
abline(v=0.5, lty=2);abline(v=1, lty=2)

##Muller A/E
plot(dat$ma.fr/dat$me.fr,dat$total.mreads, xlim=c(0.4, 1.2), col=color.gradient(dat$me.fr/dat$mb.fr), pch=16, xlab="Muller_A/Muller_E", ylab="Mapped Reads")
abline(v=0.5, lty=2);abline(v=1, lty=2)
##Muller A/F
plot(dat$ma.fr/dat$mf.fr,dat$total.mreads, xlim=c(0.2, 2.5), col=color.gradient(dat$mf.fr/dat$mb.fr), pch=16, xlab="Muller_A/Muller_F", ylab="Mapped Reads")
abline(v=0.5, lty=2);abline(v=1, lty=2)



#check if 
 d=read.csv("~/Desktop/hybridswarm/nasutaXabomincans.all.simple.csv") #match generation information for each individual
 d$Gender[which(d$Gender=="F")]="female"; d$Gender=factor(d$Gender, levels=c("female", "male", "ND"));table(d$Gender)
 d28=d[which(d$Generation==28),];table(d$Gender,d$Generation)
bg={}
for( i in 1:length(names))
	{row=which(names[i]==d$prefix)
	bg=rbind(bg, d[row,])
	}
bgd=data.frame(names, bg)
droplevels(bgd$Gender)
plot(sex.auto.fr~bgd$Gender)
col.to.rgb=function(col.vect, fr)
	{rgbcol={}
	for(col in col.vect)
		{color=col2rgb(col)/255
		rgbcol=c(rgbcol, rgb(color[1,1],color[2,1],color[3,1], fr))
		}
	return(rgbcol)
	}
boxplot(sex.auto.fr~bgd$Gender, col=col.to.rgb(c( "gold","forestgreen", "grey"),0.5),ylab="Muller_A/(Muller_B.DC.E)", xlab="Gender")
stripchart(sex.auto.fr~bgd$Gender,  vertical = TRUE, method = "jitter", pch = 16, col =col.to.rgb(c( "gold","forestgreen", "grey"),0.4),  bg=as.numeric(bgd$Gender), add = TRUE)

dat=data.frame(dat, sex.auto.fr)
dat$sex.g=rep(NA, nrow(dat))
dat$sex.g[which(dat$sex.auto.fr>0.85)]="F"
dat$sex.g[which(dat$sex.auto.fr<0.75)]="M"

#write.csv(dat, "~/Desktop/hybridswarm/old.lib.mullers.mapreads.csv")
ss=read.csv("~/Desktop/hybridswarm/old.lib.mullers.mapreads.csv")
ss$id[which(is.na(ss$sex.g))]


