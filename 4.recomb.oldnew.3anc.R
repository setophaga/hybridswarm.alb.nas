library(vioplot)
f=read.csv("~/Desktop/hybridswarm/sum.switch.oldnewlib.mullercd.3anc.female.csv")
m=read.csv("~/Desktop/hybridswarm/sum.switch.3anc.mullerced.oldnewlib.male.csv")
f[which(f$mullercd.rec>15),]
f=f[-which(f$mullercd.rec>15),]
f$gen[is.na(f$gen)]=0
m$gen[is.na(m$gen)]=0
f=f[rowSums(f[,5:10])==1,]
m=m[rowSums(m[,5:10])==1,]
f$alb.fra=(f$albX.albX+f$albX.albY+f$albX.nas/2+f$albY.nas/2)
f$nas.fra=(f$nas.nas+f$albX.nas/2+f$albY.nas/2)
plot(f$gen, f$alb.fra, pch=16, col=rgb(0,0.1,0.9,0.3),cex=2, xlab="Generations", ylab="albomicans proportion")
abline(h=0.5, lwd=2, lty=2)
hist(f$alb.fra, breaks=40, col=rgb(0,0.1,0.9,0.3))
plot(f$gen, f$nas.fra, pch=16, col=rgb(0,0.7,0.8,0.3),cex=2, xlab="Generations", ylab="nasuta proportion")
abline(h=0.5, lwd=2, lty=2)
hist(f$nas.fra, breaks=40, col=rgb(0,0.7,0.8,0.3))

m$alb.fra=(m$albX.albY+m$albX.albX+m$albX.nas/2+m$albY.nas/2)
m$nas.fra=(m$nas.nas+m$albX.nas/2+m$albY.nas/2)

plot(m$gen, m$alb.fra, pch=16, col=rgb(0,0.1,0.9,0.3),cex=2, xlab="Generations", ylab="albomicans proportion")
abline(h=0.5, lwd=2, lty=2)
hist(m$alb.fra, breaks=40, col=rgb(0,0.1,0.9,0.3))
plot(m$gen, m$nas.fra, pch=16, col=rgb(0,0.7,0.8,0.3),cex=2, xlab="Generations", ylab="nasuta proportion")
abline(h=0.5, lwd=2, lty=2)
hist(m$nas.fra, breaks=40, col=rgb(0,0.7,0.8,0.3))


fo=f[order(f$alb.fra),]
f.perc=as.matrix(fo[,5:10])
f.perct=t(f.perc[rowSums(f.perc)>0,])
barplot(f.perct, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),legend =NULL)

mo=m[order(m$alb.fra),]
m.perc=as.matrix(mo[,5:10])
m.perct=t(m.perc[rowSums(m.perc)>0,])
barplot(m.perct, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),legend =NULL)

#sorted by generation
fg=f[order(f$gen),]
f.perc=as.matrix(fg[,5:10])
f.perct=t(f.perc[rowSums(f.perc)>0,])
colnames(f.perct)=fg$gen
barplot(f.perct, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),cex.name=0.3)

#sorted by generation
mg=m[order(m$gen),]
m.perc=t(as.matrix(mg[,5:10]))
colnames(m.perc)=mg$gen
barplot(m.perc, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),cex.name=0.5)

plot(x=NULL, y=NULL)
legend("topright", legend=colnames(f.perc), fill=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"))

plot(f$alb.fra, f$mullercd.rec, pch=16, col=rgb(0.8, 0, 0.9,0.2), cex=2, xlab="albomicans admixture proportion", ylab="ancestry switches", ylim=c(0, 12))
points(m$alb.fra, m$mullercd.rec, pch=16, col=rgb(0, 0.9, 0.4,0.2), cex=2, xlab="albomicans admixture proportion", ylab="ancestry switches")

plot(f$gen,f$mullercd.rec, pch=16, col=rgb(0.8, 0, 0.9,0.2), cex=2)
plot(m$gen,m$mullercd.rec, pch=16, col=rgb(0, 0.9, 0.4,0.2), cex=2)

rgb.convert=function(col.vect,fr)
{cols={};for(i in col.vect){cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}
muller.element=c("MullerA", "MullerB","MullerDC", "MullerE",  "MullerF" )
muller.col=c( "cornflowerblue","gold", "turquoise","coral", "lightgreen") #color for muller a, cd, b, e

plot(d$Generation,d$hi.a, pch=16,col=rgb.convert(muller.col[1], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller A", ylim=c(0,1))
m=lm(d$hi.a~d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);

plot(d$Generation,d$hi.b, pch=16,col=rgb.convert(muller.col[2], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller B", ylim=c(0,1))
m=lm(d$hi.b~d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);

plot(d$Generation,d$hi.dc, pch=16,col=rgb.convert(muller.col[3], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller CD", ylim=c(0,1))
m=lm(d$hi.dc~d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);

plot(d$Generation,d$hi.e, pch=16,col=rgb.convert(muller.col[4], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller E", ylim=c(0,1))
m=lm(d$hi.e~d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);

plot(d$Generation,d$hi.f, pch=16,col=rgb.convert(muller.col[5], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller F", ylim=c(0,1))
m=lm(d$hi.f~d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);

plot(d$Generation,d$HI, pch=16,col=rgb(0,0, 0, 0.1), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Genome-wide", ylim=c(0,1))
m=lm(d$HI~d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);


d=read.csv("/Users/siluwang/Desktop/hybridswarm/plot.sep2019/alb03Xnas00.hybrids.3anc.backgroundinfo.gw.newoldlib.sep2020.csv")
d=d[-which(d$mullercd.rec>60),]
d=d[-which(is.na(d$Species)==1),]
d$Generation[is.na(d$Generation)]=0
d=d[-which(d$Generation==0),]
plot(d$hi.a,d$het.a, pch=16, col=rgb.convert(muller.col[1], 0.3), cex=1.5, xlab="Admixture proportion, Muller A", ylab="Heterozygosity, Muller A", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(d$hi.b,d$het.b, pch=16, col=rgb.convert(muller.col[2], 0.3), cex=1.5, xlab="Admixture proportion, Muller B", ylab="Heterozygosity, Muller B", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(d$hi.dc,d$het.dc, pch=16, col=rgb.convert(muller.col[3], 0.3), cex=1.5, xlab="Admixture proportion, Muller CD", ylab="Heterozygosity, Muller CD", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(d$hi.e,d$het.e, pch=16, col=rgb.convert(muller.col[4], 0.3), cex=1.5, xlab="Admixture proportion, Muller E", ylab="Heterozygosity, Muller E", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(d$hi.f,d$het.f, pch=16, col=rgb.convert(muller.col[5], 0.3), cex=1.5, xlab="Admixture proportion, Muller F", ylab="Heterozygosity, Muller F", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(d$HI,d$het, pch=16, col=rgb(0,0,0, 0.1), cex=1.5, xlab="Admixture proportion", ylab="Heterozygosity", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)

plot(d$Generation, d$het, pch=16, col=rgb(0,0,0,0.2),cex=2)
m=lm(d$het~d$Generation); summary(m)
abline(m)

plot(d$HI, d$hi.dc, xlim=c(0,1), ylim=c(0,1), pch=16, col=rgb(0,0.3, 0.9, 0.3), cex=2, xlab="Genomic ancestry proportion", ylab="Muller CD ancestry proportion")
abline(0,1)

alb.intro=function(p, h)
	{print((1-h))
	y=(2*p-1)/(1-h)
	zero=which((1-h)==0)
	y[zero]=0
	return(y)} #where p = hi, and h=het
d$alb.intro=alb.intro(d$HI, d$het)
plot(d$HI,d$alb.intro)
d$alb.intro.cd=alb.intro(d$hi.dc, d$het.dc)
d$alb.intro.a=alb.intro(d$hi.a, d$het.a)
d$alb.intro.b=alb.intro(d$hi.b, d$het.b)
d$alb.intro.e=alb.intro(d$hi.e, d$het.e)
d$alb.intro.f=alb.intro(d$hi.f, d$het.f)
plot(d$hi.dc, d$mullercd.rec*d$het)

plot(d$Generation,d$alb.intro.a, pch=16,col=rgb.convert(muller.col[1], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller A", ylim=c(-1,1))
m=lm(d$alb.intro.a~d$Generation+0); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

plot(d$Generation,d$alb.intro.b, pch=16,col=rgb.convert(muller.col[2], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller B", ylim=c(-1,1))
m=lm(d$alb.intro.b~d$Generation+0); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

plot(d$Generation,d$alb.intro.cd, pch=16,col=rgb.convert(muller.col[3], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller CD", ylim=c(-1,1))
m=lm(d$alb.intro.cd~d$Generation+0); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

plot(d$Generation,d$alb.intro.e, pch=16,col=rgb.convert(muller.col[4], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller E", ylim=c(-1,1))
m=lm(d$alb.intro.e~d$Generation+0); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

plot(d$Generation,d$alb.intro.f, pch=16,col=rgb.convert(muller.col[5], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller F", ylim=c(-1,1))
m=lm(d$alb.intro.f~d$Generation-1); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

plot(d$Generation,d$alb.intro, pch=16,col=rgb(0,0, 0, 0.1), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Genome-wide", ylim=c(-1,1))
m=lm(d$alb.intro~0+d$Generation); summary(m)
abline(m, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

##albomicans introgression
vioplot(d$alb.intro.a,d$alb.intro.b,d$alb.intro.cd,d$alb.intro.e,d$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(d$alb.intro.a,d$alb.intro.b,d$alb.intro.cd,d$alb.intro.e,d$alb.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(d$alb.intro.a,d$alb.intro.b,d$alb.intro.cd,d$alb.intro.e,d$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)


#recombination rates
vioplot(d$mullera.rec,d$mullerb.rec,d$mullercd.rec,d$mullere.rec,d$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(d$mullera.rec,d$mullerb.rec,d$mullercd.rec,d$mullere.rec,d$mullerf.rec ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(d$mullera.rec,d$mullerb.rec,d$mullercd.rec,d$mullere.rec,d$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)

#recombination rates corrected by admixture
vioplot(d$mullera.rec*d$het,d$mullerb.rec*d$het,d$mullercd.rec*d$het,d$mullere.rec*d$het,d$mullerf.rec*d$het ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(d$mullera.rec*d$het,d$mullerb.rec*d$het,d$mullercd.rec*d$het,d$mullere.rec*d$het,d$mullerf.rec*d$het ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(d$mullera.rec*d$het,d$mullerb.rec*d$het,d$mullercd.rec*d$het,d$mullere.rec*d$het,d$mullerf.rec*d$het ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)

#number of sites
#mullerA, mullerB, mullerCD, mullerE, muller F informative sites
info.sites=c(128670, 120223, 188404, 126724, 3036)
info.fr=c(128670, 120223, 188404, 126724, 3036)/sum(c(128670, 120223, 188404, 126724, 3036))