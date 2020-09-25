library(vioplot)
f=read.csv("~/Desktop/hybridswarm/sum.switch.oldnewlib.mullercd.3anc.female.csv")
m=read.csv("~/Desktop/hybridswarm/sum.switch.3anc.mullerced.oldnewlib.male.csv")
n=read.csv("~/Desktop/hybridswarm/1.pipeline/coverage.muller.hybrids.newlib/newlib96wells.coverage.csv")
o=read.csv("~/Desktop/hybridswarm/old.lib.mullers.mapreads.csv")
d=rbind(n, o[,1:8]); d$id=as.character(d$id)
mapq.match=function(sexd){
sexd$p=as.character(sexd$p);bg={}
for(i in sexd$p){row=which(d$id==i); bg=rbind(bg,d[row,])}
sexd=data.frame(sexd, bg)
return(sexd)}
f=mapq.match(f)
m=mapq.match(m)
hist(f$total.mreads)
f=f[-which(f$mullercd.rec>30),]
#remove low coverage individuals (the ones with < 50,000reaads were dropec)
f=f[-which(f$total.mreads< 50000),]
m=m[-which(m$total.mreads< 50000),]
#filter for quality
f$gen[is.na(f$gen)]=0
m$gen[is.na(m$gen)]=0
#calculate albomicans freq within muller CD with 5 haplotypes
f$alb.fra=(f$albX.albX+f$albX.albY+f$albX.nas/2+f$albY.nas/2)
f$nas.fra=(f$nas.nas+f$albX.nas/2+f$albY.nas/2)
m$alb.fra=(m$albX.albY+m$albX.albX+m$albX.nas/2+m$albY.nas/2)
m$nas.fra=(m$nas.nas+m$albX.nas/2+m$albY.nas/2)
#match genome-wide info
d=read.csv("/Users/siluwang/Desktop/hybridswarm/plot.sep2019/alb03Xnas00.hybrids.3anc.backgroundinfo.gw.newoldlib.sep2020.csv")
d$prefix=as.character(d$prefix)
#Function to add genome-wide info to f and m
gwd.match=function(f)
{bgi={}
for(i in f$p)
{row=which(i==d$prefix)
bgi=rbind(bgi,d[row,])}
f=data.frame(f, bgi)
return(f)}
f=gwd.match(f)
m=gwd.match(m)
#######################################################
###################################################.   completion of f and m matching and filter
mullercd.col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2")
mullercd.hap=c("albX.albX" ,"albX.albY" ,"albX.nas" , "albY.albY", "albY.nas" , "nas.nas"  )
rgb.convert=function(col.vect,fr)
{cols={};for(i in col.vect){cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}
#######color plotting function
#test asasociation between haploltype fraction and haplotype switcches --> expect albX and albY fraction correspond to less switches
mullercd.hap.switch=function(m){
y=m$alb.intro#mullercd.rec*m$het 
plot(m$albX.albX,y, pch=16, col=rgb.convert(mullercd.col[1], 0.5), cex=2, ylim=c(-1,1), xlim=c(0,1), xlab='haplotype fraction', ylab="Alb intro") #"MullerCD haplotype switching"
md=lm(y~m$albX.albX);print(summary(md));abline(md, lwd=2, col=mullercd.col[1]);
points(m$albX.albY,y, pch=16, col=rgb.convert(mullercd.col[2], 0.3), cex=2)
md=lm(y~m$albX.albY);abline(md, lwd=3, col=mullercd.col[2]);print(summary(md));
points(m$albX.nas,y, pch=16, col=rgb.convert(mullercd.col[3], 0.3), cex=2)
md=lm(y~m$albX.nas);abline(md, lwd=3, col=mullercd.col[3]);print(summary(md));
points(m$albY.nas,y, pch=16, col=rgb.convert(mullercd.col[5], 0.3), cex=2)
md=lm(y~m$albY.nas);abline(md, lwd=3, col=mullercd.col[5]);print(summary(md));
points(m$nas.nas,y, pch=16, col=rgb.convert(mullercd.col[6], 0.3), cex=2)
md=lm(y~m$nas.nas);abline(md, lwd=3, col=mullercd.col[6]);print(summary(md));}

mullercd.hap.switch(m)
mullercd.hap.switch(f)
plot(m$albY.nas,y, pch=16, col=rgb.convert(mullercd.col[5], 0.3), cex=2)
md=lm(y~m$albY.nas);abline(md, lwd=3, col=mullercd.col[5]);print(summary(md));
y=f$mullercd.rec*f$het
plot(f$albY.nas,f$mullercd.rec*f$het, pch=16, col=rgb.convert(mullercd.col[5], 0.3), cex=2)
md=lm(y~f$albY.nas);abline(md, lwd=3, col=mullercd.col[5]);print(summary(md));


plot.rec=function(fx, mx){
	plot(mx,m$mullercd.rec, pch=16, col=rgb(0,0.2, 0.8, 0.2), cex=2, xlim=c(0,1),ylim=c(0,21))
points(fx,f$mullercd.rec, pch=16, col=rgb(0.8,0, 0.2, 0.2), cex=2)}
plot.rec(f$albX.nas, m$X3..albX..nas.)
plot.rec(f$albX.albX, m$X5..albY..nas.)
plot(m$X5..albY..nas.,m$mullercd.rec, pch=16, col=rgb(0,0.2, 0.8, 0.2), cex=2, xlim=c(0,1),ylim=c(0,21))
points( m$X6..nas..nas.,m$mullercd.rec, pch=16, col=rgb(0,0.2, 0.8, 0.2),cex=2)

### look at muller CD alb freq and nas freq over generations 
plot(f$gen, f$alb.fra, pch=16, col=rgb(0,0.1,0.9,0.3),cex=2, xlab="Generations", ylab="albomicans proportion")
abline(h=0.5, lwd=2, lty=2)
hist(f$alb.fra, breaks=40, col=rgb(0,0.1,0.9,0.3))
plot(f$gen, f$nas.fra, pch=16, col=rgb(0,0.7,0.8,0.3),cex=2, xlab="Generations", ylab="nasuta proportion")
abline(h=0.5, lwd=2, lty=2)
hist(f$nas.fra, breaks=40, col=rgb(0,0.7,0.8,0.3))

plot(m$gen, m$alb.fra, pch=16, col=rgb(0,0.1,0.9,0.3),cex=2, xlab="Generations", ylab="albomicans proportion")
abline(h=0.5, lwd=2, lty=2)
hist(m$alb.fra, breaks=40, col=rgb(0,0.1,0.9,0.3))
plot(m$gen, m$nas.fra, pch=16, col=rgb(0,0.7,0.8,0.3),cex=2, xlab="Generations", ylab="nasuta proportion")
abline(h=0.5, lwd=2, lty=2)
hist(m$nas.fra, breaks=40, col=rgb(0,0.7,0.8,0.3))

###plot muller CD haplotype fractions for each indv, female
fo=f[order(f$alb.fra),]
f.perc=as.matrix(fo[,5:10])
f.perct=t(f.perc[rowSums(f.perc)>0,])
barplot(f.perct, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),legend =NULL)
###plot muller CD haplotype fractions for each indv, male
mo=m[order(m$alb.fra),]
m.perc=as.matrix(mo[,5:10])
m.perct=t(m.perc[rowSums(m.perc)>0,])
barplot(m.perct, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),legend =NULL)
###plot muller CD haplotype fractions for each indv, female sorted by generation
#sorted by generation
fg=f[order(f$gen),]
f.perc=as.matrix(fg[,5:10])
f.perct=t(f.perc[rowSums(f.perc)>0,])
colnames(f.perct)=fg$gen
barplot(f.perct, xlab="Generations", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),cex.name=0.3)
###plot muller CD haplotype fractions for each indv, male sorted by generation
#sorted by generation
mg=m[order(m$gen),]
m.perc=t(as.matrix(mg[,5:10]))
colnames(m.perc)=mg$gen
barplot(m.perc, xlab="Generations", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),cex.name=0.5)
###legend to decode colors
plot(x=NULL, y=NULL)
legend("topright", legend=colnames(f.perc), fill=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"))

#relationship between muller cd alb fraciton and ancestry switches, controled by admixture
plot(f$alb.fra, f$mullercd.rec*f$het, pch=16, col=rgb(0.8, 0, 0.9,0.2), cex=2, xlab="albomicans admixture proportion [Muller CD]", ylab="ancestry switches in Muller CD", ylim=c(0, 12))
points(m$alb.fra, m$mullercd.rec*m$het, pch=16, col=rgb(0, 0.2, 0.4,0.2), cex=2, xlab="albomicans admixture proportion", ylab="ancestry switches")
legend("topleft", title="sex",c("F", "M"), pch=16, col=c(rgb(0.8, 0, 0.9,0.2),rgb(0, 0.2, 0.4,0.2)), pt.cex=2)

plot(f$gen,f$mullercd.rec, pch=16, col=rgb(0.8, 0, 0.9,0.2), cex=2, xlab="Generations", ylab="Haplotype switches in Muller CD")
points(m$gen,m$mullercd.rec, pch=16, col=rgb(0, 0.2, 0.4,0.2), cex=2)
legend("topleft", title="sex",c("F", "M"), pch=16, col=c(rgb(0.8, 0, 0.9,0.2),rgb(0, 0.2, 0.4,0.2)), pt.cex=2)

muller.element=c("MullerA", "MullerB","MullerDC", "MullerE",  "MullerF" )
muller.col=c( "cornflowerblue","gold", "turquoise","coral", "lightgreen") #color for muller a, cd, b, e

#Muller A
plot(f$Generation,f$hi.a, pch=16,col=rgb.convert(muller.col[1], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller A", ylim=c(0,1))
mod=lm(f$hi.a~f$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);
points(m$Generation,m$hi.a, pch=16,col=rgb.convert(muller.col[1], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller A", ylim=c(0,1))
mod=lm(m$hi.a~m$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0.5, lwd=2, lty=2);

#Muller B
plot(f$Generation,f$hi.b, pch=16,col=rgb.convert(muller.col[2], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller B", ylim=c(0,1))
mod=lm(f$hi.b~f$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);
points(m$Generation,m$hi.b, pch=16,col=rgb.convert(muller.col[2], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller B", ylim=c(0,1))
mod=lm(m$hi.b~m$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0.5, lwd=2, lty=2);

#Muller CD
plot(f$Generation,f$hi.dc, pch=16,col=rgb.convert(muller.col[3], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller CD", ylim=c(0,1))
mod=lm(f$hi.dc~f$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);
points(m$Generation,m$hi.dc, pch=16,col=rgb.convert(muller.col[3], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller CD", ylim=c(0,1))
mod=lm(m$hi.dc~m$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0.5, lwd=2, lty=2);
#Muller E
plot(f$Generation,f$hi.e, pch=16,col=rgb.convert(muller.col[4], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller E", ylim=c(0,1))
mod=lm(f$hi.e~f$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);
points(m$Generation,m$hi.e, pch=16,col=rgb.convert(muller.col[4], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller E", ylim=c(0,1))
mod=lm(m$hi.e~m$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0.5, lwd=2, lty=2);
#Muller F
plot(f$Generation,f$hi.f, pch=16,col=rgb.convert(muller.col[5], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller F", ylim=c(0,1))
mod=lm(f$hi.f~f$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);
points(m$Generation,m$hi.f, pch=16,col=rgb.convert(muller.col[5], 0.3), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Muller F", ylim=c(0,1))
mod=lm(m$hi.f~m$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0.5, lwd=2, lty=2);

plot(f$Generation,f$HI, pch=16,col=rgb(0,0, 0, 0.1), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Genome-wide", ylim=c(0,1))
mod=lm(f$HI~f$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0.5, lwd=2, lty=2);
points(m$Generation,m$HI, pch=16,col=rgb(0,0, 0, 0.1), cex=2, xlab="Generation",ylab="alb admixture proportion", main="Genome-wide", ylim=c(0,1))
mod=lm(m$HI~m$Generation); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0.5, lwd=2, lty=2);

#TRIANGLE PLOTS FEMALE
plot(f$hi.a,f$het.a, pch=16, col=rgb.convert(muller.col[1], 0.3), cex=1.5, xlab="Admixture proportion, Muller A", ylab="Heterozygosity, Muller A", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(f$hi.b,f$het.b, pch=16, col=rgb.convert(muller.col[2], 0.3), cex=1.5, xlab="Admixture proportion, Muller B", ylab="Heterozygosity, Muller B", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(f$hi.dc,f$het.dc, pch=16, col=rgb.convert(muller.col[3], 0.3), cex=1.5, xlab="Admixture proportion, Muller CD", ylab="Heterozygosity, Muller CD", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(f$hi.e,f$het.e, pch=16, col=rgb.convert(muller.col[4], 0.3), cex=1.5, xlab="Admixture proportion, Muller E", ylab="Heterozygosity, Muller E", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(f$hi.f,f$het.f, pch=16, col=rgb.convert(muller.col[5], 0.3), cex=1.5, xlab="Admixture proportion, Muller F", ylab="Heterozygosity, Muller F", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)
plot(f$HI,f$het, pch=16, col=rgb(0,0,0, 0.1), cex=1.5, xlab="Admixture proportion", ylab="Heterozygosity", xlim=c(0,1), ylim=c(0,1))
abline(0,2);abline(2, -2);abline(h=0)

#TRIANGLE PLOTS MALE
plot(m$hi.a,m$het.a, pch=16, col=rgb.convert(muller.col[1], 0.3), cex=1.5, xlab="Admixture proportion, Muller A", ylab="Heterozygosity, Muller A", xlim=c(0,1), ylim=c(0,1));abline(0,2);abline(2, -2);abline(h=0)
plot(m$hi.b,m$het.b, pch=16, col=rgb.convert(muller.col[2], 0.3), cex=1.5, xlab="Admixture proportion, Muller B", ylab="Heterozygosity, Muller B", xlim=c(0,1), ylim=c(0,1));abline(0,2);abline(2, -2);abline(h=0)
plot(m$hi.dc,m$het.dc, pch=16, col=rgb.convert(muller.col[3], 0.3), cex=1.5, xlab="Admixture proportion, Muller CD", ylab="Heterozygosity, Muller CD", xlim=c(0,1), ylim=c(0,1));abline(0,2);abline(2, -2);abline(h=0)
plot(m$hi.e,m$het.e, pch=16, col=rgb.convert(muller.col[4], 0.3), cex=1.5, xlab="Admixture proportion, Muller E", ylab="Heterozygosity, Muller E", xlim=c(0,1), ylim=c(0,1));abline(0,2);abline(2, -2);abline(h=0)
plot(m$hi.f,m$het.f, pch=16, col=rgb.convert(muller.col[5], 0.3), cex=1.5, xlab="Admixture proportion, Muller F", ylab="Heterozygosity, Muller F", xlim=c(0,1), ylim=c(0,1));abline(0,2);abline(2, -2);abline(h=0)
plot(m$HI,m$het, pch=16, col=rgb(0,0,0, 0.1), cex=1.5, xlab="Admixture proportion", ylab="Heterozygosity", xlim=c(0,1), ylim=c(0,1));abline(0,2);abline(2, -2);abline(h=0)


plot(f$Generation, f$het, pch=16, col=rgb(0.2,0,0.1,0.2),cex=2)
mod=lm(f$het~f$Generation); summary(m)
abline(mod,lwd=2, lty=2, col=rgb(0.2,0,0.1,0.6))
points(m$Generation, m$het, pch=16, col=rgb(0,0.2,0.2,0.2),cex=2)
mod=lm(m$het~m$Generation); summary(m)
abline(mod,lwd=2, lty=2, col=rgb(0,0.2,0.2,0.6))

plot(f$HI, f$alb.fra, xlim=c(0,1), ylim=c(0,1), pch=16, col=rgb(0.6,0, 0.2, 0.2), cex=2, xlab="Genomic ancestry proportion", ylab="Muller CD ancestry proportion")
abline(0,1)
points(m$HI, m$alb.fra, xlim=c(0,1), ylim=c(0,1), pch=16, col=rgb(0,0.2, 0.6, 0.2), cex=2, xlab="Genomic ancestry proportion", ylab="Muller CD ancestry proportion")
legend("topleft", title="sex",c("F", "M"), pch=16, col=c(rgb(0.6, 0, 0.2,0.2),rgb(0, 0.2, 0.6,0.2)), pt.cex=2)
bgc.curve=function(h, p, col)
{m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm="port", lower=c(a=-1, b=-1), upper=c(a=1,b=1),control=nls.control(maxiter = 100, warnOnly=TRUE))
a=coef(m)[1]; b=coef(m)[2]; se.a=coef(summary(m))[1, "Std. Error"]; se.b=coef(summary(m))[2, "Std. Error"]
print(summary(m))
x<-data.frame(h=seq(0, 1, length.out=1000));
y<-predict(m,x)
lines(x[,1], y, col=col, lwd=2, lty=2)}
bgc.curve(f$HI, f$alb.fra, col=rgb(0.6,0, 0.2, 1))
bgc.curve(m$HI, m$alb.fra, col=rgb(0,0.2, 0.6, 1))
bgc.curve(c(f$HI,m$HI), c(f$alb.fra,m$alb.fra), col="black")
plot(m$alb.fra, m$hi.dc)

alb.intro=function(p, h)
	{print((1-h))
	y=(2*p-1)/(1-h)
	zero=which((1-h)==0)
	y[zero]=0
	return(y)} #where p = hi, and h=het
f$alb.intro=alb.intro(f$HI, f$het)
plot(f$HI,f$alb.intro)

f$alb.intro.cd=alb.intro(f$hi.dc, f$het.dc)
f$alb.intro.a=alb.intro(f$hi.a, f$het.a)
f$alb.intro.b=alb.intro(f$hi.b, f$het.b)
f$alb.intro.e=alb.intro(f$hi.e, f$het.e)
f$alb.intro.f=alb.intro(f$hi.f, f$het.f)
plot(f$hi.dc, f$mullercd.rec*f$het)
m$alb.intro.cd=alb.intro(m$hi.dc, m$het.dc)
m$alb.intro.a=alb.intro(m$hi.a, m$het.a)
m$alb.intro.b=alb.intro(m$hi.b, m$het.b)
m$alb.intro.e=alb.intro(m$hi.e, m$het.e)
m$alb.intro.f=alb.intro(m$hi.f, m$het.f)

plot(f$Generation,f$alb.intro.a, pch=16,col=rgb.convert(muller.col[1], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller A", ylim=c(-1,1))
mod=lm(f$alb.intro.a~f$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);
points(m$Generation,m$alb.intro.a, pch=16,col=rgb.convert(muller.col[1], 0.3), cex=2)
mod=lm(m$alb.intro.a~m$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0, lwd=2, lty=2);


plot(f$Generation,f$alb.intro.b, pch=16,col=rgb.convert(muller.col[2], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller B", ylim=c(-1,1))
mod=lm(f$alb.intro.b~f$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);
points(m$Generation,m$alb.intro.b, pch=16,col=rgb.convert(muller.col[2], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller B", ylim=c(-1,1))
mod=lm(m$alb.intro.b~m$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0, lwd=2, lty=2);

plot(f$Generation,f$alb.intro.cd, pch=16,col=rgb.convert(muller.col[3], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller CD", ylim=c(-1,1))
mod=lm(f$alb.intro.cd~f$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);
points(m$Generation,m$alb.intro.cd, pch=16,col=rgb.convert(muller.col[3], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller CD", ylim=c(-1,1))
mod=lm(m$alb.intro.cd~m$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0, lwd=2, lty=2);

plot(f$Generation,f$alb.intro.e, pch=16,col=rgb.convert(muller.col[4], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller E", ylim=c(-1,1))
mod=lm(f$alb.intro.e~f$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);
points(m$Generation,m$alb.intro.e, pch=16,col=rgb.convert(muller.col[4], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller E", ylim=c(-1,1))
mod=lm(m$alb.intro.e~m$Generation+0); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0, lwd=2, lty=2);

plot(f$Generation,f$alb.intro.f, pch=16,col=rgb.convert(muller.col[5], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller F", ylim=c(-1,1))
mod=lm(f$alb.intro.f~f$Generation-1); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);
points(m$Generation,m$alb.intro.f, pch=16,col=rgb.convert(muller.col[5], 0.3), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Muller F", ylim=c(-1,1))
mod=lm(m$alb.intro.f~m$Generation-1); summary(m)
abline(mod, lwd=4, lty=1, col="blue"); abline(h=0, lwd=2, lty=2);

plot(f$Generation,f$alb.intro, pch=16,col=rgb(0,0, 0, 0.1), cex=2, xlab="Generation",ylab="alb-biased introgression", main="Genome-wide", ylim=c(-1,1))
mod=lm(f$alb.intro~0+f$G eneration); summary(m)
abline(mod, lwd=4, lty=1, col="red"); abline(h=0, lwd=2, lty=2);

##albomicans introgression
vioplot(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)

vioplot(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)

#recombination rates
vioplot(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)
vioplot(m$mullera.rec,m$mullerb.rec,m$mullercd.rec,m$mullere.rec,m$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(m$mullera.rec,m$mullerb.rec,m$mullercd.rec,m$mullere.rec,m$mullerf.rec ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(m$mullera.rec,m$mullerb.rec,m$mullercd.rec,m$mullere.rec,m$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)

#recombination rates corrected by admixture, female
vioplot(f$mullera.rec*f$het,f$mullerb.rec*f$het,f$mullercd.rec*f$het,f$mullere.rec*f$het,f$mullerf.rec*f$het ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches (admixture-corrected)", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(f$mullera.rec*f$het,f$mullerb.rec*f$het,f$mullercd.rec*f$het,f$mullere.rec*f$het,f$mullerf.rec*f$het ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(f$mullera.rec*f$het,f$mullerb.rec*f$het,f$mullercd.rec*f$het,f$mullere.rec*f$het,f$mullerf.rec*f$het ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)
#recombination rates corrected by admixture, male
vioplot(m$mullera.rec*m$het,m$mullerb.rec*m$het,m$mullercd.rec*m$het,m$mullere.rec*m$het,m$mullerf.rec*m$het ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches (admixture-corrected)", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(m$mullera.rec*m$het,m$mullerb.rec*m$het,m$mullercd.rec*m$het,m$mullere.rec*m$het,m$mullerf.rec*m$het ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(m$mullera.rec*m$het,m$mullerb.rec*m$het,m$mullercd.rec*m$het,m$mullere.rec*m$het,m$mullerf.rec*m$het ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)



#number of sites
#mullerA, mullerB, mullerCD, mullerE, muller F informative sites
info.sites=c(128670, 120223, 188404, 126724, 3036)
info.fr=c(128670, 120223, 188404, 126724, 3036)/sum(c(128670, 120223, 188404, 126724, 3036))