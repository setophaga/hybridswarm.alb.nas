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

f$hi.dc.corrected=f$alb.fra 
f$het.dc.corrected=(f$albX.nas+f$albY.nas)
plot(f$hi.dc.corrected, f$het.dc.corrected)
m$hi.dc.corrected=m$alb.fra 
m$het.dc.corrected=(m$albX.nas+m$albY.nas)
plot(m$hi.dc.corrected, m$het.dc.corrected)

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

correct.rec=function(f){f$mullera.rec=f$mullera.rec-1; f$mullerb.rec=f$mullerb.rec-1; f$mullercd.rec=f$mullercd.rec-1; f$mullere.rec=f$mullere.rec-1; f$mullerf.rec=f$mullerf.rec-1; return(f)}
f=correct.rec(f); m=correct.rec(m)
alb.intro=function(p, h)
	{print((1-h))
	y=(2*p-1)/(1-h)
	zero=which((1-h)==0)
	y[zero]=0
	return(y)} 
nocd.het.hi=function(f){
f$hi.noCD=rowMeans(cbind(f$hi.a, f$hi.b, f$hi.e, f$hi.f));f$het.noCD=rowMeans(cbind(f$het.a, f$het.b, f$het.e, f$het.f)); 
f$alb.intro.noCD=alb.intro(f$hi.noCD, f$het.noCD);return(f)}
f=nocd.het.hi(f);m=nocd.het.hi(m)
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
barplot(f.perct, xlab="Individuals", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"), cex.axis=0.2,legend =NULL)
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
table(fg$gen)
barplot(f.perct, xlab="Generations", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),cex.name=0.3, cex.axis=0.2)
###plot muller CD haplotypefractions for each indv, male sorted by generation
#sorted by generation
mg=m[order(m$gen),]
m.perc=t(as.matrix(mg[,5:10]))
colnames(m.perc)=mg$gen
barplot(m.perc, xlab="Generations", col=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"),cex.name=0.5)
###legend to decode colors
plot(x=NULL, y=NULL)
legend("topright", legend=c("Neo-X, Neo-X", "Neo-X, Neo-Y", "Neo-X, nas", "NeoY, NeoY", "NeoY, nas", "nas, nas"), fill=c("coral4","deepskyblue", "gold", "chartreuse3", "lightslateblue", "aquamarine2"))

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


#TRIANGLE PLOTS FEMALE and males together
################################################################
plot(c(f$hi.a, m$hi.a),c(f$het.a, m$het.a), pch=21, bg=rgb.convert(muller.col[1], 0.3),col=rgb.convert(muller.col[1], 1), cex=1.5, xlab="Admixture proportion", ylab="Heterozygosity", xlim=c(0,1), ylim=c(0,1), lwd=1)
abline(0,2, lwd=1, lty=2);abline(2, -2, lwd=1, lty=2);abline(h=0, lwd=1, lty=2)

plot(c(f$hi.b, m$hi.b),c(f$het.b,m$het.b), pch=21, bg=rgb.convert(muller.col[2], 0.3),col=rgb.convert(muller.col[2], 1), cex=1.5, xlab="Admixture proportion, Muller B", ylab="Heterozygosity, Muller B", xlim=c(0,1), ylim=c(0,1))
abline(0,2, lwd=1, lty=2);abline(2, -2, lwd=1, lty=2);abline(h=0, lwd=1, lty=2)

plot(c(f$hi.dc,m$hi.dc), c(f$het.dc, m$het.dc),pch=21, col=rgb.convert(muller.col[3], 1), bg=rgb.convert(muller.col[3], 0.3),cex=1.5, xlab="Admixture proportion, Muller CD", ylab="Heterozygosity, Muller CD", xlim=c(0,1), ylim=c(0,1))
abline(0,2, lwd=1, lty=2);abline(2, -2, lwd=1, lty=2);abline(h=0, lwd=1, lty=2)

plot(c(f$hi.e,m$hi.e),c(f$het.e, m$het.e),pch=21, col=rgb.convert(muller.col[4], 1),bg=rgb.convert(muller.col[4], 0.3), cex=1.5, xlab="Admixture proportion", ylab="Heterozygosity", xlim=c(0,1), ylim=c(0,1))
abline(0,2, lwd=1, lty=2);abline(2, -2, lwd=1, lty=2);abline(h=0, lwd=1, lty=2)

plot(c(f$hi.f,m$hi.f), c(f$het.f,m$het.f), pch=21, col=rgb.convert(muller.col[5], 1), bg=rgb.convert(muller.col[5], 0.3),cex=1.5, xlab="Admixture proportion", ylab="", xlim=c(0,1), ylim=c(0,1))
abline(0,2, lwd=1, lty=2);abline(2, -2, lwd=1, lty=2);abline(h=0, lwd=1, lty=2)

plot(c(f$HI,m$HI),c(f$het,m$het), pch=21, col=rgb(0,0,0, 1),bg=rgb(0,0,0, 0.1), cex=1.5, xlab="Admixture proportion", ylab="", xlim=c(0,1), ylim=c(0,1))
abline(0,2, lwd=1, lty=2);abline(2, -2, lwd=1, lty=2);abline(h=0, lwd=1, lty=2); abline(v=0.5, lwd=1, lty=2)
################################################################

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


plot(f$Generation, fhet, pch=16, col=rgb(0.2,0,0.1,0.2),cex=2)
mod=lm(f$het~f$Generation); summary(m)
abline(mod,lwd=2, lty=2, col=rgb(0.2,0,0.1,0.6))
points(m$Generation, m$het, pch=16, col=rgb(0,0.2,0.2,0.2),cex=2)
mod=lm(m$het~m$Generation); summary(m)
abline(mod,lwd=2, lty=2, col=rgb(0,0.2,0.2,0.6))


#print out predominant haplotype in muller CD
mode.plot=function(f){
	cols=c(which(colnames(f)=="albX.albX"):which(colnames(f)=="nas.nas")); cd.hap={}
	for(i in 1:nrow(f))
		{if(length(which(f[i,cols]>0.5))>0)
			{cd.hap=c(cd.hap,which(f[i,cols]>0.5))}
		else{cd.hap=c(cd.hap,NA); print(f[i,cols])}}
	f$cd.hap=cd.hap
		return(f)}
		
	f=mode.plot(f)
	m=mode.plot(m)

col.mode<- c( "coral4", "deepskyblue3","gold","green", "lightslateblue","aquamarine2")
mode.hap=c("neoX,neoX", "neoX,neoY","nas,neoX","neoY, neoY", "nas,neoY", "nas,nas") #

############. ############. ############. BGC 
############. ############. ############. BGC 
plot(f$HI, f$alb.fra, xlim=c(0,1), ylim=c(0,1), pch=21,  bg=rgb.convert(col.mode[as.numeric(f$cd.hap)], 0.6), cex=1.3, xlab="Genomic ancestry proportion", ylab="Muller CD ancestry proportion")
abline(0,1, lty=2, lwd=2)
points(m$HI, m$alb.fra, xlim=c(0,1), ylim=c(0,1), pch=21, bg=rgb.convert(col.mode[as.numeric(m$cd.hap)], 0.6), cex=1.3, xlab="Genomic ancestry proportion", ylab="Muller CD ancestry proportion")
legend("topleft", title="Muller_CD haplotype",mode.hap[-4], pch=21, pt.bg=rgb.convert(col.mode[-4], 0.6), pt.cex=1.3, cex=0.7)
bgc.curve=function(h, p, col)
{m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=1),  control=nls.control(maxiter = 100, warnOnly=TRUE)) #algorithm="plinear", or "port"
a=coef(m)[1]; b=coef(m)[2]; se.a=coef(summary(m))[1, "Std. Error"]; se.b=coef(summary(m))[2, "Std. Error"]
print(summary(m))
x<-data.frame(h=seq(0, 1, length.out=1000));
y<-predict(m,x)
lines(x[,1], y, col=col, lwd=2, lty=1)}
#bgc.curve(f$HI[-which(f$gen==0)], f$alb.fra[-which(f$gen==0)], col=rgb(0.6,0, 0.2, 1))
#bgc.curve(m$HI[-which(m$gen==0)], m$alb.fra[-which(m$gen==0)], col=rgb(0,0.2, 0.6, 1))
bgc.curve(c(f$HI,m$HI), c(f$alb.fra,m$alb.fra), col="black")

plot(m$alb.fra, m$hi.dc)

alb.intro=function(p, h)
	{print((1-h))
	y=(2*p-1)/(1-h)
	zero=which((1-h)==0)
	y[zero]=0
	return(y)} #where p = hi, and h=het
f$alb.intro=alb.intro(f$HI, f$het)
m$alb.intro=alb.intro(m$HI, m$het)
plot(f$HI,f$alb.intro)

f$alb.intro.cd=alb.intro(f$hi.dc.corrected, f$het.dc.corrected)
f$alb.intro.a=alb.intro(f$hi.a, f$het.a)
f$alb.intro.b=alb.intro(f$hi.b, f$het.b)
f$alb.intro.e=alb.intro(f$hi.e, f$het.e)
f$alb.intro.f=alb.intro(f$hi.f, f$het.f)
plot(f$hi.dc, f$mullercd.rec*f$het)
m$alb.intro.cd=alb.intro(m$hi.dc.corrected, m$het.dc.corrected)
m$alb.intro.a=alb.intro(m$hi.a, m$het.a)
m$alb.intro.b=alb.intro(m$hi.b, m$het.b)
m$alb.intro.e=alb.intro(m$hi.e, m$het.e)
m$alb.intro.f=alb.intro(m$hi.f, m$het.f)

f=f[-which(f$gen==0),]
m=m[-which(m$gen==0),]
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
intro.f=data.frame(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f );colnames(intro.f)=c("alb.intro.a","alb.intro.b","alb.intro.cd","alb.intro.e","alb.intro.f")
intro.m=data.frame(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f);colnames(intro.m)=c("alb.intro.a","alb.intro.b","alb.intro.cd","alb.intro.e","alb.intro.f")
intro=rbind(intro.f, intro.m)
vioplot(intro$alb.intro.a,intro$alb.intro.b,intro$alb.intro.cd,intro$alb.intro.e,intro$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="alb Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(intro$alb.intro.a,intro$alb.intro.b,intro$alb.intro.cd,intro$alb.intro.e,intro$alb.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(intro$alb.intro.a,intro$alb.intro.b,intro$alb.intro.cd,intro$alb.intro.e,intro$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="alb Introgression", xlab="", lwd=2, names= muller.element, add=T)
alb.intro.y=c(intro$alb.intro.a,intro$alb.intro.b,intro$alb.intro.cd,intro$alb.intro.e,intro$alb.intro.f)
muller.x=c(rep("MullerA", nrow(intro)),rep("MullerB", nrow(intro)),rep("MullerCD", nrow(intro)),rep("MullerE", nrow(intro)),rep("MullerF", nrow(intro)))
pairwise.t.test(alb.intro.y,muller.x, p.adj="fdr")

vioplot(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(f$alb.intro.a,f$alb.intro.b,f$alb.intro.cd,f$alb.intro.e,f$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)

vioplot(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(m$alb.intro.a,m$alb.intro.b,m$alb.intro.cd,m$alb.intro.e,m$alb.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)


#recombination rates
rec.f=data.frame(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ); colnames(rec.f)=c("mullera.rec","mullerb.rec", "mullercd.rec","mullere.rec","mullerf.rec")
rec.m=data.frame(m$mullera.rec,m$mullerb.rec,m$mullercd.rec,m$mullere.rec,m$mullerf.rec ); colnames(rec.m)=c("mullera.rec","mullerb.rec", "mullercd.rec","mullere.rec","mullerf.rec")

rec=rbind(rec.f,rec.m)
sum=128670+120223+188404+126724+3036
f$lambda=sqrt((f$HI-0.5)^2+f$het^2); m$lambda=sqrt((m$HI-0.5)^2+m$het^2); 
rec$lambda=c(f$labmda, m$labmda)
#both female and male
vioplot(rec$mullera.rec/128670,rec$mullerb.rec/120223,rec$mullercd.rec/188404,rec$mullere.rec/126724,col=rgb.convert(muller.col[-5], 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="", xlab="", lwd=2, names= muller.element[-5])
stripchart(data.frame(rec$mullera.rec/128670,rec$mullerb.rec/120223,rec$mullercd.rec/188404,rec$mullere.rec/126724), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col[-5], 1))  
vioplot(rec$mullera.rec/128670,rec$mullerb.rec/120223,rec$mullercd.rec/188404,rec$mullere.rec/126724 ,col=rgb.convert(muller.col[-5], 0.1) ,xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element[-5], add=T)
rec.corrected=c(rec$mullera.rec/128670,rec$mullerb.rec/120223,rec$mullercd.rec/188404,rec$mullere.rec/126724)
mul.reps=c(rep("mullerA", nrow(rec)),rep("mullerB", nrow(rec)),rep("mullerCD", nrow(rec)),rep("mullerE", nrow(rec)))
summary(aov(rec.corrected~mul.reps))
pairwise.t.test(rec.corrected, mul.reps)
##female only
vioplot(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Ancestry switches", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(f$mullera.rec,f$mullerb.rec,f$mullercd.rec,f$mullere.rec,f$mullerf.rec ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= muller.element, add=T)
##male only
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

#genomewide introgression ~muller cd happlotype female and male
alb.intro=c(f$alb.intro.noCD, m$alb.intro.noCD)
cd.hap=c(f$cd.hap, m$cd.hap)
vioplot(alb.intro~cd.hap,col=rgb.convert(col.mode[as.numeric(sort(unique(cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(cd.hap)))],horizontal = F, las = 1, cex=1.2, ylab="alb introgression (admixture-corrected)", xlab="Muller CD haplotypes", lwd=2, names= mode.hap[as.numeric(sort(unique(cd.hap)))])
stripchart( alb.intro~cd.hap, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(cd.hap)))], 1))  
vioplot(alb.intro~cd.hap, col=rgb.convert(col.mode[as.numeric(sort(unique(cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(cd.hap)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(cd.hap)))], add=T)
print(summary(aov((alb.intro)~cd.hap)))
pairwise.t.test((alb.intro), cd.hap, p.adj = "fdr")
#albintro cd only 
alb.intro.cd=c(f$alb.intro.cd, m$alb.intro.cd)
cd.hap=c(f$cd.hap, m$cd.hap)
vioplot(alb.intro.cd~cd.hap,col=rgb.convert(col.mode[as.numeric(sort(unique(cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(cd.hap)))],horizontal = F, las = 1, cex=1.2, ylab="alb introgression (admixture-corrected)", xlab="Muller CD haplotypes", lwd=2, names= mode.hap[as.numeric(sort(unique(cd.hap)))])
stripchart( alb.intro.cd~cd.hap, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(cd.hap)))], 1))  
vioplot(alb.intro.cd~cd.hap, col=rgb.convert(col.mode[as.numeric(sort(unique(cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(cd.hap)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(cd.hap)))], add=T)
print(summary(aov((alb.intro.cd)~cd.hap)))
pairwise.t.test((alb.intro.cd), cd.hap, p.adj = "fdr")

#genomewide introgression ~muller cd happlotype female and male seperate
plot.cdhap.gwalbintro=function(f)
{vioplot(f$alb.intro~f$cd.hap,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(f$cd.hap)))],horizontal = F, las = 1, cex=1.2, ylab="alb introgression (admixture-corrected)", xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(f$cd.hap)))])
stripchart( alb.intro~cd.hap, data=f,vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 1))  
vioplot(alb.intro~cd.hap,data=f,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(f$cd.hap)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(f$cd.hap)))], add=T)
print(summary(aov((f$alb.intro)~f$cd.hap)))
pairwise.t.test((f$alb.intro), f$cd.hap, p.adj = "fdr")
}
plot.cdhap.gwalbintro(f[-which(f$gen==0),])
plot.cdhap.gwalbintro(m[-which(m$gen==0),]) 


plot.cdhap.gwancswitch=function(f)
{f$rec=f$mullera.rec-1+f$mullerb.rec-1+f$mullercd.rec-1  +f$mullere.rec-1 +f$mullerf.rec -1#
y=(f$rec*sqrt((f$het-0.5)^2+f$HI^2)); x=f$cd.hap 
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(f$cd.hap)))],horizontal = F, las = 1, cex=1.2, ylab="Haplotype switches (admixture-corrected)", xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(f$cd.hap)))])
stripchart( y~x, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 1))  
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(f$cd.hap)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(f$cd.hap)))], add=T)
print(summary(aov(y~x)))
pairwise.t.test(y,x, p.adj="fdr")
}
plot.cdhap.gwancswitch(f[-which(f$gen==0),]); 
plot.cdhap.gwancswitch(m[-which(m$gen==0),]); 

#Plot genomwise ancestry switches with males and females together
f$rec=f$mullera.rec+f$mullerb.rec +f$mullercd.rec +f$mullere.rec+f$mullerf.rec
m$rec=m$mullera.rec+m$mullerb.rec +m$mullercd.rec +m$mullere.rec +m$mullerf.rec
f$rec.nocd=f$mullera.rec+f$mullerb.rec  +f$mullere.rec+f$mullerf.rec
m$rec.nocd=m$mullera.rec+m$mullerb.rec  +m$mullere.rec +m$mullerf.rec
rec.gw=c(f$rec, m$rec)
rec.nocd=c(f$rec.nocd, m$rec.nocd)
x=c(f$cd.hap,m$cd.hap)
lambda=c(sqrt((f$het-0.5)^2+f$HI^2), sqrt((m$het-0.5)^2+m$HI^2))
y=rec.gw*lambda
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(x)))],horizontal = F, las = 1, cex=1.2, ylab="Haplotype switches (admixture-corrected)", xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(x)))])
stripchart( y~x, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 1))  
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(x)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(x)))], add=T)
print(summary(aov(y~x)))
pairwise.t.test(y,x, p.adj="fdr")

###switches only CD
rec.cd=c(f$mullercd.rec, m$mullercd.rec)
y=rec.cd*lambda
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(x)))],horizontal = F, las = 1, cex=1.2, ylab="Haplotype switches (admixture-corrected)", xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(x)))])
stripchart( y~x, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 1))  
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(x)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(x)))], add=T)
print(summary(aov(y~x)))
pairwise.t.test(y,x, p.adj="fdr")

###switches without CD
y=rec.nocd*lambda
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(x)))],horizontal = F, las = 1, cex=1.2, ylab="Haplotype switches (admixture-corrected)", xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(x)))], cex.axis=1.3)
stripchart( y~x, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 1))  
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(x)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(x)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(x)))], add=T)
print(summary(aov(y~x)))
pairwise.t.test(y,x, p.adj="fdr")

plot.cdhap.rawancswitch=function(f)
{f$rec=f$mullera.rec+f$mullerb.rec +f$mullercd.rec +f$mullere.rec +f$mullerf.rec
y=(f$rec); x=f$cd.hap
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(f$cd.hap)))],horizontal = F, las = 1, cex=1.2, ylab="Genome-wide ancestry switches", xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(f$cd.hap)))])
stripchart( y~x, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 1))  
vioplot(y~x,col=rgb.convert(col.mode[as.numeric(sort(unique(f$cd.hap)))], 0.1), xaxt="n",border = col.mode[as.numeric(sort(unique(f$cd.hap)))],horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap[as.numeric(sort(unique(f$cd.hap)))], add=T)
summary(aov(y~x))
pairwise.t.test(y,x, p.adj="fdr")
}

plot.cdhap.rawancswitch(f[-which(f$gen==0),]); 
plot.cdhap.rawancswitch(m[-which(m$gen==0),]); 

