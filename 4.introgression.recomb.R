library(RColorBrewer)
library(vioplot)
rgb.convert=function(col.vect,fr)
{cols={};for(i in col.vect){cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}

d=read.csv("~/Desktop/hybridswarm/1.pipeline/posterior.alb03.nas00.all.July22th/hi.recomb.csv")
plot(d$info.sites.cd, d$info.sites.neosex)
hist(d$info.sites.neosex, breaks=50)
d=d[-which(d$info.sites.neosex<(max(d$info.sites.neosex, na.rm=T)*0.5)),]
d$nnn[which(d$mode==0.75)]
d=d[-which(d$Species=="albomicans"),] #remove parentals
d=d[-which(d$Species=="nasuta"),]
 
d$HI=1-d$HI #CAREFUL!! ONLY run it once, so that 0 for nasuta, 1 for albomicans
d$hi.a=1-d$hi.a; d$hi.b=1-d$hi.b;d$hi.e=1-d$hi.e;d$hi.dc=1-d$hi.dc;d$hi.f=1-d$hi.f
asy.intro=(2*d$HI-1)/(1-d$het)
asy.intro.a=(2*d$hi.a-1)/(1-d$het.a)
asy.intro.b=(2*d$hi.b-1)/(1-d$het.b)
asy.intro.e=(2*d$hi.e-1)/(1-d$het.e)
asy.intro.dc=(2*d$hi.dc-1)/(1-d$het.dc)
asy.intro.f=(2*d$hi.f-1)/(1-d$het.f)
col.mode<- c("aquamarine2","darkslateblue", "gold", "deepskyblue3", "coral4")
mode.hap=c("nas,nas","nas,neoY","nas,neoX","neoX,neoY", "neoX,neoX") #

cd=d[which(d$info.sites.cd>0),]
plot( log((d$recomb.mullercd*d$het)/d$info.sites.neosex, base=10),d$HI, ylab="Ancestry Proprotion (alb)", xlab=expression(paste(log[10], "(admixture-corrected recombination rate Muller CD)")),col=rgb.convert(col.mode[as.factor(d$mode)], 0.8), pch=16, cex=1.5) #divided by het
m=lm(d$HI~log((d$recomb.mullercd*d$het)/d$info.sites.neosex, base=10))
abline(m, lwd=2, lty=2)
col.mode<- c("aquamarine2","lightslateblue", "gold", "deepskyblue3", "coral4")
triplot=function(hi, het, main)
{plot(hi, het, xlim=c(0,1), ylim=c(0,1),bg=rgb.convert(col.mode[as.factor(d$mode)], 0.5), pch=21, cex=1.5, xlab="Ancestry", ylab="Heterozygosity", main=main)
#legend("topright",title="Muller_CD mode", mode.hap, pch=21, pt.bg=rgb.convert(col.mode, 0.5), cex=0.6, pt.cex=1.5)
abline(0,2);abline(2,-2)}

triplot(d$hi.dc, d$het.dc, "Muller_DC")
triplot(d$hi.a, d$het.a, "Muller_A")
triplot(d$hi.b, d$het.b, "Muller_B")
triplot(d$hi.e, d$het.e, "Muller_E")
triplot(d$hi.f, d$het.f, "Muller_F")
triplot(d$HI, d$het, "Whole genome")

plot(x=NULL, y=NULL)
legend("topright",mode.hap, fill=col.mode, cex=1.5)
d$Species=as.character(d$Species)

#recombination frequency in different muller cd ancestry background
recomb.cd=log(,base=10)
vioplot(((d$recomb.mullercd*d$het)/d$info.sites.neosex)~d$mode,col=rgb.convert(col.mode, 0.1), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Admixture-corrected Recombination", xlab="", lwd=2, names= mode.hap)
stripchart(((d$recomb.mullercd*d$het)/d$info.sites.neosex)~d$mode, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode, 0.8))  
vioplot((d$recomb.mullercd*d$het)/d$info.sites.neosex~d$mode,col=rgb.convert(col.mode, 0.2), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2,  xlab="", lwd=2, names= mode.hap, add=T)

d$recomb=sum(d$recomb.mullera, d$recomb.mullerb, d$recomb.mullercd, d$recomb.mullere)
d$info.sites.atoe=sum(d$info.sites.a, d$info.sites.b, d$info.sites.cd,d$info.sites.e)

plot((d$recomb*d$het)/d$info.sites.atoe,asy.intro, col=col.neosex.sp[as.factor(d$mode)], pch=16)
vioplot(asy.intro~as.character(d$mode),col=rgb.convert(col.mode, 0.1), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Introgression Asymmetry", xlab="", lwd=2, names= mode.hap)
stripchart(asy.intro~d$mode, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode, 0.8))  
vioplot(asy.intro~as.character(d$mode),col=rgb.convert(col.mode, 0.2), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Introgression Asymmetry", xlab="", lwd=2, names= mode.hap, add=T)

vioplot(asy.intro.dc~as.character(d$mode),col=rgb.convert(col.mode, 0.1), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Introgression Asymmetry (Muller CD)", xlab="", lwd=2, names= mode.hap)
stripchart(asy.intro.dc~d$mode, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode, 0.8))  
vioplot(asy.intro.dc~as.character(d$mode),col=rgb.convert(col.mode, 0.2), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Introgression Asymmetry", xlab="", lwd=2, names= mode.hap, add=T)


muller.element=c("MullerA", "MullerDC", "MullerB","MullerE",  "MullerF" )
muller.col=c( "cornflowerblue","gold", "turquoise","coral", "lightgreen") #color for muller a, cd, b, e
vioplot(asy.intro.a,asy.intro.dc,asy.intro.b,asy.intro.e,asy.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(asy.intro.a,asy.intro.dc,asy.intro.b,asy.intro.e,asy.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(asy.intro.a,asy.intro.dc,asy.intro.b,asy.intro.e,asy.intro.f ,col=rgb.convert(muller.col, 0.2), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)

#Figure chromosome recombination rate
muller.element.nof=c("MullerA", "MullerDC", "MullerB","MullerE");muller.col.nof=c( "cornflowerblue","gold", "turquoise","coral")
vioplot((d$recomb.mullera*d$het/d$info.sites.a),(d$recomb.mullercd*d$het/d$info.sites.cd),(d$recomb.mullerb*d$het/d$info.sites.b),(d$recomb.mullere*d$het/d$info.sites.e), col=rgb.convert(muller.col.nof, 0.1), xaxt="n",border = muller.col.nof,horizontal = F, las = 1, cex=1.2, ylab="Recombination (admiture-corrected)", xlab="", lwd=2, names= muller.element.nof)
#,(d$recomb.mullerf*d$het/(d$info.sites.f+1)),
stripchart(data.frame((d$recomb.mullera*d$het/d$info.sites.a),(d$recomb.mullercd*d$het/d$info.sites.cd),(d$recomb.mullerb*d$het/d$info.sites.b),(d$recomb.mullere*d$het/d$info.sites.e)), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col.nof, 1))  
vioplot((d$recomb.mullera*d$het/d$info.sites.a),(d$recomb.mullercd*d$het/d$info.sites.cd),(d$recomb.mullerb*d$het/d$info.sites.b),(d$recomb.mullere*d$het/d$info.sites.e), xaxt="n",border = muller.col.nof,horizontal = F, las = 1, cex=1.2, ylab="", xlab="", lwd=2, names=muller.element.nof, add=T,col=rgb.convert(muller.col.nof, 0.4))


plot((d$recomb.mullercd*d$het/d$info.sites.neosex),asy.intro.dc, col=rgb.convert(col.mode[as.factor(d$mode)], 0.5), pch=16, ylab="Asymmetrical Introgression (Muller CD)", xlab=expression(paste(log[10],"(admixture-corrected recombination frequency Muller CD)")), cex=2)

x=log(d$recomb.mullercd*d$het/(d$info.sites.cd),base=10)
plot(x,asy.intro, col=rgb.convert(col.mode[as.factor(d$mode)], 0.7), pch=16, ylab="Asymmetrical Introgression", xlab=expression("log(Muller CD admixture-corrected recombination frequency)"), cex=2)
m=lm(asy.intro~x);summary(m)
abline(m, lwd=2, lty=2)

plot(log((d$recomb.mullercd*d$het/(d$info.sites.cd)), base=10),asy.intro.dc, col=rgb.convert(col.mode[as.factor(d$mode)], 0.8), pch=16, ylab="Asymmetrical Introgression (Muller CD)", xlab=expression(paste(log[10],"(admixture-corrected recombination frequency)")), cex=2)
x1=log((d$rec*d$het/(d$info.sites)), base=10); y=asy.intro.dc
lms=loess(y~x1, span=1); plx=predict(lms,x1)
lines(x1,plx, col="grey")


#find the weirdos
female.weird=as.character(d$names[c(intersect(which(d$sex.g=="F"), which(d$mode==0.75)),intersect(which(d$sex.g=="F"), which(d$mode==0.25)))])
male.weird=as.character(d$names[intersect(which(d$sex.g=="M"), which(d$mode==1))])
ss=read.csv("~/Desktop/hybridswarm/old.lib.mullers.mapreads.csv")
nasex=as.character(ss$id[which(is.na(ss$sex.g))]) 
ha=intersect(c(male.weird, female.weird), nasex)
rownames(d)=d$names; rownaems(ss)=as.character(ss$id)
d[ha,]; ss[which(ss$id==ha[1]),]
d[c(male.weird, female.weird),]
conflict=c(intersect(which(d$Gender=="male"), which(d$sex.g=="F")),intersect(which(d$Gender=="female"), which(d$sex.g=="M")))
d[conflict,]
mwhat=intersect(which(d$mode==1),which(d$sex.g=="M"))
fwhat=c(intersect(which(d$sex.g=="F"), which(d$mode==0.75)),intersect(which(d$sex.g=="F"), which(d$mode==0.25)))
plot(d$sex.auto.fr,d$total.mreads, pch=16, col="grey")
points(d$sex.auto.fr[mwhat],d$total.mreads[mwhat],bg="forestgreen", pch=21, cex=1.3)
points(d$sex.auto.fr[fwhat],d$total.mreads[fwhat], bg="gold", pch=21, cex=1.3)
abline(v=c(0.7,0.9), lwd=2, lty=2)


