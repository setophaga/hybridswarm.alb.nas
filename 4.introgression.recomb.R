library(RColorBrewer)
library(vioplot)

d=read.csv("~/Desktop/hybridswarm/1.pipeline/posterior.alb03.nas00.all.July22th/hi.recomb.csv")
max(d$info.sites, na.rm=T)
hist(d$info.sites, breaks=50)
d=d[-which(d$info.sites<159848*0.3),]
d=d[-which(is.na(d$Generation)==1),]
d$HI=1-d$HI #CAREFUL!! ONLY run it once, so that 0 for nasuta, 1 for albomicans
d$hi.a=1-d$hi.a; d$hi.b=1-d$hi.b;d$hi.e=1-d$hi.e;d$hi.dc=1-d$hi.dc;d$hi.f=1-d$hi.f
asy.intro=(2*d$HI-1)/(1-d$het)
asy.intro.a=(2*d$hi.a-1)/(1-d$het.a)
asy.intro.b=(2*d$hi.b-1)/(1-d$het.b)
asy.intro.e=(2*d$hi.e-1)/(1-d$het.e)
asy.intro.dc=(2*d$hi.dc-1)/(1-d$het.dc)
asy.intro.f=(2*d$hi.f-1)/(1-d$het.f)

plot(d$HI, log(d$rec/(d$info.sites*d$het), base=10)) #divided by het
boxplot((d$rec/(d$info.sites*d$het))~d$mode)

col.mode<- c("aquamarine2","darkolivegreen3", "gold", "deepskyblue3", "darkslateblue")
triplot=function(hi, het, main)
{plot(hi, het, xlim=c(0,1), ylim=c(0,1),col=rgb.convert(col.mode[as.factor(d$mode)], 0.5), pch=16, cex=1.5, xlab="Ancestry", ylab="Heterozygosity", main=main)
abline(0,2);abline(2,-2)}

triplot(d$hi.dc, d$het.dc, "Muller_DC")
triplot(d$hi.a, d$het.a, "Muller_A")
triplot(d$hi.b, d$het.b, "Muller_B")
triplot(d$hi.e, d$het.e, "Muller_E")
triplot(d$hi.f, d$het.f, "Muller_F")


mode.hap=c("nas,nas","nas,neoY","nas,neoX","neoX,neoY", "neoX,neoX")
plot((d$rec/(d$info.sites*d$het)),asy.intro, col=col.neosex.sp[as.factor(d$mode)], pch=16)

vioplot(asy.intro~d$mode,col=rgb.convert(col.mode, 0.1), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Introgression Asymmetry", xlab="", lwd=2, names= mode.hap)
stripchart(asy.intro~d$mode, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(col.mode, 0.8))  
vioplot(asy.intro~d$mode,col=rgb.convert(col.mode, 0.2), xaxt="n",border = col.mode,horizontal = F, las = 1, cex=1.2, ylab="Introgression Asymmetry", xlab="", lwd=2, names= mode.hap, add=T)



muller.element=c("MullerA", "MullerDC", "MullerB","MullerE",  "MullerF" )
muller.col=c( "cornflowerblue","gold", "turquoise","coral", "lightgreen") #color for muller a, cd, b, e
rgb.convert=function(col.vect,fr)
{cols={};for(i in col.vect){cl=col2rgb(i)/255; cols=c(cols, rgb(cl[1,1],cl[2,1],cl[3,1], fr))}; return(cols)}
vioplot(asy.intro.a,asy.intro.dc,asy.intro.b,asy.intro.e,asy.intro.f ,col=rgb.convert(muller.col, 0.1), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(asy.intro.a,asy.intro.dc,asy.intro.b,asy.intro.e,asy.intro.f ), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 1))  
vioplot(asy.intro.a,asy.intro.dc,asy.intro.b,asy.intro.e,asy.intro.f ,col=rgb.convert(muller.col, 0.2), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Asymmetrical Introgression", xlab="", lwd=2, names= muller.element, add=T)

plot(log((d$rec/(d$info.sites*d$het)), base=10),asy.intro.dc, col=rgb.convert(col.mode[as.factor(d$mode)], 0.5), pch=16, ylab="Asymmetrical Introgression", xlab=expression(paste(log[10],"(recombination rate / heterozygosity)")), cex=2)

plot(d$rec/(d$info.sites*d$het),asy.intro.dc, col=rgb.convert(col.mode[as.factor(d$mode)], 0.5), pch=16, ylab="Asymmetrical Introgression", xlab=expression(paste(log[10],"(recombination rate / heterozygosity)")), cex=2)


vioplot(dat$ma.fr, dat$mdc.fr,dat$mb.fr,dat$me.fr,dat$mf.fr,col=rgb.convert(muller.col, 0.5), xaxt="n",border = muller.col,horizontal = F, las = 1, cex=1.2, ylab="Mapped reads/Chr length", xlab="", lwd=2, names= muller.element)
stripchart(data.frame(dat$ma.fr, dat$mdc.fr,dat$mb.fr,dat$me.fr,dat$mf.fr), vertical = TRUE,method = "jitter", add = TRUE, pch = 16, cex=1, col=rgb.convert(muller.col, 0.5))  

