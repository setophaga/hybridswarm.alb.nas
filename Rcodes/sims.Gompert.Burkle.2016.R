## adapted from Doro's scripts for entropy and Populus, which in turn
## were adapted from simulations that Alex wrote.

## do simulations with specified parameters
## keep track of locus ancestries and genotypes
## for different numbers of moms, progeny, dads, reference samples code below must be changed

## contrast parentals F1, BC1, F2, and F3.  In particular, a) show
## range of variation within each class--note that these are
## expectations for individuals, not population expectations, except
## under certain conditions [namely equal abundance of
## parentals]--make connection to NewHybrids papers, b) show the
## effect of shared polymorphism within each species on the range of
## variation in genetic composition (so examine genetic composition in
## addition to ancestry variation).  This applies even to F1s



## do we want to consider what do population mixtures look like?

##  for the piece involving allele freq rather than ancestry
##  ... histogram of allele frequency in parentals and F1

### also some comments on populations of hybrids and the inference of dynamics

### ability of taxonomists to correctly id hybrids in addition to species: Anderson, Heiser, Grant, Nabokov

### gene flow v. hybridization continuum

##----------------------------------

make.map<- function(nloci,nchr){
  start.chr<-seq(1,nloci,nloci/nchr)
  end.chr<-seq(nloci/nchr,nloci,nloci/nchr)
  locusids<- paste(rep(1:nchr,each=nloci/nchr),rep(1:(nloci/nchr),times=nchr),sep=":")
  list(nchr=nchr,lociperchr= nloci/nchr, startchr=start.chr, endchr=end.chr, locusids=locusids)
}

## simulate cluster allele frequencies
sim.afreq<- function(alpha, fst, nloci){
  sim.pi <- rbeta(nloci, alpha, alpha)
  sim.pi[sim.pi==1] <- 0.9999999
  sim.pi[sim.pi==0] <- 0.0000001

  sim.p <- matrix(nrow=nloci, ncol=2) # 2 clusters <<<< ---- note we are only considering k=2
  for(k in 1:2){
    sim.p[,k] <- rbeta(nloci,
                       shape1 = sim.pi     * (-1 + 1/fst),
                       shape2 = (1-sim.pi) * (-1 + 1/fst))
  }

  sim.p[sim.p>0.9999999] <- 0.9999999
  sim.p[sim.p<0.0000001] <- 0.0000001

  return(list(pi=sim.pi, p=sim.p))
}


## sample genotypes for par0, par1, F1
sam.parentals<- function(afreq,nparind){
  nloci <- dim(afreq)[1]
  # genotypes
  par0.g<-array(dim=c(nloci,nparind,2))
  par1.g<-array(dim=c(nloci,nparind,2))#,dimnames=c("locus","ind","copy")
  F1.g<-array(dim=c(nloci,nparind,2))#,dimnames=as.list(c("locus","ind","copy"))

  par0.g[,,1]<- matrix(rbinom(nloci*nparind, 1, prob=afreq[,1]),nrow=nloci, ncol=nparind)
  par0.g[,,2]<- matrix(rbinom(nloci*nparind, 1, prob=afreq[,1]),nrow=nloci, ncol=nparind)

  par1.g[,,1]<- matrix(rbinom(nloci*nparind, 1, prob=afreq[,2]),nrow=nloci, ncol=nparind)
  par1.g[,,2]<- matrix(rbinom(nloci*nparind, 1, prob=afreq[,2]),nrow=nloci, ncol=nparind)

  F1.g[,,1]<- matrix(rbinom(nloci*nparind, 1, prob=afreq[,1]),nrow=nloci, ncol=nparind)
  F1.g[,,2]<- matrix(rbinom(nloci*nparind, 1, prob=afreq[,2]),nrow=nloci, ncol=nparind)

  # ancestries
  par0.a<-array(data=0,dim=c(nloci,nparind,2))#,dimnames=as.list(c("locus","ind","copy"))
  par1.a<-array(data=1,dim=c(nloci,nparind,2))#,dimnames=as.list(c("locus","ind","copy"))
  F1.a<-array(dim=c(nloci,nparind,2))#,dimnames=as.list(c("locus","ind","copy"))
  F1.a[,,1]<- par0.a[,,1]
  F1.a[,,2]<- par1.a[,,1]
  return(list(par0=par0.g, par1=par1.g, F1=F1.g, par0anc=par0.a, par1anc=par1.a, F1anc=F1.a))
}


## sample gametes (make crossover)
make.gamete<- function(genotype, ancestry, map){ # only take 1 individual[loci,copy]
  nloci<- map$nchr * map$lociperchr
  chrindex<-numeric(nloci)
  chrbreak<-map$startchr + floor(runif(n=map$nchr,min=0,max=(map$lociperchr)-1)) # crossover after 1st and before last locus on  chromosome
  gen.gamete<-numeric(nloci) # genotype
  anc.gamete<-numeric(nloci) # ancestry

  # identify crossover breakpoints
  for (chr in 1:map$nchr){
    copy<-sample(c(1,2),1) # start of chromosome block comes from first or second copy at random
    if (copy == 1) {
      chrindex[map$startchr[chr]:chrbreak[chr]]<- 1
      chrindex[(chrbreak[chr]+1):map$endchr[chr]]<- 2
    }
    else if (copy == 2){
      chrindex[map$startchr[chr]:chrbreak[chr]]<- 2
      chrindex[(chrbreak[chr]+1):map$endchr[chr]]<- 1
    }
  }
  gen.gamete[which(chrindex==1)]<- genotype[which(chrindex==1),1] # take allele from first copy
  gen.gamete[which(chrindex==2)]<- genotype[which(chrindex==2),2] # take allele from second copy

  anc.gamete[which(chrindex==1)]<- ancestry[which(chrindex==1),1] # take ancestry from first copy
  anc.gamete[which(chrindex==2)]<- ancestry[which(chrindex==2),2] # take ancestry from second copy

  return(list(gamete=gen.gamete,ancestry=anc.gamete))
}


## make next generation (repeat sample gametes and make offspring for 2nd generation, i.e. F2, BC0, BC1)
make.nextgen<- function(parent1, parent2, parent1anc, parent2anc, nparind, map){
  ## take parental genotype and ancestry arrays [locus,ind,copy], number of desired progeny
  nloci <- dim(parent1)[1]
  off.g<-array(dim=c(nloci,nparind,2)) # genotype #,dimnames=as.list(c("locus","ind","copy"))
  off.a<-array(dim=c(nloci,nparind,2)) # ancestry #,dimnames=as.list(c("locus","ind","copy"))

  for (i in 1:nparind){
    # choose parents at random
    par1<- sample(1:dim(parent1)[2],1)
    par2<- ifelse(dim(parent2)[2]>= par1,sample((1:dim(parent2)[2])[-par1],1),sample(1:dim(parent2)[2],1)) #(-> no selfing)

    gamete1<- make.gamete(genotype=parent1[,par1,], ancestry=parent1anc[,par1,],map) # first copy from parent 1
    gamete2<- make.gamete(genotype=parent2[,par2,], ancestry=parent2anc[,par2,],map) # second copy from parent 2

    off.g[,i,1] <- gamete1$gamete
    off.g[,i,2] <- gamete2$gamete
    off.a[,i,1] <- gamete1$ancestry
    off.a[,i,2] <- gamete2$ancestry
  }
  return(list(gen=off.g,anc=off.a))
}

require(RColorBrewer)

runsim<-function(){
  nlocikeep<- 100
  nloci<- nlocikeep*10
  ## simulate more loci than kept for output, to ensure equal number
  ## of loci among simulations after filtering loci based on minor
  ## allele frequency (maf)

  nchr<- 2
  alpha<- 0.8
  fst<- 0.5
  maf<- 0.05
  nindparents<- 200 # 100
  nindnextgen<- 200  #  50

  ##-------------------
  ## specify number of moms, dads, and reference samples from each category
  ## "BC0" is F1 x par0;  "BC1" is F1 x par1
  categories<- c("par0","par1","F1","F2","F5","F20","BC0","BC1")
  catrefs<- c(20,20,20,10,10,10,10,10)
  ids<-rep(categories, catrefs)
  tmpcols<-brewer.pal(7, "RdYlBu")
  cols<-rep(c(tmpcols[1], tmpcols[7], tmpcols[3], tmpcols[3], tmpcols[3], tmpcols[3], tmpcols[2], tmpcols[6]), catrefs)
  symbols<-rep(c(21,21,21,22,24,25,21,21), catrefs)
  nrefs<- sum(catrefs)

  map<- make.map(nloci,nchr)
  allelefreq<- sim.afreq(alpha,fst,nloci)
  parentalsams<- sam.parentals(afreq=allelefreq$p,nparind=nindparents)
  F2<- make.nextgen(parent1=parentalsams$F1, parent2=parentalsams$F1,
                    parent1anc=parentalsams$F1anc, parent2anc=parentalsams$F1anc,
                    nparind=nindnextgen, map=map)
  BC0<- make.nextgen(parent1=parentalsams$par0, parent2=parentalsams$F1,
                     parent1anc=parentalsams$par0anc, parent2anc=parentalsams$F1anc,
                     nparind=nindnextgen, map=map)
  BC1<- make.nextgen(parent1=parentalsams$par1, parent2=parentalsams$F1,
                     parent1anc=parentalsams$par1anc, parent2anc=parentalsams$F1anc,
                     nparind=nindnextgen, map=map)

  #### make F5
  curgen<-F2
  for(i in 1:3){
      FN<- make.nextgen(parent1=curgen$gen, parent2=curgen$gen,
                        parent1anc=curgen$anc, parent2anc=curgen$anc,
                        nparind=nindnextgen, map=map)
      curgen<-FN}
  F5<-FN
  for(i in 1:15){
      FN<- make.nextgen(parent1=curgen$gen, parent2=curgen$gen,
                        parent1anc=curgen$anc, parent2anc=curgen$anc,
                        nparind=nindnextgen, map=map)
      curgen<-FN}
  F20<-FN

  ## add population reference  together
  ## (change subscripts to match cat*s above, if number of inds in each category changed)
  ## 40 par0, 40 par1, 20 F1, 10 F2, 10 F5, 10 BC0 (BC1 to parent 0) and 10 BC1

  ## genotypes (a1 is allele 1, and a2 is allele 2)
  all.genotypes.a1<- cbind(parentalsams$par0[,1:20,1],parentalsams$par1[,1:20,1],
                           parentalsams$F1[,1:20,1],
                           F2$gen[,1:10,1], F5$gen[,1:10,1], F20$gen[,1:10,1],
                           BC0$gen[,1:10,1], BC1$gen[,1:10,1])

  all.genotypes.a2<- cbind(parentalsams$par0[,1:20,2], parentalsams$par1[,1:20,2],
                           parentalsams$F1[,1:20,2],
                           F2$gen[,1:10,2], F5$gen[,1:10,2],  F20$gen[,1:10,2],
                           BC0$gen[,1:10,2], BC1$gen[,1:10,2])

  all.genotypes<- array(dim=c(dim(all.genotypes.a1),2))
  all.genotypes[,,1]<- all.genotypes.a1
  all.genotypes[,,2]<- all.genotypes.a2

  ## ancestries
  all.ancestries.a1<- cbind(parentalsams$par0anc[,1:20,1],parentalsams$par1anc[,1:20,1],
                            parentalsams$F1anc[,1:20,1],
                            F2$anc[,1:10,1], F5$anc[,1:10,1],  F20$anc[,1:10,1],
                            BC0$anc[,1:10,1], BC1$anc[,1:10,1])
  all.ancestries.a2<- cbind(parentalsams$par0anc[,1:20,2], parentalsams$par1anc[,1:20,2],
                            parentalsams$F1anc[,1:20,2],
                            F2$anc[,1:10,2], F5$anc[,1:10,2], F20$anc[,1:10,2],
                            BC0$anc[,1:10,2], BC1$anc[,1:10,2])
  all.ancestries<- array(dim=c(dim(all.ancestries.a1),2))
  all.ancestries[,,1]<- all.ancestries.a1
  all.ancestries[,,2]<- all.ancestries.a2

  ## filter loci based on minor allele frequencies across all inds
  af<-rowSums( all.genotypes[,,1] +  all.genotypes[,,2] ) /(2*nrefs) # allele frequencies in sample

  cat("   maf loci:", sum(af >= maf & af <= 1-maf), "  rare loci:",
      sum(af < maf | af > 1-maf), fill=T)

  keep.maf<- sort(sample( which(af >= maf & af <= 1-maf), nlocikeep))
  keep.rare<- sort(sample( which(af < maf | af > 1-maf), nlocikeep))

  ## true intra-source ancestry (q; only for loci in keep.maf), mean ancestry
  anc.pointest.maf<- colSums(all.ancestries[keep.maf,,1] + all.ancestries[keep.maf,,2]) /
    (length(keep.maf)*2)
  anc.pointest.rare<- colSums(all.ancestries[keep.rare,,1] + all.ancestries[keep.rare,,2]) /
    (length(keep.rare)*2)

  ## true inter-source ancestry (Q12; only for loci in keep.maf)
  het.maf <- apply(all.ancestries[keep.maf,,1]+all.ancestries[keep.maf,,2], 2, function(x) sum(x==1)) / length(keep.maf)
  het.rare <- apply(all.ancestries[keep.rare,,1]+all.ancestries[keep.rare,,2], 2, function(x) sum(x==1)) / length(keep.rare)


  return( list(ids=ids, cols=cols, symbols=symbols, anc.maf=anc.pointest.maf, anc.rare=anc.pointest.rare,
               het.maf=het.maf, het.rare=het.rare, geno.maf=all.genotypes[keep.maf,,],
               geno.rare=all.genotypes[keep.rare,,]))
}


tmp<-runsim()

do.pca<-function(gmat, write.gcov=FALSE, inds=""){
    gmn<-apply(gmat,1,mean, na.rm=T)
    gmnmat<-matrix(gmn,nrow=nrow(gmat),ncol=ncol(gmat))
    gprime<-gmat-gmnmat ## remove mean (i.e., center)

    gcovarmat<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
    for(i in 1:ncol(gmat)){
        for(j in i:ncol(gmat)){
            if (i==j){
                gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
            }
            else{
                gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
                gcovarmat[j,i]<-gcovarmat[i,j]
            }
        }
    }
    if(write.gcov==TRUE){
        inds<-ifelse(inds == "", paste("i", 1:ncol(gmat)), inds)
        write.table(round(gcovarmat,5),file="gcovarmat.txt",
                    quote=F,row.names=F,col.names=inds)
    }
    prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
}

pca.result<-do.pca(tmp$geno.maf[,,1]  + tmp$geno.maf[,,2])
pcSummary<-summary(pca.result)


par(mar=c(4.5,4.5,0.5,1.5), pty="s", mfrow=c(1,2))
plot(tmp$anc.maf, tmp$het.maf, xlab="admixture proportion (q)",
     ylab=expression(inter-population ~ ancestry ~ (Q[12])), type="n", axes=F)
axis(1, at=c(0,0.5,1))
axis(2, at=c(0,0.5,1))
lines(c(0,0.5,1), c(0,1,0))
points(tmp$anc.maf, tmp$het.maf, col="darkgray", bg=tmp$cols, pch=tmp$symbols)
text(0.02, 0, "Taxon 1", pos=4)
text(0.98, 0, "Taxon 2", pos=2)
text(0.5, 0.97, expression(F[1]), pos=1)
text(0.23, 0.5, expression(BC[0]), pos=3, srt=67.5)
text(0.77, 0.5, expression(BC[1]), pos=3, srt=-67.5)
text(0.5, 0.35, expression(F[2] ~ "," ~ F[5] ~ "&" ~ F[20]), pos=1)

plot(pca.result$x[,'PC1'], pca.result$x[,'PC2'],
     col="darkgray", bg=tmp$cols, pch=tmp$symbols,
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""))

dev.print(pdf, "qQ_PCA.pdf", width=8, height=4.5)

summary(tmp)
##take genotype matrix
geno=tmp$geno.maf[,,1]  + tmp$geno.maf[,,2]
p1.geno=geno[,1:20]
p2.geno=geno[,21:40]
h.geno=geno[,41:110]
hist(p1.geno)
#find species-diagnostic loci
an.loci=which(rowSums(abs(p1.geno-p2.geno))>(1.5*ncol(p1.geno))) #max =2*20=40. but not all the sites have to be fixed, so cut off 30 would be good enough
p1.g=p1.geno[an.loci,]
p2.g=p2.geno[an.loci,]
h.g=h.geno[an.loci,]
p1.an=p1.g; p2.an=p2.g; h.an=h.g
for(i in 1:length(an.loci))
{if(sum(p1.g[i,])>(1.2*ncol(p1.geno))) #2 is p1 genotype should be changed to 0
	{p1.an[i,which(p1.g[i,]==2)]=0
	p1.an[i,which(p1.g[i,]==1)]=0.5
	p1.an[i,which(p1.g[i,]==0)]=1
	p2.an[i,which(p2.g[i,]==0)]=1
	p2.an[i,which(p2.g[i,]==1)]=0.5
	p2.an[i,which(p2.g[i,]==2)]=0
	h.an[i,which(h.g[i,]==2)]=0
	h.an[i,which(h.g[i,]==1)]=0.5
	h.an[i,which(h.g[i,]==0)]=1
	}
else{p1.an[i,which(p1.g[i,]==2)]=1 #2 is p2 ancestry should be changed to 1
	p1.an[i,which(p1.g[i,]==1)]=0.5
	p2.an[i,which(p2.g[i,]==1)]=0.5
	p2.an[i,which(p2.g[i,]==2)]=1
	h.an[i,which(h.g[i,]==2)]=1
	h.an[i,which(h.g[i,]==1)]=0.5
	h.an[i,which(h.g[i,]==0)]=0} 
}
col <- colorRampPalette(c("turquoise","grey","deepskyblue"))(100)
heatmap(as.matrix(t(h.an)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)
heatmap(as.matrix(t(p1.an)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)
heatmap(as.matrix(t(p2.an)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)
heatmap(as.matrix(t(p2.geno)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)
heatmap(as.matrix(t(p1.geno)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)

h.hi=tmp$anc.maf[41:110]
plot(h.hi, colMeans(h.an))

#fit BGC
phi=function(a, b, h)
{phi=h+2*(h-h^2)*(a+b*(2*h-1))
return(phi)	}

# h=colMeans(h.an)
h=h.hi
coefi={}
par(mfrow = c(3, 7))
for(i in 1:nrow(h.an))
{p=h.an[i,]
plot(h, p)
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm="port", lower=c(a=-1, b=-1), upper=c(a=1,b=1),control=nls.control(maxiter = 100, warnOnly=TRUE))
a=coef(m)[1]; b=coef(m)[2]; se.a=coef(summary(m))[1, "Std. Error"]; se.b=coef(summary(m))[2, "Std. Error"]
coefi=rbind(coefi, c(a,se.a, b, se.b))
x<-data.frame(h=seq(0, 1, length.out=1000));
y<-predict(m,x)
abline(0,1, lty=2, lwd=1.5)
lines(x[,1], y, col="blue", lwd=2)
text(0.2,1, paste("alpha= ", round(a, 3)), cex=0.7)
text(0.2,0.95, paste("beta= ", round(b, 3)), cex=0.7)
}

