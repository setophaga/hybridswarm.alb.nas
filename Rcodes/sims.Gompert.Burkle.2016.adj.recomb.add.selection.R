source('~/Desktop/hybridswarm/predictNLS_function.R', chdir = TRUE)
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
  chrbreak<-map$startchr + floor(runif(n=map$nchr,min=round((map$lociperchr)*0.1),max=(map$lociperchr)-1)) # crossover after 1st and before last locus on  chromosome
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
make.nextgen<- function(parent1, parent2, parent1anc, parent2anc, nparind, map,fitness.p1, fitness.p2){
  ## take parental genotype and ancestry arrays [locus,ind,copy], number of desired progeny
  nloci <- dim(parent1)[1]
  off.g<-array(dim=c(nloci,nparind,2)) # genotype #,dimnames=as.list(c("locus","ind","copy"))
  off.a<-array(dim=c(nloci,nparind,2)) # ancestry #,dimnames=as.list(c("locus","ind","copy"))

  for (i in 1:nparind){
    # choose parents at random
 	par1<- sample(1:dim(parent1)[2],1, prob=fitness.p1)
    par2<- ifelse(dim(parent2)[2]>= par1,sample((1:dim(parent2)[2])[-par1],1, prob=fitness.p2[-par1]),sample(1:dim(parent2)[2],1, prob=fitness.p2)) 
    				#(-> no selfing)
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

#add selection
zygote.selection=function(anc, locus1, locus2, s) #specify loci of interest for BDMI
	{fitness=rep(1, length(anc[1,,1])) #start fitness vector with 1 for each individual
	for(i in 1:length(anc[1,, 1]))
		{if(anc[locus1,i,1]==anc[locus1,i, 2] && anc[locus2,i, 1]==anc[locus2,i, 2] && anc[locus1,i, 1]==anc[locus2,i, 1])
			{fitness[i]=1}
		else{fitness[i]=1-s}	  #for heterozygotes and all the homozygotes of alternative ancestry, fitness reduction
		}
	return(fitness)
	}

#########
#########
#########
###make specifications for simulation run 
  nlocikeep<- 100
  nloci<- nlocikeep*3
  ## simulate more loci than kept for output, to ensure equal number
  ## of loci among simulations after filtering loci based on minor
  ## allele frequency (maf)
  #locus ID for BDMI
  locus1=90; locus2=250
  nchr<- 2
  alpha<- 0.8
  fst<- 0.99
  maf<- 0.05
  nindparents<- 200 # 100
  nindnextgen<- 200  #  50
	s= 0.1 #selection against incompatible BDMI combos s= #selection coefficient for incompatibility reduction hybrid fitness =1-s


  map<- make.map(nloci,nchr)
  allelefreq<- sim.afreq(alpha,fst,nloci)
  parentalsams<- sam.parentals(afreq=allelefreq$p,nparind=nindparents)

  P0.fitness=zygote.selection(parentalsams$par0anc, locus1, locus2, s)
  P1.fitness=zygote.selection(parentalsams$par1anc, locus1, locus2, s)
  F1.fitness=zygote.selection(parentalsams$F1anc, locus1, locus2, s)

  F2<- make.nextgen(parent1=parentalsams$F1, parent2=parentalsams$F1,
                    parent1anc=parentalsams$F1anc, parent2anc=parentalsams$F1anc,
                    nparind=nindnextgen, map=map, fitness.p1=F1.fitness, fitness.p2=F1.fitness)
  F2.fitness=zygote.selection(anc=F2$anc, locus1, locus2, s)
           
  BC0<- make.nextgen(parent1=parentalsams$par0, parent2=parentalsams$F1,
                     parent1anc=parentalsams$par0anc, parent2anc=parentalsams$F1anc,
                     nparind=nindnextgen, map=map,fitness.p1=P0.fitness, fitness.p2=F1.fitness)
                     
  BC1<- make.nextgen(parent1=parentalsams$par1, parent2=parentalsams$F1,
                     parent1anc=parentalsams$par1anc, parent2anc=parentalsams$F1anc,
                     nparind=nindnextgen, map=map,fitness.p1=P1.fitness, fitness.p2=F1.fitness)

  #### make F5
  curgen<-F2
  for(i in 1:3){
      FN<- make.nextgen(parent1=curgen$gen, parent2=curgen$gen,
                        parent1anc=curgen$anc, parent2anc=curgen$anc,
                        nparind=nindnextgen, map=map, fitness.p1=zygote.selection(anc=curgen$anc, locus1, locus2, s), fitness.p2=zygote.selection(anc=curgen$anc, locus1, locus2, s))
      curgen<-FN}
  F5<-FN
  for(i in 1:16){
      FN<- make.nextgen(parent1=curgen$gen, parent2=curgen$gen,
                        parent1anc=curgen$anc, parent2anc=curgen$anc,
                        nparind=nindnextgen, map=map, fitness.p1=zygote.selection(anc=curgen$anc, locus1, locus2, s), fitness.p2=zygote.selection(anc=curgen$anc, 	       locus1, locus2, s))
      curgen<-FN}
  F21<-FN
  for(i in 1:6){
      FN<- make.nextgen(parent1=curgen$gen, parent2=curgen$gen,
                        parent1anc=curgen$anc, parent2anc=curgen$anc,
                        nparind=nindnextgen, map=map, fitness.p1=zygote.selection(anc=curgen$anc, locus1, locus2, s), fitness.p2=zygote.selection(anc=curgen$anc, locus1, locus2, s))
      curgen<-FN}
   F27<-FN
    for(i in 1:1){
      FN<- make.nextgen(parent1=curgen$gen, parent2=curgen$gen,
                        parent1anc=curgen$anc, parent2anc=curgen$anc,
                        nparind=nindnextgen, map=map, fitness.p1=zygote.selection(anc=curgen$anc, locus1, locus2, s), fitness.p2=zygote.selection(anc=curgen$anc, locus1, locus2, s))
      curgen<-FN}
   F28<-FN

                    
  ## add population reference  together
  ## (change subscripts to match cat*s above, if number of inds in each category changed)
  ## 40 par0, 40 par1, 20 F1, 10 F2, 10 F5, 10 BC0 (BC1 to parent 0) and 10 BC1

  ## genotypes (a1 is allele 1, and a2 is allele 2)
  all.genotypes.a1<- cbind(parentalsams$par0[,1:20,1],parentalsams$par1[,1:20,1],
                           parentalsams$F1[,1:10,1],
                           F2$gen[,1:10,1], F5$gen[,1:10,1], F21$gen[,1:50,1],
                           BC0$gen[,1:10,1], BC1$gen[,1:10,1],F27$gen[,1:50,1], F28$gen[,1:50,1])

  all.genotypes.a2<- cbind(parentalsams$par0[,1:20,2], parentalsams$par1[,1:20,2],
                           parentalsams$F1[,1:10,2],
                           F2$gen[,1:10,2], F5$gen[,1:10,2],  F21$gen[,1:50,2],
                           BC0$gen[,1:10,2], BC1$gen[,1:10,2], F27$gen[,1:50,2], F28$gen[,1:50,2])

  all.genotypes<- array(dim=c(dim(all.genotypes.a1),2))
  all.genotypes[,,1]<- all.genotypes.a1
  all.genotypes[,,2]<- all.genotypes.a2

  ## ancestries
  all.ancestries.a1<- cbind(parentalsams$par0anc[,1:20,1],parentalsams$par1anc[,1:20,1],
                            parentalsams$F1anc[,1:10,1],
                            F2$anc[,1:10,1], F5$anc[,1:10,1],  F21$anc[,1:50,1],
                            BC0$anc[,1:10,1], BC1$anc[,1:10,1], F27$anc[,1:50,1], F28$anc[,1:50,1])
  all.ancestries.a2<- cbind(parentalsams$par0anc[,1:20,2], parentalsams$par1anc[,1:20,2],
                            parentalsams$F1anc[,1:10,2],
                            F2$anc[,1:10,2], F5$anc[,1:10,2], F21$anc[,1:50,2],
                            BC0$anc[,1:10,2], BC1$anc[,1:10,2],F27$anc[,1:50,2], F28$anc[,1:50,2])
  all.ancestries<- array(dim=c(dim(all.ancestries.a1),2))
  all.ancestries[,,1]<- all.ancestries.a1
  all.ancestries[,,2]<- all.ancestries.a2
  ##-------------------
  ## specify number of moms, dads, and reference samples from each category
  ## "BC0" is F1 x par0;  "BC1" is F1 x par1
  categories<- c("par0","par1","F1","F2","F5","F21","BC0","BC1","F27", "F28")
  catrefs<- c(20,20,10,10,10,50,10,10, 50, 50)
  ids<-rep(categories, catrefs)
  nrefs<- sum(catrefs)

    anc.all= colSums(all.ancestries[,,1] + all.ancestries[,,2]) /(nloci*2)
    het.all=apply(all.ancestries[,,1]+all.ancestries[,,2], 2, function(x) sum(x==1)) / nloci
   
tmp=list(ids=ids, cols=cols, symbols=symbols, anc.indv=anc.all, het=het.all, geno=all.genotypes[,,], anc=all.ancestries[,,])
tmp$ids

#########
#########
#########.                  ANCESTRY ANALYSIS AS EMPIRICAL DATA
#########.   
#########.  
##take genotype matrix
#sequence relative to 
 p1=1:20; p2=21:40;f1=41:50;f2=51:60;f5=61:70;f21=71:120;bc0=121:130;bc1=131:140;f27=141:190; f28=191:240
 
h.all=41:length(tmp$anc[1,,1])
anc=(tmp$anc[,,1]  + tmp$anc[,,2])/2
p1.anc=anc[,p1]
p2.anc=anc[,p2]
h.an=anc[,h.all]
 hist(p1.anc)
colnames(h.an)=tmp$ids[-c(p1, p2)]
hybrid.cat=tmp$ids[-c(p1, p2)]
col <- colorRampPalette(c("turquoise","grey","deepskyblue"))(100)
heatmap(as.matrix(t(h.an)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(h.an[c(locus1,locus2), which(hybrid.cat=="F27")])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(h.an[c(1,151), which(hybrid.cat=="F27")])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)

heatmap(as.matrix(t(anc[, f1])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, f2])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, f5])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, f21])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, bc0])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, bc1])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, f27])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(anc[, f28])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)


heatmap(as.matrix(t(h.an[, which(hybrid.cat=="F21")])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)
heatmap(as.matrix(t(h.an[, which(hybrid.cat=="F21")])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)

heatmap(as.matrix(t(h.an[, which(hybrid.cat=="F27")])), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, cexCol=0.5)


#check if this ancestry representative loci is consistent with known ancestry for each individual
h.hi=tmp$anc.ind[41:length(tmp$anc[1,,1])]; plot(h.hi, colMeans(h.an))
plot(tmp$anc.ind,(colMeans(tmp$anc[,,1])+colMeans(tmp$anc[,,2]))/2)

 #name relative to the full anc
p1=which(tmp$ids=="par0");p2=which(tmp$ids=="par1");f1=which(tmp$ids=="F1");f2=which(tmp$ids=="F2");f5=which(tmp$ids=="F5");f21=which(tmp$ids=="F21");bc0=which(tmp$ids=="BC0");bc1=which(tmp$ids=="BC1");f27=which(tmp$ids=="F27"); f28=which(tmp$ids=="F28")

f1=h.an[,which(hybrid.cat=="F1")];f2=h.an[,which(hybrid.cat=="F2")];f5=h.an[,which(hybrid.cat=="F5")]; f21=h.an[,which(hybrid.cat=="F21")];bc0=h.an[,which(hybrid.cat=="BC0")]; bc1=h.an[,which(hybrid.cat=="BC1")]; f27=h.an[,which(hybrid.cat=="F27")];f28=h.an[,which(hybrid.cat=="F28")]
heatmap(as.matrix(t(f28)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, main="Gen28")
heatmap(as.matrix(t(f27)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, main="Gen27")
heatmap(as.matrix(t(f21)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, main="Gen21")
heatmap(as.matrix(t(f2)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, main="Gen2")
heatmap(as.matrix(t(f5)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3, main="Gen5")

# heatmap(as.matrix(t(p2.geno)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)
# heatmap(as.matrix(t(p1.geno)), Colv = NA, Rowv = NA, col = col, scale = "none",cexRow=0.3)

#fit BGC
phi=function(a, b, h)
{phi=h+2*(h-h^2)*(a+b*(2*h-1))
return(phi)	}

# h=colMeans(h.an)
bgc.gen=function(h.hi, h.an)
{h=h.hi
coefi={}
par(mfrow = c(5, 5), mar=c(1,1,1,1))
for(i in 1:nrow(h.an))
{p=h.an[i,]
plot(h, p)
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm="port", lower=c(a=-1, b=-1), upper=c(a=1,b=1),control=nls.control(maxiter = 100, warnOnly=TRUE))
a=coef(m)[1]; b=coef(m)[2]; se.a=coef(summary(m))[1, "Std. Error"]; se.b=coef(summary(m))[2, "Std. Error"]
coefi=rbind(coefi, c(a,se.a, b, se.b))
x<-data.frame(h=seq(0, 1, length.out=1000));
y=predict(m, x)
abline(0,1, lty=2, lwd=1.5)
lines(x[,1], y, col="blue", lwd=2)
text(0.2,1, paste("alpha= ", round(a, 3)), cex=0.7)
text(0.2,0.95, paste("beta= ", round(b, 3)), cex=0.7)
}
return(coefi)
}

#all.h=bgc.gen(h.hi, h.an)
#seperate hybrid generations
f1.hi=h.hi[which(hybrid.cat=="F1")];f2.hi=h.hi[which(hybrid.cat=="F2")];f5.hi=h.hi[which(hybrid.cat=="F5")]; 
f21.hi=h.hi[which(hybrid.cat=="F21")];f27.hi=h.hi[which(hybrid.cat=="F27")];f28.hi=h.hi[which(hybrid.cat=="F28")]

#pdf("gen21.bgc.simulation.10percent.inv.pdf")
f21.bgc=bgc.gen(f21.hi, f21)
#dev.off()
#pdf("gen27.bgc.simulation.10percent.inv.pdf")
f27.bgc=bgc.gen(f27.hi, f27)
f28.bgc=bgc.gen(f28.hi, f28)
f5.bgc=bgc.gen(f5.hi, f5)
dev.off()
hist(f27.bgc[,3])
hist(f28.bgc[,3])
hist(f21.bgc[,3])
hist(f21.bgc[,3])
hist(f5.bgc[,3])


###make a BGC plot of the first locus
p=h.an[1,]; h=h.hi
plot(h, p)
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm="port", lower=c(a=-1, b=-1), upper=c(a=1,b=1),control=nls.control(maxiter = 100, warnOnly=TRUE), )
a=coef(m)[1]; b=coef(m)[2]; se.a=coef(summary(m))[1, "Std. Error"]; se.b=coef(summary(m))[2, "Std. Error"]
x<-data.frame(h=seq(0, 1, length.out=100));
ci<-predictNLS(m, newdata=x)
lci<-ci[,6]
uci<-ci[,7]
y<-ci[,1]
i.for <- order(x$h)
i.back<- order(x$h , decreasing = TRUE )
x.polygon <- c( x$h[i.for] , x$h[i.back] )
y.polygon <- c( uci[i.for] , lci[i.back])
polygon( x.polygon , y.polygon , col = rgb(0/255, 191/255, 255/255, 0.4) , border = NA )
abline(0,1, lty=2, lwd=1.5)
lines(x[,1], y, col="blue", lwd=2)

###do k.means clusters
c.an=kmeans(h.an, 2)
k=2
while(c.an$betweenss/c.an$totss<0.95)
{k=k+1;c.an=kmeans(h.an, k)}

#take the order of clusters as the sequence along the genome
dc.cls=factor(c.an$cluster, levels=unique(c.an$cluster))
####################################             check above here !!!!         
#now use this order for the cluster factor dc.cls
#now with c.an$cluster and h.an, calculate the mean ancestry for each cluster for each individual
clusters=unique(c.an$cluster)
cd.cl={}
for(j in 1:ncol(h.an)) #loop among individuals
	{cl={} #fresh vector containing mean ancestry for each cluster 
	for(c in clusters) #loop among clusters
		{cl=c(cl,mean(h.an[which(c.an$cluster==c),j]))}
	cd.cl=rbind(cd.cl, cl) #attach cluster ancestry row for each indivudual
	}
colnames(cd.cl)=paste("c", clusters)
rownames(cd.cl)=rownames(c.an)
heatmap(as.matrix(cd.cl), Colv = NA, Rowv = NA, scale="none", col=col, cexRow=0.5,cexCol=0.5)
h.c=cd.cl
#take out rows for individuals from each generate 
f21.c=h.c[which(hybrid.cat=="F21"),];f27.c=h.c[which(hybrid.cat=="F27"),];f28.c=h.c[which(hybrid.cat=="F28"),]

bgc.gen.cluster=function(h.an)
{h=rowMeans(h.an) #hybrid index of each individual
coefi={}
par(mfrow = c(5, 5), mar=c(1,1,1,1))
for(j in 1:ncol(h.an)) #loop through each cluster
{p=h.an[,j] #take the column for each genetic cluster
plot(h, p)
m=nls(p~h+2*(h-h^2)*(a+b*(2*h-1)), start=list(a=0, b=0), algorithm="port", lower=c(a=-1, b=-1), upper=c(a=1,b=1),control=nls.control(maxiter = 100, warnOnly=TRUE))
a=coef(m)[1]; b=coef(m)[2]; se.a=coef(summary(m))[1, "Std. Error"]; se.b=coef(summary(m))[2, "Std. Error"]
coefi=rbind(coefi, c(a,se.a, b, se.b))
x<-data.frame(h=seq(0, 1, length.out=1000));
y=predict(m, x)
abline(0,1, lty=2, lwd=1.5)
lines(x[,1], y, col="blue", lwd=2)
text(0.2,1, paste("alpha= ", round(a, 3)), cex=0.7)
text(0.2,0.95, paste("beta= ", round(b, 3)), cex=0.7)
}
return(coefi)
}
f27.bgcc=bgc.gen.cluster(f27.c);
f28.bgcc=bgc.gen.cluster(f28.c);
dev.off()
hist(f27.bgcc); 
hist(f28.bgcc)
n =length(clusters)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
par(mfrow = c(2, 1),mar = c(1,1,1,1) + 0.2,oma = c(5,4,1,1) + 0.2)
plot(seq(1:length(f27.bgcc[,3])),f27.bgcc[,3], pch=21, bg=color[clusters], ylab=expression(beta)); abline(h=0)
plot(f27.bgc[,3], pch=21, bg=color[c.an$cluster],ylab=expression(beta)); abline(h=0)
abline(v=locus1, lwd=2)
abline(v=locus2, lwd=2)
lines(seq(1:length(f27.bgc[,3])),f27.bgc[,3])
dev.off()

par(mfrow = c(2, 1),mar = c(1,1,1,1) + 0.2,oma = c(5,4,1,1) + 0.2)
plot(seq(1:length(f28.bgcc[,3])),f28.bgcc[,3], pch=21, bg=color[clusters], ylab=expression(beta)); abline(h=0)
plot(f28.bgc[,3], pch=21, bg=color[c.an$cluster],ylab=expression(beta)); abline(h=0)
abline(v=locus1, lwd=2)
abline(v=locus2, lwd=2)
lines(seq(1:length(f28.bgc[,3])),f28.bgc[,3])
dev.off()

par(mfrow = c(2, 1))
plot(f27.bgcc[,1], pch=21, bg=col_vector[clusters+5], ylab=expression(alpha)); abline(h=0)
plot(f27.bgc[,1], pch=21, bg=col_vector[c.an$cluster+5],ylab=expression(alpha)); abline(h=0)
dev.off()


f21.bgcc=bgc.gen.cluster(f21.c);
par(mfrow = c(2, 1))
plot(f21.bgcc[,3], pch=21, bg=color[clusters], ylab=expression(beta)); abline(h=0)
plot(f21.bgc[,3], pch=21, bg=color[c.an$cluster],ylab=expression(beta)); abline(h=0)
dev.off()

# f27.bgcc, f27.c
#test if LD change analysis would be affected by inversion
cutoff=0.5
f27.c;f21.c
outliers=which(f27.bgcc[,3]>cutoff)#take clusters with beta greater than cutoff
others=c(which(f27.bgcc[,3]<cutoff, f27.bgcc[,3]==cutoff))#take clusters with beta greater than cutoff
#run bootstrap 
N=10000;  r2.diff={}; r2.diff.other={}
for(per in 1:N)
	{g.c=sample(outliers, 2)
	g.o=sample(others, 2)
	r2.27=cor(f27.c[,g.c[1]], f27.c[,g.c[2]])^2
	r2.21=cor(f21.c[,g.c[1]], f21.c[,g.c[2]])^2
	r2.diff=c(r2.diff, (r2.27-r2.21))
	#calculate change for genome-wide clusters
	r2.27.other=cor(f27.c[,g.o[1]], f27.c[,g.o[2]])^2
	r2.21.other=cor(f21.c[,g.o[1]], f21.c[,g.o[2]])^2
	r2.diff.other=c(r2.diff.other, (r2.27.other-r2.21.other))
	}
t.test(r2.diff.other, r2.diff)
hist(r2.diff, xlim=c(-1, 0.8),col=rgb(1,0.843137,0,0.3), breaks=100, border=rgb(1,0.843137,0,0.6), main="simulation", xlab=bquote(Delta~r^2),ylim=c(0, 1500))
hist(r2.diff.other, xlim=c(-1, 0.8),add=T,col=rgb(0,0,0.8,0.2), breaks=100, border=rgb(0,0,0.8,0.2))
hist(r2.diff, xlim=c(-1, 0.8),col=rgb(1,0.843137,0,0.3), breaks=100, border=rgb(1,0.843137,0,0.3), add=T)
legend(-1, 1500,  c("Barrier clusters", "Control clusters"), fill=c(rgb(1,0.843137,0,0.3), rgb(0,0,0.8,0.2)))

