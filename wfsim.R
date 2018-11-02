# Code to replicate the simulations described in:
# Karjus et al. 2018, "Challenges in detecting evolutionary forces in language change using diachronic corpora" (ArXiv e-print)
# (run Wright-Fisher simulations with given parameters, bin the 
# output into bins akin to those used in diachronic text corpora,
# and apply the FIT (Feder et al. 2014), following Newberry et al. (2017))
# 
# Code: Andres Karjus

# Required packages
library(dplyr)    
library(RColorBrewer) # just for plot colors
library(compiler)     # recommended: speeds up loops using JIT compilation
enableJIT(3)

##### Functions ####

FITfast = function(af.vec){
  tp.vec = (1:length(af.vec)) # all values present here, no gaps in time
  if(sum(af.vec >= 1 | af.vec <= 0)){ 
    # warning(paste("min, max values are", min(af.vec), max(af.vec), ", must be in (0,1), force-fixing" ))
    af.vec=ifelse(af.vec <= 0, 0.001, ifelse(af.vec >= 1, 0.999, af.vec))
  }
  Yi <- (lead(af.vec) - af.vec)/sqrt(2*af.vec*(1 - af.vec)*(lead(tp.vec) - tp.vec))
  #Remove the NA caused by the lead function
  Yi <- Yi[!is.na(Yi)]
  fitp = t.test(Yi)$p.value
  swp = shapiro.test(Yi)$p.value
  #kolm =  ks.test(Yi,"pnorm")$p.value # kolmogorov-smirnov test
  if(length(Yi)>4){ kolm = lillie.test(Yi)$p.value} else {kolm=NA} # Lilliefors (Kolmogorov-Smirnov) test for the composite hypothesis of normality
  if(length(Yi)>7){ad = ad.test(Yi)$p.value} else { ad=NA } # Anderson-Darling test, min 8
  return( list(fitp=fitp, swp=swp, kolm=kolm, ad=ad) )
}

wfsim = function(s,         # selection coefficient (e.g., 0, 0.1)
                 N=1000,    # population size
                 start=500, # start mutant pop value (e.g. 500 out of 1000 if 0.5 case)
                 len=200    # n generations to simulate (time series length)
){
  # based on some course materials: https://public.wsu.edu/~gomulki/mathgen/materials/wrightfisher.R
  j=c(start, rep(NA, len-1))   # initialize series with starting value
  for(i in 2:len)	{  # simulate W-F model
    p.star = j[i-1]*(1+s)/(j[i-1]*s + N)          # post-selection expected frequency
    j[i]=rbinom(n=1, size=N, prob=min(p.star,1))  # generates random deviates from the binomial distr (yields new number of mutants)
    # n=1 just one  value, size=number of trials (=popsize), prob=probability of success on each trial
    # (floating point errors, need to take care with min() so prob is <=1)
  }
  return(j) # returns the series
}


dofitm = function(N=1000, start=500, reps=100, maxs=5, gens=200,ns=200){
  ss = c(0,exp(seq(log(0.001), log(maxs), length.out = ns-1)))
  maxlen = round(median(table(cut(1:gens, 4)))) # length of bin if minimal n bins
  nx= (1:maxlen)[!duplicated((sapply(1:maxlen, function(x) length(seq(1,gens,x)))))]   # bin lengths, 25
  nbins = sapply(nx, function(x) length(seq(1,gens, x) )) # 25
  
  fitm = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  swm  = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  kolmm  = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  adm  = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  fitratio = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  fitratioswm = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  fitswm = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  for(n in seq_along(nbins)){
    for(s in seq_along(ss)){
      rmvec = rep(NA,reps)
      swvec = rep(NA,reps)
      kolmvec = rep(NA,reps)
      advec   = rep(NA,reps)
      for(rp in 1:reps){ 
        #try({
        j = wfsim(s=ss[s], N=N, start=start, len=gens)
        j = data.frame(j=j, year=1:gens)
        if(n<gens){
          sum.df = j %>% group_by(cut(year, nbins[n])
          ) %>% summarize(value=sum(j), count=n())
          xlin = (sum.df$value / sum.df$count)/N
        } else { # if nothing to bin
          xlin=j$j/N
        }
        tmp = FITfast(xlin)  # FIT function fixes 0,1 events by making them +-0.01 from the boundary
        #})
        rmvec[rp]  = tmp$fitp
        swvec[rp]  = tmp$swp
        kolmvec[rp]=tmp$kolm
        advec[rp]  =tmp$ad
      }
      
      fitm[s, n] = mean(rmvec, na.rm=T) # mean FIT p-value 
      fitratio[s,n] = length(which(rmvec<0.05))/length(which(!is.na(rmvec)))
      swm[s, n] = mean(swvec, na.rm=T)     # mean Shapiro-Wilk p-value
      kolmm[s, n] = mean(kolmvec, na.rm=T) # mean Kolmogorov
      adm[s, n] = mean(advec, na.rm=T)     # mean Andersson
      
      rmvec[swvec<0.1] = NA # set to NA if normality assumption violated
      # record after filtering:
      fitswm[s,n] = mean(rmvec, na.rm=T)   #  mean FIT p-value, but sw-filtered (so if all NA then NA!)
      fitratioswm[s,n] = length(which(rmvec<0.05))/length(which(!is.na(rmvec)))
    }
    cat(".")
  }
  return(list(fitm=fitm, nx=nx, nbins=nbins, fitratio = fitratio, swm=swm, fitswm = fitswm,fitratioswm=fitratioswm, kolmm=kolmm, adm=adm))
}

#### Run ####


# The following calls both return a list with 4 elements:
# [[1]] the matrix of FIT mean p-values, across 200 s-values (rows) and 26 increasing bin lengths/decreasing numbers of bins (columns)
# [[2]] the bin lenghts
# [[3]] the number of bins for each bin length
# [[4]] another matrix: over 100 replicates of each combination, % of FIT p<0.05
# This is the space of s values explored by default: c(0, exp(seq(log(0.01), log(5), length.out = 200-1))

Sys.time() # running these might take a while
# N = pop size; start = mutant starting n; reps = how many replications of each combo; maxs = largest s value to test; gens=how many generations (~years in a corpus); ns = how many s values to test (taken from a log sequence by default, starting with log(0.001), with s=0 prepended to the start)
fitm1 = dofitm(N=1000, start=500, reps=100, maxs=5, gens = 200, ns = 200)
Sys.time()
fitm2 = dofitm(N=1000, start=50, reps=100, maxs=5, gens=200, ns=200)
Sys.time()


#### Some quick plots ####
par(mfrow=c(3,2), mar=c(3,3,1,1))
# Space of selection strength values s explored:
plot( c(0, exp(seq(log(0.001), log(5), length.out = 200-1))), ylab="s", cex=0.5 )

# 50% start, 3 example slices
desat = function(cols, sat=0.5) {
  x = diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(x[1,], x[2,], x[3,])
}
plot(NA, ylim=c(-0.001,0.3), xlim=c(1,25), xlab="", ylab="",yaxs="i", xaxt="n", yaxt="n")
cols = desat(RColorBrewer::brewer.pal(5, "RdYlBu"))
for(i in 1:5){
  rect(0, c(-0.001, 0.01,0.05, 0.1, 0.2)[i], 27,0.4, col= cols[i],border=NA )
}
axis(1, at = 1:25, labels = fitm1$nbins, las=1) # n bins
axis(2, at = c(0.01,0.05,0.1,0.2,0.3), labels=c(0.01,0.05,"0.1","0.2","0.3"), las=1) # pval
pc=c(1,4,17); p=1 
sx = as.numeric(rownames(fitm1[[1]]))
x3  = sapply(c(0.01,0.02,0.1), function(x) which(abs(sx-x)==min(abs(sx-x)))[1] )
for(i in x3){ # 0.0102 0.0203 0.0998
  lines(1:25,fitm1[[1]][i,], type="b", cex=0.8, lwd=1, pch=pc[p])
  p=p+1
}
legend(x=0, y=0.3, legend = c("s=0.1","s=0.02", "s=0.01"), border=NULL, bg=NA,bty = "n", pch=rev(pc)); mtext("number of bins",1, 2.1, cex=m); mtext(expression("Mean FIT "~italic(p)*"-value"),2, 2,cex=0.8)



# Full heatmaps
par(mar=c(1,1,1,1))
# 50% start, power:
image((fitm1[[4]]), breaks=seq(-0.00001,1, length.out = 1000), col=gray.colors(999, start=1,end=0), xaxt="n",yaxt="n")
# 50% start, FIT means:
image((fitm1[[1]]),breaks=c(-1,0.01,0.05,0.1,0.2,1), col=RColorBrewer::brewer.pal(5, "RdYlBu"), xaxt="n",yaxt="n")
# # 5% start, power:
image((fitm2[[4]]), breaks=seq(-0.00001,1, length.out = 1000), col=gray.colors(999, start=1,end=0), xaxt="n",yaxt="n")
# 5% start, FIT means:
image((fitm2[[1]]),breaks=c(-1,0.01,0.05,0.1,0.2,1), col=RColorBrewer::brewer.pal(5, "RdYlBu"), xaxt="n",yaxt="n")
mtext(text="<- low s   high s ->", outer=T,side = 1, line=-1)
mtext(text="<- many short bins   few long bins -> ", outer=T,side = 2, line=-1)


