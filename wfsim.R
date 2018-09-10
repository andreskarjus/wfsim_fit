# Code to run Wright-Fisher simulations with given parameters, bin the 
# output into bins akin to those used in diachronic text corpora,
# and apply the FIT (Feder et al. 2014), following Newberry et al. (2017).
# 
# Code: Andres Karjus

# Required packages
library(dplyr)
library(RColorBrewer)
library(compiler) # recommended: speeds up loops using JIT compilation
enableJIT(3)

##### Functions ####
FITfast = function(af.vec){
  # based on p.c. with author
  tp.vec = (1:length(af.vec))
  if(sum(af.vec >= 1 | af.vec <= 0)){ 
    # warning(paste("min, max values are", min(af.vec), max(af.vec), ", must be in (0,1), force-fixing" ))
    af.vec=ifelse(af.vec <= 0, 0.001, ifelse(af.vec >= 1, 0.999, af.vec))
  }
  Yi <- (lead(af.vec) - af.vec)/sqrt(2*af.vec*(1 - af.vec)*(lead(tp.vec) - tp.vec))
  #Remove the NA caused by the lead function
  Yi <- Yi[!is.na(Yi)]
  return(t.test(Yi)$p.value)
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


dofitm = function(N=1000, start=500, reps=100, maxs=5, gens=200){
  ss = c(0,exp(seq(log(0.01), log(maxs), length.out = gens-1)))
  maxlen = round(median(table(cut(1:gens, 3)))) # length of bin if minimal n bins
  nx= (1:maxlen)[!duplicated((sapply(1:maxlen, function(x) length(seq(1,gens,x)))))]   # bin lengths
  nbins = sapply(nx, function(x) length(seq(1,gens, x) ))
  
  fitm = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  fitratio = matrix(NA, length(ss), length(nbins), dimnames=list(as.character(ss), as.character(nbins)))
  for(n in seq_along(nbins)){
    for(s in seq_along(ss)){
      rmvec = rep(NA,reps)
      for(rp in 1:reps){ 
        #try({
        j = wfsim(s=ss[s], N=N, start=start, len=gens)
        j = data.frame(j=j, year=1:gens)
        if(n<gens){
          sum.df = j %>% group_by(cut(year, nbins[n])
          ) %>% summarize(value=sum(j), count=n())
          xlin = (sum.df$value / sum.df$count)/N
          #lines(xlin, col=rgb(0,0,0, 0.05))
          rmvec[rp]  = FITfast(xlin )  # FIT function fixes 0,1 events by making them -0.01 from the boundary
        } else { 
          xlin=j$j/N
          rmvec[rp] = FITfast(xlin)
          #lines(xlin, col=rgb(0,0,0, 0.05))
        }
        #})
      }
      fitm[s, n] = mean(rmvec, na.rm=T)
      fitratio[s,n] = length(which(rmvec<0.05))/length(which(!is.na(rmvec)))
    }
    cat(".")
  }
  return(list(fitm=fitm, nx=nx, nbins=nbins, fitratio = fitratio))
}


#### Run ####


# The following calls both return a list with 4 elements:
# [[1]] the matrix of FIT mean p-values, across 200 s-values (rows) and 26 increasing bin lengths/decreasing numbers of bins (columns)
# [[2]] the bin lenghts
# [[3]] the number of bins for each bin length
# [[4]] another matrix: over 100 replicates of each combination, % of FIT p<0.05
# This is the space of s values explored by default: c(0, exp(seq(log(0.01), log(5), length.out = 200-1))

Sys.time() # running these might take a while
fitm1 = dofitm(N=1000, start=500, reps=100, maxs=5, gens=200) # starting at 50%
fitm2 = dofitm(N=1000, start=50, reps=100, maxs=5, gens=200)  # starting at 5%
Sys.time()


#### Some quick plots ####

# Space of selection strength values s explored:
plot( c(0, exp(seq(log(0.01), log(5), length.out = 200-1))), ylab="s", cex=0.5 )

# 50% start, 3 example slices
plot(NA, ylim=c(-0.001,0.3), xlim=c(1,26), xlab="", ylab="",yaxs="i", xaxt="n", yaxt="n")
axis(1, at = 1:26, labels = fitm1$nbins, las=1); axis(2, at = c(0.01,0.05,0.1,0.2,0.3), labels=c(0.01,0.05,"0.1","0.2","0.3"), las=1)
pc=c(1,4,17); p=1 
for(i in c(2,24, 75)){ # exact: 0.01, 0.01815483, 0.09887234
  lines(1:26,fitm1[[1]][i,], type="b", cex=0.8, lwd=1, pch=pc[p])
  p=p+1
}
legend(x=0, y=0.28, legend = c("s=0.1","s=0.02", "s=0.01"), border=NULL, bg=NA,bty = "n",pch=rev(pc))
mtext("number of bins",1, 2.1, cex=0.9); mtext(expression("Mean FIT "~italic(p)*"-value"),2, 2,cex=0.9)

# Full heatmaps
par(mfrow=c(2,2), mar=c(1,1,1,1))
image((fitm1[[4]]), breaks=seq(-0.00001,1, length.out = 1000), col=gray.colors(999, start=1,end=0), xaxt="n",yaxt="n")
image((fitm1[[1]]),breaks=c(-1,0.01,0.05,0.1,0.2,1), col=RColorBrewer::brewer.pal(5, "RdYlBu"), xaxt="n",yaxt="n")
image((fitm2[[4]]), breaks=seq(-0.00001,1, length.out = 1000), col=gray.colors(999, start=1,end=0), xaxt="n",yaxt="n")
image((fitm2[[1]]),breaks=c(-1,0.01,0.05,0.1,0.2,1), col=RColorBrewer::brewer.pal(5, "RdYlBu"), xaxt="n",yaxt="n")
mtext(text="<- low s   high s ->", outer=T,side = 1, line=-1)
mtext(text="<- many short bins   few long bins -> ", outer=T,side = 2, line=-1)


