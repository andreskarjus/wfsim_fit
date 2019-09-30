# Code to replicate the simulations described in:
# Karjus et al. 2018, "Challenges in detecting evolutionary forces in language change using diachronic corpora" (ArXiv e-print)
# (run Wright-Fisher simulations with given parameters, bin the 
# output into bins akin to those used in diachronic text corpora,
# and apply the FIT (Feder et al. 2014), following Newberry et al. (2017))
# 
# Code: Andres Karjus

# Required packages
library(dplyr)    
library(nortest)
library(compiler)     # recommended: speeds up loops using JIT compilation
enableJIT(3)



FITfast = function(af.vec, alltests=T){
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
  if(alltests){
    if(length(Yi)>4){ kolm = lillie.test(Yi)$p.value} else {kolm=NA} # Lilliefors (Kolmogorov-Smirnov) test for the composite hypothesis of normality
    if(length(Yi)>7){ad = ad.test(Yi)$p.value} else { ad=NA } # Anderson-Darling test, min 8
  } else {
    kolm=0;ad=0
  }
  
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
}; wfsim=cmpfun(wfsim,options=list(optimize=3))


dofitm = function(N=1000, start=500, reps=1000, maxs=5, gens=200,ns=200, minshare=10){
  ss = c(0,exp(seq(log(0.001), log(maxs), length.out = ns-1)))
  maxlen = round(median(table(cut(1:gens, 4)))) # length of bin if minimal n bins
  nx= (1:maxlen)[!duplicated((sapply(1:maxlen, function(x) length(seq(1,gens,x)))))]   # bin lengths, 25
  nbins = sapply(nx, function(x) length(seq(1,gens, x) )) # 25
  minper = reps/minshare
  
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
        tmp = FITfast(xlin, alltests = F)  
        # FIT function fixes 0,1 events by making them +-0.001 from the boundary (with N=1000, it's basically +-1)
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
      if(length(which(!is.na(rmvec))) < minper ){  # if too few valid reps, then don't calculate % & mean
        fitswm[s,n] = NA
        fitratioswm[s,n] = NA
      } else {
        fitswm[s,n] = mean(rmvec, na.rm=T)   #  mean FIT p-value, but sw-filtered (so if all NA then NA!)
        fitratioswm[s,n] = length(which(rmvec<0.05))/length(which(!is.na(rmvec)))
      }
    }
    cat(".")
  }
  return(list(fitm=fitm, nx=nx, nbins=nbins, fitratio = fitratio, swm=swm, fitswm = fitswm,fitratioswm=fitratioswm, kolmm=kolmm, adm=adm))
}



#### Time series length and selection strength for appendix ####

dofitm_lengths = function(N=1000, start=500, reps=1000, maxs=5, lens=4:200,ns=200, minshare=10){
  ss = c(0,exp(seq(log(0.001), log(maxs), length.out = ns-1)))
  minper=reps/minshare
  fitratio = matrix(NA, length(ss), length(lens), dimnames=list(as.character(ss), as.character(lens)))
  fitratioswm = matrix(NA, length(ss), length(lens), dimnames=list(as.character(ss), as.character(lens)))
  
  for(n in seq_along(lens)){
    for(s in seq_along(ss)){
      rmvec = rep(NA,reps)
      swvec = rep(NA,reps)
      for(rp in 1:reps){ 
        j = wfsim(s=ss[s], N=N, start=start, len=lens[n]) # difference here
        try({ # once got  Error in shapiro.test(Yi) : all 'x' values are identical for 4-length
          tmp = FITfast(j/N, alltests = F)  
          rmvec[rp]  = tmp$fitp
          swvec[rp]  = tmp$swp
        })
      }
      
      fitratio[s,n] = length(which(rmvec<0.05))/length(which(!is.na(rmvec)))
      # filter:
      rmvec[swvec<0.1] = NA # set to NA if normality assumption violated
      # record after filtering:
      if(length(which(!is.na(rmvec))) < minper ){  # if too few valid reps, then don't calculate % & mean
        fitratioswm[s,n] = NA
      } else {
        fitratioswm[s,n] = length(which(rmvec<0.05))/length(which(!is.na(rmvec)))
      }
    }
    cat(".")
  }
  return(list(lens=lens, fitratio = fitratio, fitratioswm=fitratioswm))
}



#### run ####
# running these might take a while

Sys.time()
fitm1 = dofitm(N=1000, start=500, reps=1000, maxs=5, gens = 200, ns = 200)
Sys.time()
fitm2 = dofitm(N=1000, start=50, reps=1000, maxs=5, gens=200, ns=200)
Sys.time()



#### run length model ####
lens=unique(round(exp(seq(log(4), log(200), length.out = 200))))
Sys.time()
fitm1_len = dofitm_lengths(N=1000, start=500, reps=1000, maxs=5, lens = lens, ns = 100)
Sys.time()
fitm2_len = dofitm_lengths(N=1000, start=50, reps=1000, maxs=5, lens= lens, ns=100)
Sys.time()


# The calls both return a list with these elements:
# [[1]] the matrix of FIT mean p-values, across 200 s-values (rows) and 26 increasing bin lengths/decreasing numbers of bins (columns)
# [[2]] the bin lenghts
# [[3]] the number of bins for each bin length
# [[4]] another matrix: over 1000 replicates of each combination, % of FIT p<0.05
# [[5]] Shapiro-Wilk p
# [[6]] same as 1, but with SW p<0.1 cases removed before calculating mean
# [[7]] same as 1, but with SW p<0.1 cases removed before calculating %

# This is the space of s values explored by default: c(0, exp(seq(log(0.001), log(5), length.out = 200-1))

# N = pop size; start = mutant starting n; reps = how many replications of each combo; maxs = largest s value to test; gens=how many generations (~years in a corpus); ns = how many s values to test (taken from a log sequence by default, starting with log(0.001), with s=0 prepended to the start)




#### plotting  ####

library(ggplot2)
library(viridisLite)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(reshape2)


thm = theme(
  panel.border = element_rect(color="darkgray",fill=NA), 
  legend.background = element_rect(fill="white", color=NA),
  legend.box.spacing = unit(0,"mm"),
  legend.box.margin = margin(0,0,0,0),
  axis.ticks= element_line(color="darkgray", size = 0.1), 
  axis.ticks.length = unit(-1.4, "mm"), 
  axis.text.x = element_text(margin=unit(c(t = 2, r = 0, b = 0, l = 0),"mm")),
  axis.text.y.left = element_text(margin=unit(c(t = 0, r = 2, b = 0, l = 0),"mm"))
)


figure4 = function(fitm2, swm=F, ratio=F, lenplot=F, legend=T, theme=thm, lines=NULL) {
  sx=as.numeric(rownames(fitm2$fitratio)); 
  sb = sapply(c(0,0.01,0.1,1:5), function(x) mean(which(abs(sx-x)==min(abs(sx-x)))) )
  atsmall=sort(union(sb,sapply(c(
    seq(0,0.01,0.001),
    seq(0.01,0.1,0.01),
    seq(0.1,1,0.1)
  ), function(x) which(abs(sx-x)==min(abs(sx-x)))[1] )))
  sl=rep("",length(atsmall));sl[which(atsmall %in%sb)]=c(0,0.01,0.1,1:5)
  if(legend) lp = c(0.98,0.98)
  if(!legend) lp = "none"
  tickcol="black"
  ticksize=0.1
  rectcol = rgb(0.95,0.95,0.95,0.5)
  
  thm = theme + theme(
    axis.ticks= element_line(color=tickcol, size = 0.1),
    panel.grid = element_line(color="grey90", size=0.1),
    legend.position = lp, 
    legend.justification = c(1,1), 
    legend.background = element_rect(fill=rgb(1,1,1,0.8), color=NA), 
    legend.margin=margin(c(0,0,0,0)) )
  
  ylb = ifelse(lenplot, "length (n generations)", "number of bins")
  
  if(ratio){
    if(swm) d = fitm2$fitratioswm
    if(!swm) d = fitm2$fitratio
    slabs = 1:nrow(d)
    ylabslen=1:ncol(d)
    
    cols = c(viridis(5, end=0.4)[1], plasma(5, begin=0.5,end = 1, direction = 1), viridis(5, begin=0.55,end=1,  direction=-1)) #; plot(1:11, col=cols, pch=20, cex=5)
    labs=paste0("\u2264", c(0.05,seq(0.1,1, 0.1))*100)
    d = d %>% reshape2::melt(varnames=c("s", "nbins")) %>%  
      mutate(value=cut(value, breaks=c(0,0.05, seq(0.1,1, 0.1)),labels=labs,include.lowest = T)) %>%  
      mutate(s=factor(s, labels=slabs))
    if(!lenplot){
      d = d %>% mutate(nbins=factor(nbins, levels=rev(levels(as.factor(nbins)))))
    } else {
      d = d %>% mutate(nbins=factor(nbins, labels = ylabslen))
    }
    
    fig = ggplot(d, aes(y = nbins, x = s, fill = value)) +
      geom_tile(color=NA)+
      theme_minimal() +
      scale_fill_manual(values= cols, na.value=NA, drop=F, name= "% p<0.05", labels=labs, na.translate = F) +
      scale_x_discrete(expand=c(0,0), breaks=atsmall, labels=sl) +
      labs(x="selection strength parameter",  y=ylb) 
    #geom_vline(xintercept = 0.01, color="black") +
    
    if(!lenplot){
      fig = fig +
        scale_y_discrete(expand=c(0,0), position="left") +
        #geom_rect(mapping =  aes(xmin=0.9,xmax=98,  ymin=16,ymax=23), inherit.aes = F, color=rgb(1,1,1,0.5), fill=NA) +
        # white rectangle:
        annotate("segment", x=0.8, xend=98, y=16-0.4, yend=16-0.4, color=rectcol)+
        annotate("segment", x=0.8, xend=98, y=23+0.4, yend=23+0.4, color=rectcol) +
        annotate("segment", x=98.4, xend=98.4, y=16-0.2, yend=23+0.2, color=rectcol) +
        
        # ticks
        annotate("segment", x=199, y=1:25, yend=1:25,xend=200.5, color=tickcol, size=ticksize ) +        # right side ticks
        annotate("segment", x=atsmall, y=25.5, yend=25.2,xend=atsmall, color=tickcol, size=ticksize ) +  # topside ticks
        annotate("segment", x=atsmall[which(sl!="")], y=0.5, yend=1.1,xend=atsmall[which(sl!="")] , color=tickcol, size=ticksize )  # logticks
      if(!is.null(lines)) fig=fig+annotate("segment", x=c(lines), xend=c(lines), y=1.3, yend=25, color=rgb(0,0,0,0.4))  # example lines
    }
    if(lenplot){
      lx = unique(round(exp(seq(log(4), log(200), length.out = 200))))
      atsmall2 = sapply(c(
        4:10,
        seq(20,200,10)
      ), function(x) (which(abs(lx-x)==min(abs(lx-x))) )[1])
      atbig = sapply(c(4, 10,20,50,100,200), function(x) (which(abs(lx-x)==min(abs(lx-x))) )[1])
      
      fig=fig+
        annotate("segment", x=100, y=atsmall2, yend=atsmall2,xend=100.5, color=tickcol, size=ticksize) +        # right side ticks
        annotate("segment", x=atsmall, y=118.5, yend=118.2,xend=atsmall, color=tickcol, size=ticksize) +  # topside ticks
        annotate("segment", x=0.5, y=atsmall2, yend=atsmall2,xend=1.2, color=tickcol, size=ticksize ) +  # left logticks
        annotate("segment", x=0.5, y=atbig, yend=atbig,xend=1.8, color=tickcol, size=ticksize ) +  # left big logticks
        annotate("segment", x=99.5, y=atbig, yend=atbig,xend=100.5, color=tickcol, size=ticksize ) +  # right big logticks
        annotate("segment", x=atsmall[which(sl!="")], y=0.5, yend=2.1,xend=atsmall[which(sl!="")], color=tickcol, size=ticksize ) +  # logticks
        geom_vline(xintercept = 55 ,alpha=0.5, ) + geom_hline(yintercept = 47, alpha=0.5) +
        scale_y_discrete(expand=c(0,0), position="left", breaks=(atbig), labels = c(4, 10,20,50,100,200)  )
      #scale_y_discrete(expand=c(0,0), position="left", breaks=1, labels = "!!!!!!!!!!!") 
      #+ annotate("text", x=2,y=1:47,label=1:47, size=3, color="white")
    }
    fig=fig+thm
  }
  
  if(!ratio){
    if(swm) d = fitm2$fitswm
    if(!swm) d = fitm2$fitm
    
    cols = RColorBrewer::brewer.pal(5, "RdYlBu")
    labs=c("<0.01", "<0.05", "<0.1", "<0.2", "\u2265 0.2")
    fig =
      d %>% reshape2::melt(varnames=c("s", "nbins")) %>% mutate(value=cut(value, breaks=c(0,0.01,0.05, 0.1,0.2,1),include.lowest = T, right=F)) %>% 
      mutate(s=factor(s, labels=1:200)) %>% 
      mutate(nbins=factor(ceiling(200/nbins))) %>%  # levels=rev(levels(as.factor(nbins))))) %>%  # right-excluding for <0.05 etc
      
      ggplot(aes(y = nbins, x = s, fill = value)) +
      geom_tile(color=NA)+
      theme_minimal() +
      scale_fill_manual(values= cols, na.value=NA, drop=F, name= "mean p", labels=labs,na.translate = F) +
      scale_y_discrete(expand=c(0,0), position="right") + # sec.axis = sec_axis(~as.numeric(.)/200, name="bin length"))+
      scale_x_discrete(expand=c(0,0), breaks=atsmall, labels=sl)+
      theme(axis.ticks= element_line(), 
            axis.ticks.length = unit(-1.4, "mm"), 
            axis.text.x = element_text(margin=unit(c(t = 2, r = 0, b = 0, l = 0),"mm")),
            axis.text.y.right = element_text(margin=unit(c(t = 0, r = 1, b = 0, l = 2),"mm"))
      ) + thm +
      labs(x="selection strength parameter",  y="bin length") +
      
      #geom_rect(mapping =  aes(xmin=0.9,xmax=98,  ymin=16,ymax=23), inherit.aes = F, color=rgb(1,1,1,0.5), fill=NA) +
      # white rectangle:
      annotate("segment", x=0.8, xend=98, y=16-0.4, yend=16-0.4, color=rectcol)+
      annotate("segment", x=0.8, xend=98, y=23+0.4, yend=23+0.4, color=rectcol) +
      annotate("segment", x=98.4, xend=98.4, y=16-0.2, yend=23+0.2, color=rectcol) +
      
      # ticks
      annotate("segment", x=0,5, y=1:25, yend=1:25,xend=1.5, color=tickcol, size=ticksize ) + # ! left side ticks
      annotate("segment", x=atsmall, y=25.5, yend=25.2,xend=atsmall, color=tickcol, size=ticksize ) +  # topside
      annotate("segment", x=atsmall[which(sl!="")], y=0.5, yend=1.1,xend=atsmall[which(sl!="")] , 
               color=tickcol, size=ticksize )  # logticks
    if(!is.null(lines)) fig=fig+annotate("segment", x=c(lines), xend=c(lines), y=1.3, yend=25, color=rgb(0,0,0,0.4))  # example lines
  }
  return(fig)
}



# figure 4 (50%)
plot_grid(figure4(fitm1, swm=F, ratio=T, lines=56),
          figure4(fitm1, swm=F, ratio=F, lines=56),
          figure4(fitm1, swm=T, ratio=T, lines=56),
          figure4(fitm1, swm=T, ratio=F, lines=56), 
          ncol=2, labels = c("a.1", "b.1", "a.2", "b.2"), label_y = 0.99, label_x = -0.008 )

# figure 6 (4.2) 5%
plot_grid(figure4(fitm2, swm=F, ratio=T, lines=88),
          figure4(fitm2, swm=F, ratio=F, lines=88),
          figure4(fitm2, swm=T, ratio=T, lines=88),
          figure4(fitm2, swm=T, ratio=F, lines=88), 
          ncol=2, labels = c("a.1", "b.1", "a.2", "b.2"), label_y = 0.99,  label_x = -0.008 )

# figure, S supplement for lengths
plot_grid(figure4(fitm2_len, swm=F, ratio=T, lenplot = T, legend=F),
          figure4(fitm1_len, swm=F, ratio=T,lenplot = T, legend=F),
          figure4(fitm2_len, swm=T, ratio=T,lenplot = T,  legend=F),
          figure4(fitm1_len, swm=T, ratio=T, lenplot = T, legend=T), 
          ncol=2, labels = c("a.1 (5%) ", "b.1 (50%)", "a.2 (5%) ", "b.2 (50%)"), label_y = 0.99, label_x = 0.08 )






