rm(list=ls())

require(forecast)
require(lubridate)
require(pomp)
require(MASS)
require(tsiR)
require(plyr)
require(magrittr)
require(reshape2)
require(ggplot2)
require(Rmisc)
require(scales)
inSLURM <- c(Sys.getenv('SLURM_JOB_ID') != "")

if(!inSLURM){

  setwd('/Users/adbecker/Dropbox/LondonRegionalAnalysis///')
}

nms <- c('OR.north','OR.south','OR.east','OR.west')

comp.pr <- function(ts,p=F){

  spec <- spectrum(ts,plot=F)
  options <- 1/spec$freq/13
  power <- spec$spec

  if(p){
    plot(options,power,type='l',xlim=c(0,3))
  }

  freq <- options[which.max(power)]
  bi.ind.ll <-which.min(abs(options - 1.75))
  bi.ind.ul <-which.min(abs(options - 2.25))

  an.ind.ll <-which.min(abs(options - 0.75))
  an.ind.ul <-which.min(abs(options - 1.25))

  yr2 <- power[bi.ind.ll:bi.ind.ul]
  yr1 <- power[an.ind.ll:an.ind.ul]

  pr <- sum(yr2)/sum(yr1)

  pr.dat <- data.frame(options,power)

  return(pr.dat)
}

nsim <- 100
for(nm in nms){
  print(nm)
  load(paste0('Data/InSilico/',nm,'_monthly_sine.RData'))

  mf1 <- mf1[order(mf1$loglik),]

  best <- which.max(mf1$loglik)

  parms <- as.numeric(mf1[best,1:length(coef(m1))])
  parm.names <- names(coef(m1))
  names(parms) <- parm.names

  max.ll <- max(mf1$loglik)

  mif1 <-m1
  coef(mif1) <- parms

  data <- data.frame(m1)

  true <- comp.pr(data$deaths)

  pr.vec <- rep(NA,nsim)
  res <- matrix(NA,nrow(true),nsim)


  for(n in 1:nsim){
    print(n/nsim)
    sim <-data.frame(simulate(mif1))

    x <- sim$deaths

    pr.dat <- comp.pr(x)

    res[,n] <- pr.dat$power / sum(pr.dat$power)

  }

  ci <- apply(res, 1, function(x){quantile(x,probs=c(.025,.975))})

  high <- ci[2,]
  low <- ci[1,]
  mean <- apply(res, 1, function(x){mean(x)})
  res <- data.frame(res)

  res$options <- true$options
  res$mean <- mean
  res$low <- low
  res$high <- high


  q <- ggplot()+
    geom_line(data=true,
              aes(options,power/sum(power)),colour='black',size=0.9,alpha=0.8)+xlim(0,3.5)+
    geom_line(data=res,aes(options,mean),colour='red',alpha=0.8,size=1)+
    geom_ribbon(data=res,aes(x = options,ymin=low,ymax=high),alpha=0.3,fill='red')+
    theme_bw()+
    xlab('periodicity (years)')+
    ylab('normalized power')+
    #scale_y_continuous(breaks=seq(0,0.2,length=3))+
    theme_bw(base_size = 12)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position = 'none')+
    ggtitle(nm)

  assign(paste0('q.',nm),q)

}

q.OR.south <- q.OR.south + ylab('normalized power')
q.OR.east <- q.OR.east + ylab('normalized.power')

q.OR.north <- q.OR.north + ggtitle('outer north')
q.OR.west <- q.OR.west + ggtitle('outer west')
q.OR.east <- q.OR.east + ggtitle('outer east')
q.OR.south <- q.OR.south + ggtitle('outer south')

require(cowplot)
p <- plot_grid(q.OR.south, q.OR.west ,q.OR.east , q.OR.north, align = 'hv',labels = 'AUTO')

#jpeg(file='London_sine_OR_spectra_fits.jpg',height=6,width=8,units="in",res=300)
print(p)
#dev.off()
