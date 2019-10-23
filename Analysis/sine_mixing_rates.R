rm(list=ls())

require(Rmpi)
require(doRNG)
require(WaveletComp)
require(forecast)
require(lubridate)
require(pomp)
require(tsiR)
require(plyr)
require(magrittr)
require(reshape2)
require(ggplot2)

args<-commandArgs(TRUE)

inSLURM <- c(Sys.getenv('SLURM_JOB_ID') != "")

if(!inSLURM){
  args <- c('north','south')
  setwd('~/Dropbox/LondonRegionalAnalysis////')

}

r1 <- c('north','north','north','north','south','south','south','east','east','west','north','east','west','south')
r2 <- c('west','east','south','central','west','east','central','west','central','central','OR.north','OR.east','OR.west','OR.south')

for(it in 1:length(r1)){


  file <- paste0('Data/InSilico/',r1[it],'_',r2[it],'_sine_0_0.5_monthly.RData')
  load(file)

  ggplot(mf1,aes(eps,loglik))+geom_point()+
    geom_hline(yintercept = max(mf1$loglik) - 1.92)+
    ggtitle(file)-> p1
  print(p1)

  best <- which.max(mf1$loglik)

  max.ll <- max(mf1$loglik)
  CI.range <- which(max.ll - mf1$loglik  < 1.92)
  eps.range87 <- range(mf1$eps[CI.range])

  nm2 <- r2[it]
  if(nm2 == 'OR.north') nm2 <- 'outer north'
  if(nm2 == 'OR.south') nm2 <- 'outer south'
  if(nm2 == 'OR.east') nm2 <- 'outer east'
  if(nm2 == 'OR.west') nm2 <- 'outer west'

  df <- data.frame(
    p1 = r1[it],
    p2 = nm2,
    eps = mf1$eps[best],
    low87 = eps.range87[1],
    high87 = eps.range87[2],
    sum.iota = mf1$iota1[best] + mf1$iota2[best]
  )

  if(it == 1){

    mixing <- df
  } else{
    mixing <- rbind(mixing,df)
  }
}

mixing$combo <- paste0(mixing$p1,'-',mixing$p2)

ggplot()+
  geom_point(data=mixing,aes(combo,eps))+
  geom_errorbar(data=mixing,aes(combo,ymin=low87,ymax=high87))+
  theme_bw(base_size = 13)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  #theme(legend.position = 'none')+
  #theme(legend.position = c(0.8, 0.8))+
  ylab('coupling rate') + xlab('pairwise region')->
  p2


stop()
#jpeg(file='London_sine_mixing.jpeg',height=6,width=8,units="in",res=300)
print(p2)
#dev.off()

