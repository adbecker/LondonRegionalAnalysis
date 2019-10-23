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
  setwd('/Users/adbecker/Dropbox/LondonRegionalAnalysis//')

}
load('Data/InSilico/sine_region_results.RData')

log.ir <- res %>%
  dplyr::mutate(amplitude = amp,
                phase = phase) %>%
  dplyr::select(CBR,amplitude,phase) %>% log()


# log.ir <- res %>%
#   dplyr::select(-c(CBRlow,CBRhigh,rho,cohort,alpha,pop,
#                    sigma,gamma,mu,loglik,loglik.se,
#                    nfail.min,nfail.max,etime,roundfreq,
#                    R_0,region,R0,period)) %>% log()
#


ir.species <- res$roundfreq

# apply PCA - scale. = TRUE is highly
# advisable, but default is FALSE.
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE)

# print method
print(ir.pca)

summary(ir.pca)

library(devtools)

library(ggbiplot)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1,
              groups = factor(ir.species), ellipse = TRUE,
              circle = TRUE) +
  scale_color_discrete(name = 'periodicity')+
  theme(legend.direction = 'horizontal',
        legend.position = 'top')+
  theme_classic()

#jpeg(file='sine_CBR_PCA.jpg',height=4,width=6,units="in",res=300)

print(g)
#dev.off()


