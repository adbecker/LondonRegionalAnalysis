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
  args <- 'north'
  setwd('/Users/adbecker/Dropbox/LondonRegionalAnalysis//')

}
load('Data/London_Regions_data.RData')

PR.fun <- function(x){
  x[is.na(x)] <- 0
  spec <- spectrum(x,plot=F,na.rm=T)
  options <- 1/spec$freq/13
  power <- spec$spec

  freq <- options[which.max(power)]
  bi.ind.ll <-which.min(abs(options - 1.85))
  bi.ind.ul <-which.min(abs(options - 2.25))

  an.ind.ll <-which.min(abs(options - 0.85))
  an.ind.ul <-which.min(abs(options - 1.25))

  yr2 <- power[bi.ind.ll:bi.ind.ul]
  yr1 <- power[an.ind.ll:an.ind.ul]

  pr <- sum(yr2)/sum(yr1)
}

period.fun <- function(x){
  x[is.na(x)] <- 0
  spec <- spectrum(x,plot=F,na.rm=T)
  options <- 1/spec$freq/13
  power <- spec$spec

  period <- round(options[which.max(power)])

  return(period)
}


regions <- c('all','south','west','central','east','north','OR.west','OR.east','OR.south','OR.north')

contact.df <- data.frame(matrix(NA,13,length(regions)))
names(contact.df) <- regions

for(args in regions){

  rdata <- subset(data,scale == args)

  load(paste0('Data/InSilico//',args,'_monthly_sine.RData'))

  mf1 <- mf1[order(mf1$loglik),]

  best <- which.max(mf1$loglik)

  mif1 <-m1
  coef(mif1) <- mf1[best,]

  sim <- data.frame(simulate(mif1))

  sim %>%
    subset(time >= 1900 & time < 1901) %>%
    dplyr::select(contact)->beta

  contact.df[[as.character(args)]] <- beta

  new.df <- mf1[best,]
  new.df$region <- args

  period <- period.fun(rdata$deaths)
  new.df$period <- period

  rdata %>%
    dplyr::mutate(year = floor(time)) %>%
    #subset(year != 1897) %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(CBR = mean(birthrate) * 1000 / mean(pop)) -> r.CBR


  new.df$pop <- mean(rdata$pop)
  new.df$CBR <- mean(r.CBR$CBR)
  new.df$CBRlow <- range(r.CBR$CBR)[1]
  new.df$CBRhigh <- range(r.CBR$CBR)[2]

  if(args == regions[1]){
    df <- new.df
  }else{
    df <- rbind(df,new.df)
  }

}

stop()

df$R0 <- df$beta0 / df$gamma

res <- df
res$roundfreq <- res$period


log.ir <- res %>%
  dplyr::mutate(amplitude = amp,
                phase = phase) %>%
  dplyr::select(CBR,amplitude,phase) %>% log()

# log.ir <- res %>%
#   dplyr::select(-c(CBRlow,CBRhigh,rho,cohort,alpha,nzeros,freq,PR,year,sigma,gamma,mu,loglik,loglik.se,nfail.min,nfail.max,etime,
#                    district,roundfreq,I0.high,I0.low,S0.high,S0.low,cos.amp,cos.phase,amp,lon,lat,R_0)) %>% log()

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

g

save(res,file='Data/InSilico/sine_region_results.RData')

grid.size <- 150

source('Analysis/phase_amp_grid.R')

stop()

grid.file <- 'Data/InSilico/sine_bifurcation_092319.RData'
if(file.exists(grid.file)){
  load(grid.file)
}else{
  d <- phase_amp_grid(stochastic = FALSE ,grid.size = grid.size,R0=round(mean(res$R0[-1])),
                      phase.min=0,phase.max=2,amp.min=0.03,amp.max=.35)

  save(d,file = grid.file)
}

d$roundfreq <- as.character(d$roundfreq)
d$roundfreq[which(d$roundfreq == '3' | d$roundfreq == '3+ ')] = '2+'

res$alpha  <- c(rep(T,6),rep(F,4))

OR.alpha <- 1

p1 <-  ggplot(data=d, aes(amp, phase)) + geom_tile(aes(fill = roundfreq),colour='white') +
  geom_point(data=subset(res,region == 'all'),aes(amp,phase,fill=as.factor(roundfreq),alpha=alpha),
             size=6,colour='grey90',pch=22,stroke=2)+
  geom_point(data=subset(res,region != 'all'),aes(amp,phase,fill=as.factor(roundfreq),alpha=alpha),
             size=6,colour='grey90',pch=21,stroke=2)+
  scale_alpha_discrete(range = c(OR.alpha, 1),guide=FALSE)+
  scale_colour_manual(values = c("black", "dodgerblue",'grey75'),
                      name= "period", guide = guide_legend(reverse = F))+
  scale_fill_manual(values = c("black", "dodgerblue",'grey75'),
                    name= "period", guide = guide_legend(reverse = F))+
  theme(text = element_text(size=10))+
  xlim(c(0.03,.3))+
  ylim(range(d$phase))+
  xlab(expression(seasonal~forcing~'amplitude,'~alpha))+
  ylab(expression(seasonal~forcing~'phase,'~phi))+
  theme_classic()+
  theme(legend.position="top",
        legend.text = element_text(size=10),
        legend.key.size = unit(3,'point'),
        legend.margin  = unit(-2,'cm'))


p1


date_num <- as.numeric(sim$time)
year <- floor(date_num)
year_beginning <- as.POSIXct(paste0(year, '-01-01'))
year_end <- as.POSIXct(paste0(year+1, '-01-01'))
date <- year_beginning + (date_num %% 1) * (year_end - year_beginning)
months <- month.abb[month(format(date, format='%Y-%m-%d'))]
month.letter <- substring(months, 1, 1)

scaled.contact <- data.frame(apply(contact.df,2,function(x) x / mean(x,na.rm=T)))
scaled.contact$time <- seq(1,13,1)

region.names <- as.character(res$region)
region.names[7:10] <- c('outer west', 'outer east', 'outer south','outer north')

melt.contact <- melt(scaled.contact,id.vars = 'time')
melt.contact$variable <- rep(region.names,each=13)

melt.contact$period <- rep(res$roundfreq,each = 13)

melt.contact$roundfreq <- rep(res$roundfreq,each = 13)
melt.contact$ID <- rep(res$ID,each=13)
melt.contact$region <- paste0(melt.contact$ID,': ',melt.contact$variable)
melt.contact$roundfreq <- melt.contact$region

melt.contact$alpha <- rep(c(rep(T,6),rep(F,4)),each=13)

label.region.names <- names(table(melt.contact$region))

std <- function(x) sd(x)/sqrt(length(x))

melt.contact %>%
  subset(variable != 'all')%>%
  dplyr::select(time,value,variable,period) %>%
  dplyr::group_by(time,period) %>%
  dplyr::summarize(mean = mean(value),
                   low = mean(value) - std(value),
                   high =mean(value) + std(value)) -> summary.contact


require(ggrepel)
ggplot()+
  geom_point(data = subset(melt.contact,variable != 'all'),aes(time,value,colour=factor(period),shape=factor(period),group=variable,alpha = 1))+
  geom_point(data = subset(melt.contact,variable == 'all'),aes(time,value,colour=factor(period),group=variable,alpha = 1),pch=15,stroke=2)+
  geom_line(data = summary.contact,aes(time,mean,color = factor(period)),size=1)+
  # geom_errorbar(data = summary.contact, aes(time,ymin=low,ymax=high,color = period),alpha=0.8)+
  geom_ribbon(data = summary.contact,
              aes(time,ymin=low,ymax=high,color = factor(period),fill = factor(period)),alpha=.8)+
  #scale_alpha_discrete(range = c(OR.alpha, 1),guide=FALSE)+
  #scale_x_continuous(breaks=1:13)+
  scale_x_continuous(breaks = 1:12,
                     labels = month.letter[1:12])+
  #geom_text_repel(data=subset(melt.contact,time == 1),
  #         aes(label=ID),colour='black',size=4)+
  #  geom_text_repel(data=subset(melt.contact,time == 13),
  #           aes(label=ID),colour='black',size=4)+
  #annotate('text',x=9,y=1.3,label='1: south\n2: west\n3: central\n4: west\n5: north',hjust=0,size=2)+
  scale_colour_manual(
    name= NULL,
    labels = NULL,
    #values = c("dodgerblue", rep("black",2),rep('dodgerblue',3),'black','dodgerblue','black','dodgerblue'),
    values = c("black", rep("dodgerblue",2),rep('black',3),'dodgerblue','black','dodgerblue','black'),
    guide = guide_legend(reverse = F)
  )+
  scale_fill_manual(
    name= NULL,
    labels = NULL,
    #values = c("dodgerblue", rep("black",2),rep('dodgerblue',3),'black','dodgerblue','black','dodgerblue'),
    values = c("black", rep("dodgerblue",2),rep('black',3),'dodgerblue','black','dodgerblue','black'),
    guide = guide_legend(reverse = F)
  )+
  # scale_linetype_manual(
  #   name= "region",
  #   labels=label.region.names,
  #   values= 1:10,
  #   guide = guide_legend(reverse = F)
  # )+
  #guides(colour=FALSE)+
  #ylab(bquote(normalized~'seasonality,'~beta~'/mean('~beta~')'))+
  #ylab(bquote(beta~'/mean('~beta~')'))+
  ylab('normalized seasonality')+
  xlab('month')+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size=10))->
  p3

p3

melt.contact$int <- paste(melt.contact$roundfreq, melt.contact$region, sep=".")

stop()
require(cowplot)

jpeg(file='sine_OR_no_labels_phase_amp.jpeg',height=8,width=8,units="in",res=300)
#png("ts_plots.png",width=12,height=9,units="in",res=300)
#pdf(file='OR_phase_amp.pdf',height=8,width=8)

plot_grid(p3, p1, labels = c("A", "B"),rel_heights  = c(1,2),nrow=2,align= 'v',axis='l')

dev.off()

