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

  mif1 %>%
    simulate(nsim=200,as.data.frame=TRUE,include.data=TRUE,na.rm=T) %>%
    subset(select=c(time,sim,deaths)) ->sim.data

  sim.data %>%
    mutate(data=sim=="data")%>%
    ddply(~time+data,summarize,
          p=c(0.05,0.5,0.95),q=quantile(deaths,prob=p,names=FALSE,na.rm=T)) %>%
    mutate(p=mapvalues(p,from=c(0.05,0.5,0.95),to=c("lo","med","hi")),
           data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
    dcast(time+data~p,value.var='q') %>%
    ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
    #geom_line(aes(time,med),linetype=2)+
    geom_ribbon(alpha=0.2)+
    ggtitle(nm)+
    #annotate(geom = 'label',x = 1904,y=Inf,nm)
    theme_bw(base_size = 12)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    ylab(NULL)+xlab(NULL)+theme(legend.position = 'none')+
    geom_vline(xintercept = seq(1898,1906,by=1),alpha=0.4,linetype=2)+
    scale_x_continuous(breaks=seq(1898,1906,by=4))->q

  assign(paste0('q.',nm),q)

}

q.OR.west <- q.OR.west +theme(legend.title = element_text(NULL),legend.position = 'top')+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
  scale_fill_discrete(name=NULL)+
  scale_color_discrete(name=NULL)+
  ylim(0,80)


q.OR.south <- q.OR.south + ylab('deaths')
q.OR.east <- q.OR.east + ylab('deaths')

q.OR.north <- q.OR.north + ggtitle('outer north')
q.OR.west <- q.OR.west + ggtitle('outer west')
q.OR.east <- q.OR.east + ggtitle('outer east')
q.OR.south <- q.OR.south + ggtitle('outer south')

require(cowplot)
p <- plot_grid(q.OR.south, q.OR.west ,q.OR.east , q.OR.north, align = 'hv',labels = 'AUTO')

jpeg(file='London_sine_OR_fits.jpg',height=6,width=8,units="in",res=300)
print(p)
dev.off()
