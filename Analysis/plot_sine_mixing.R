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

make.pomp.object <- function(args){

  save.name <- paste0('Data/InSilico/',args[1],'_',args[2],'_monthly_sine.RData')

  data %>%
    subset(scale == args[1])  %>%
    mutate(pop1 = pop, birthrate1 = birthrate,deaths1=deaths) -> data1

  data %>%
    subset(scale == args[2])  %>%
    mutate(pop2 = pop, birthrate2 = birthrate,deaths2=deaths) -> data2

  ## ----rprocess------------------------------------------------------------
  rproc <- Csnippet("
                    // PATCH ONE

                    double beta1, beta2, br1, seas1, seas2, foi1, dw1, births1;
                    double rate1[7], trans1[7];

                    // cohort effect
                    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
                    br1 = cohort*birthrate1/dt + (1-cohort)*birthrate1;
                    else
                    br1 = (1.0-cohort)*birthrate1;

                    seas1 = beta10*(1 + amp1*sin(2*M_PI*t + phase1));
                    seas2 = beta20*(1 + amp2*sin(2*M_PI*t + phase2));


                    // transmission rate
                    beta1 = seas1;

                    // transmission rate
                    beta2 = seas2;

                    // force of infection
                    foi1 = beta1*((I1+iota1) * (1 - eps) + eps*I2)/pop1 ;
                    // white noise (extra-demographic stochasticity)
                    dw1 = rgammawn(sigmaSE1,dt);

                    rate1[0] = foi1*dw1/dt;  //infection rate (stochastic)
                    rate1[1] = mu;  		    // natural S death
                    rate1[2] = sigma;		  // rate of ending of latent stage
                    rate1[3] = mu;			    // natural E death
                    rate1[4] = gamma;		  // recovery
                    rate1[5] = mu;			    // natural I death
                    rate1[6] = phi1;

                    // Poisson births
                    births1 = rpois(br1*dt);

                    // transitions between classes
                    reulermultinom(2,S1,&rate1[0],dt,&trans1[0]);
                    reulermultinom(2,E1,&rate1[2],dt,&trans1[2]);
                    reulermultinom(2,I1,&rate1[4],dt,&trans1[4]);
                    reulermultinom(1,C1,&rate1[6],dt,&trans1[6]);

                    S1 += births1 - trans1[0] - trans1[1];
                    E1 += trans1[0] - trans1[2] - trans1[3];
                    I1 += trans1[2] - trans1[4] - trans1[5];
                    W1 += (dw1 - dt)/sigmaSE1;  // standardized i.i.d. white noise
                    C1 += trans1[4];           // true incidence
                    D1 = rpois(phi1*C1);
                    contact1 = beta1;

                    // PATCH TWO

                    double br2, foi2, dw2, births2;
                    double rate2[7], trans2[7];

                    // cohort effect
                    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
                    br2 = cohort*birthrate2/dt + (1-cohort)*birthrate2;
                    else
                    br2 = (1.0-cohort)*birthrate2;

                    // force of infection
                    foi2 = beta2*((I2+iota2) * (1 - eps) + eps*I1)/pop2 ;

                    // white noise (extra-demographic stochasticity)
                    dw2 = rgammawn(sigmaSE2,dt);

                    rate2[0] = foi2*dw2/dt;  //infection rate (stochastic)
                    rate2[1] = mu;  		    // natural S death
                    rate2[2] = sigma;		  // rate of ending of latent stage
                    rate2[3] = mu;			    // natural E death
                    rate2[4] = gamma;		  // recovery
                    rate2[5] = mu;			    // natural I death
                    rate2[6] = phi2;

                    // Poisson births
                    births2 = rpois(br2*dt);

                    // transitions between classes
                    reulermultinom(2,S2,&rate2[0],dt,&trans2[0]);
                    reulermultinom(2,E2,&rate2[2],dt,&trans2[2]);
                    reulermultinom(2,I2,&rate2[4],dt,&trans2[4]);
                    reulermultinom(1,C2,&rate2[6],dt,&trans2[6]);

                    S2 += births2 - trans2[0] - trans2[1];
                    E2 += trans2[0] - trans2[2] - trans2[3];
                    I2 += trans2[2] - trans2[4] - trans2[5];
                    W2 += (dw2 - dt)/sigmaSE2;  // standardized i.i.d. white noise
                    C2 += trans2[4];           // true incidence
                    D2 = rpois(phi2*C2);
                    contact2 = beta2;

                    ")
  ## ----initializer---------------------------------------------------------
  initlz <- Csnippet("
                     double m1 = pop1/(S1_0+E1_0+I1_0+R1_0);
                     double m2 = pop2/(S2_0+E2_0+I2_0+R2_0);

                     S1 = nearbyint(m1*S1_0);
                     E1 = nearbyint(m1*E1_0);
                     I1 = nearbyint(m1*I1_0);
                     W1 = 0;
                     C1 = 0;
                     D1 = 0;


                     S2 = nearbyint(m2*S2_0);
                     E2 = nearbyint(m2*E2_0);
                     I2 = nearbyint(m2*I2_0);
                     W2 = 0;
                     C2 = 0;
                     D2 = 0;

                     ")

  ## ----dmeasure------------------------------------------------------------
  dmeas <- Csnippet("
                    double m1 = rho*D1;
                    double v1 = m1*(1.0-rho+psi1*psi1*m1);
                    double tol = 1.0e-18;
                    double m2 = rho*D2;
                    double v2 = m2*(1.0-rho+psi2*psi2*m2);
                    double lik1, lik2;

                    if (deaths1 > 0.0) {
                    lik1 = pnorm(deaths1+0.5,m1,sqrt(v1)+tol,1,0)-pnorm(deaths1-0.5,m1,sqrt(v1)+tol,1,0)+tol;
                    } else {
                    lik1 = pnorm(deaths1+0.5,m1,sqrt(v1)+tol,1,0)+tol;
                    }
                    if(ISNA(deaths1)) {
                    lik1 = (give_log) ? 0 : 1;
                    }


                    if (deaths1 > 0.0) {
                    lik2 = pnorm(deaths2+0.5,m2,sqrt(v2)+tol,1,0)-pnorm(deaths2-0.5,m2,sqrt(v2)+tol,1,0)+tol;
                    } else {
                    lik2 = pnorm(deaths2+0.5,m2,sqrt(v2)+tol,1,0)+tol;
                    }
                    if(ISNA(deaths2)) {
                    lik2 = (give_log) ? 0 : 1;
                    }

                    lik = lik1 + lik2;

                    ")

  ## ----rmeasure------------------------------------------------------------
  rmeas <- Csnippet("
                    double m1 = rho*D1;
                    double v1 = m1*(1.0-rho+psi1*psi1*m1);
                    double tol = 1.0e-18;
                    deaths1 = rnorm(m1,sqrt(v1)+tol);

                    if (deaths1 > 0.0) {
                    deaths1 = nearbyint(deaths1);
                    } else {
                    deaths1 = 0;
                    }

                    double m2 = rho*D2;
                    double v2 = m2*(1.0-rho+psi2*psi2*m2);
                    deaths2 = rnorm(m2,sqrt(v2)+tol);

                    if (deaths2 > 0.0) {
                    deaths2 = nearbyint(deaths2);
                    } else {
                    deaths2 = 0;
                    }
                    ")

  ## ----transforms----------------------------------------------------------
  toEst <- Csnippet("
                    Teps = logit(eps);
                    ")

  fromEst <- Csnippet("
                      Teps = expit(eps);
                      ")



  data1 %>%
    dplyr::select(time,pop1,birthrate1) -> covar1

  data2 %>%
    dplyr::select(time,pop2,birthrate2) -> covar2


  covar <- cbind(covar1,covar2)
  covar <- covar[,-4]


  infdata <- data.frame('time'=data1$time,'deaths1'=data1$deaths1,'deaths2'=data2$deaths2)

  paramnames=c('phi1','phi2',
               "sigma","gamma","alpha",
               "iota1",'iota2',
               "rho",
               "sigmaSE1",'sigmaSE2',
               "psi1",'psi2',
               "cohort",
               "S1_0","E1_0","I1_0","R1_0",
               "S2_0","E2_0","I2_0","R2_0",
               'mu','eps',
               'beta10','amp1','phase1',
               'beta20','amp2','phase2'

  )
  load(paste0('Data/InSilico/',args[1],'_monthly_sine.RData'))
  mf1 <- mf1[order(mf1$loglik),]
  region1.mf1 <- mf1

  r1.parms <- region1.mf1[which.max(region1.mf1$loglik),]

  load(paste0('Data/InSilico/',args[2],'_monthly_sine.RData'))
  mf1 <- mf1[order(mf1$loglik),]
  region2.mf1 <- mf1

  r2.parms <- region2.mf1[which.max(region2.mf1$loglik),]
  ## ----pomp-construction---------------------------------------------------
  delta.t <- 1/365.25
  infdata %>%
    pomp(t0=with(infdata,2*time[1]-time[2]),
         time="time",
         rprocess=euler.sim(rproc,delta.t=delta.t),
         dmeasure=dmeas,
         rmeasure=rmeas,
         toEstimationScale=toEst,
         fromEstimationScale=fromEst,
         initializer=initlz,
         covar=covar,
         tcovar="time",
         zeronames=c("C1","W1","D1",'contact1',
                     "C2","W2","D2",'contact2'),
         statenames=c("S1","E1","I1","C1","W1",'D1','contact1',
                      "S2","E2","I2","C2","W2",'D2','contact2'),
         paramnames=paramnames
    ) -> m1

  parms <- c('phi1' = r1.parms$phi,'phi2'=r2.parms$phi,
             "sigma"=r1.parms$sigma,"gamma"=r1.parms$gamma,"alpha"=r1.parms$alpha,
             "iota1"=r1.parms$iota,'iota2'=r2.parms$iota,
             "rho"=r1.parms$rho,
             "sigmaSE1"=r1.parms$sigmaSE,'sigmaSE2'=r2.parms$sigmaSE,
             "psi1"=r1.parms$psi,'psi2'=r2.parms$psi,
             "cohort"=r1.parms$cohort,
             "S1_0"=r1.parms$S_0,"E1_0"=r1.parms$E_0,"I1_0"=r1.parms$I_0,"R1_0"=r1.parms$R_0,
             "S2_0"=r2.parms$S_0,"E2_0"=r2.parms$E_0,"I2_0"=r2.parms$I_0,"R2_0"=r2.parms$R_0,
             'mu'=r2.parms$mu,'eps'=0.1,
             'beta10' = r1.parms$beta0,'amp1'=r1.parms$amp,'phase1'=r1.parms$phase,
             'beta20' = r2.parms$beta0,'amp2'=r2.parms$amp,'phase2'=r2.parms$phase

  )

  names(parms) = paramnames

  m1 <- pomp(m1,params = parms)


  return(m1)
}

specfun <- function(ts){

  spec <- spectrum(ts,plot=F)
  options <- 1/spec$freq/ 13
  power <- spec$spec

  bi.ind.ll <-which.min(abs(options - 1.75))
  bi.ind.ul <-which.min(abs(options - 2.25))

  an.ind.ll <-which.min(abs(options - 0.75))
  an.ind.ul <-which.min(abs(options - 1.25))

  yr2 <- power[bi.ind.ll:bi.ind.ul]
  yr1 <- power[an.ind.ll:an.ind.ul]

  pr <- sum(yr2)/sum(yr1)

  freq <-  mean(power[bi.ind.ll:bi.ind.ul] / power[an.ind.ll:an.ind.ul])

  freq <- options[which.max(power)]
  return(c(freq))
}


make.mixing <- function(m1){

  eps.vec <- seq(0,0.5,by=0.02)

  percent.vec <- rep(NA,length(eps.vec))

  for(it in 1:length(eps.vec)){

    print(it / length(eps.vec))

    coef(m1)['eps'] <- eps.vec[it]

    nsim <- 200

    x <- simulate(m1,nsim=nsim,as.data.frame=TRUE,include.data=F)

    x %>%
      dplyr::group_by(sim) %>%
      dplyr::summarize(p1 = specfun(deaths1),
                       p2 = specfun(deaths2),
                       rp1 = round(p1),
                       rp2 = round(p2),
                       diff = abs(round(p2) - round(p1))
                       ) -> x.sum
    # ) %>%
    # subset(rp1 != 0 ) %>%
    # subset( rp2 !=0)-> x.sum

    percent.vec[it] <-  1 - length(which(x.sum$diff > 0)) / nsim

    x %>%
      subset(sim == 13) %>%
      ggplot(aes(time,deaths1))+geom_line()+geom_line(aes(time,deaths2),col='red')

  }

  diff.df <- data.frame(eps = eps.vec,
                        percent = percent.vec)

  return(diff.df)
}


make.CIs <- function(args){

  save.name <- paste0('Data/InSilico/',args[1],'_',args[2],'_sine_0_0.5_monthly.RData')
  new.save.name <- paste0('Data/InSilico/',args[2],'_',args[1],'_sine_0_0.5_monthly.RData')

  if(file.exists(save.name)){
    load(save.name)
  } else{
    load(new.save.name)
  }

  if(file.exists(save.name) && file.exists(new.save.name)){
    ind <- which.max(c(file.info(save.name)$mtime,file.info(new.save.name)$mtime))
    load(c(save.name,new.save.name)[ind])
  }

  best <- which.max(mf1$loglik)

  max.ll <- mf1$loglik[best]

  CI.ind <- which(max.ll - mf1$loglik < 1.92)
  print(length(CI.ind))
  CIs <- range(mf1$eps[CI.ind])

  lo <- CIs[1]
  hi <- CIs[2]
  eps <- mf1$eps[best]

  return(list(lo=lo,hi=hi,eps=eps,mf1=mf1,max.ll=max.ll))
}

main.region <- 'east'
r1 <- 'central'
r2 <- paste0('OR.',main.region)

args <- c(main.region,r1)

NW.m1 <- make.pomp.object(args=c(main.region,r1))

W.OW.m1 <- make.pomp.object(c(main.region,r2))

NW.mixing <- make.mixing(NW.m1)

W.OW.mixing <- make.mixing(W.OW.m1)

NW.cis <- make.CIs(args=c(main.region,r1))

W.OW.cis <- make.CIs(c(main.region,r2))



NW.sim.plot <- ggplot(NW.mixing,aes(eps,y=100*percent),fill='black')+
  geom_bar(stat='identity',alpha=0.9,col='black')+
  #ylab(paste0('percent out of phase simulations (',args[1],' - ',args[2],')'))+xlab('coupling rate')+
  ylab(paste0('percent synchronized simulations (',main.region,' - ',r1,')'))+
  xlab('coupling rate')+
  #geom_vline(xintercept = eps)+
  geom_vline(xintercept = NW.cis$lo,linetype=2)+
  geom_vline(xintercept = NW.cis$hi,linetype=2)+
  ylim(0,100)


NW.ll.plot <- ggplot(NW.cis$mf1,aes(eps,loglik))+geom_point(col='black',alpha=0.9)+
  geom_hline(yintercept = c(NW.cis$max.ll- 1.92),linetype=2)+geom_vline(xintercept = c(NW.cis$lo,NW.cis$hi),linetype=2)+
  #stat_smooth(method='loess',col='grey50')+
  xlab('coupling rate') + ylab("profile log likelihood")+
  scale_y_continuous(breaks=round(seq(NW.cis$max.ll-20,NW.cis$max.ll-1,length=2)))+
  ggtitle(NULL)#+xlim(0,0.1)



W.OW.sim.plot <- ggplot(W.OW.mixing,aes(eps,y=100*percent),fill='black')+geom_bar(stat='identity',alpha=0.9,col='black')+
  ylab(paste0('percent synchronized simulations (',main.region,' - outer ',main.region,')'))+
  xlab('coupling rate')+
  #geom_vline(xintercept = eps)+
  geom_vline(xintercept = W.OW.cis$lo,linetype=2)+
  geom_vline(xintercept = W.OW.cis$hi,linetype=2)+
  ylim(0,100)


W.OW.ll.plot <- ggplot(W.OW.cis$mf1,aes(eps,loglik))+geom_point(col='black',alpha=0.9)+
  geom_hline(yintercept = c(W.OW.cis$max.ll- 1.92),linetype=2)+geom_vline(xintercept = c(W.OW.cis$lo,W.OW.cis$hi),linetype=2)+
  #stat_smooth(method='loess',col='grey50')+
  xlab('coupling rate') + ylab("profile log likelihood")+
  scale_y_continuous(breaks=round(seq(W.OW.cis$max.ll-3,W.OW.cis$max.ll-1,length=2)))+
  ggtitle(NULL)#+xlim(0,0.1)


require(cowplot)

jpeg(file=paste0(args[1],'_',args[2],'_sim.jpg'),height=8,width=8,units='in',res=300)

first.row <- plot_grid(NW.ll.plot,W.OW.ll.plot,align='h',labels=c('A','B'))
second.row <- plot_grid(NW.sim.plot,W.OW.sim.plot,align='h',labels=c('C','D'))

p1 <-plot_grid(first.row,second.row,labels=NULL,rel_heights = c(1,2),nrow = 2)

print(p1)

dev.off()
