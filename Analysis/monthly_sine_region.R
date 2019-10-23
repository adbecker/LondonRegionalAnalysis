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

data <- subset(data,scale == args)

## ----rprocess------------------------------------------------------------
rproc <- Csnippet("
                  double beta, br, seas, foi, dw, births;
                  double rate[7], trans[7];

                  // cohort effect
                  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
                  br = cohort*birthrate/dt + (1-cohort)*birthrate;
                  else
                  br = (1.0-cohort)*birthrate;

                  seas = beta0*(1 + amp*sin(2*M_PI*t + phase));

                  // transmission rate
                  beta = seas;

                  // force of infection
                  foi = beta*pow(I+iota,alpha)/pop;
                  // white noise (extra-demographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);

                  rate[0] = foi*dw/dt;  //infection rate (stochastic)
                  rate[1] = mu;  		    // natural S death
                  rate[2] = sigma;		  // rate of ending of latent stage
                  rate[3] = mu;			    // natural E death
                  rate[4] = gamma;		  // recovery
                  rate[5] = mu;			    // natural I death
                  rate[6] = phi;

                  // Poisson births
                  births = rpois(br*dt);

                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  reulermultinom(1,C,&rate[6],dt,&trans[6]);

                  //trans[0] = rbinom(S, 1.0 - exp(-rate[0]* dt));
                  //trans[1] = rbinom(S, 1.0 - exp(-rate[1]* dt));
                  //trans[2] = rbinom(E, 1.0 - exp(-rate[2]* dt));
                  //trans[3] = rbinom(E, 1.0 - exp(-rate[3]* dt));
                  //trans[4] = rbinom(I, 1.0 - exp(-rate[4]* dt));
                  //trans[5] = rbinom(I, 1.0 - exp(-rate[5]* dt));
                  //trans[6] = rbinom(C, 1.0 - exp(-rate[6]* dt));

                  S += births - trans[0] - trans[1];
                  E += trans[0] - trans[2] - trans[3];
                  I += trans[2] - trans[4] - trans[5];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[4];           // true incidence
                  //D += trans[6];
                  D = rpois(phi*C);
                  contact = beta;
                  ")
## ----initializer---------------------------------------------------------
initlz <- Csnippet("
                   double m = pop/(S_0+E_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   E = nearbyint(m*E_0);
                   I = nearbyint(m*I_0);
                   W = 0;
                   C = 0;
                   D = 0;
                   ")

## ----dmeasure------------------------------------------------------------
dmeas <- Csnippet("
                  double m = rho*D;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (deaths > 0.0) {
                  lik = pnorm(deaths+0.5,m,sqrt(v)+tol,1,0)-pnorm(deaths-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(deaths+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  if(ISNA(deaths)) {
                  lik = (give_log) ? 0 : 1;
                  }
                  //if (give_log) lik = log(lik);
                  ")

## ----rmeasure------------------------------------------------------------
rmeas <- Csnippet("
                  double m = rho*D;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  deaths = rnorm(m,sqrt(v)+tol);
                  if (deaths > 0.0) {
                  deaths = nearbyint(deaths);
                  } else {
                  deaths = 0.0;
                  }
                  ")

## ----transforms----------------------------------------------------------
toEst <- Csnippet("
                  Tbeta0=log(beta0);
                  Tamp = logit(amp);
                  Tphase=log(phase);
                  Tmu=log(mu);
                  Tphi = logit(phi);
                  Tsigma = log(sigma);
                  Tgamma = log(gamma);
                  Talpha = log(alpha);
                  Tiota = log(iota);
                  Trho = logit(rho);
                  Tcohort = logit(cohort);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  to_log_barycentric(&TS_0, &S_0, 4);
                  ")

fromEst <- Csnippet("
                    Tbeta0=exp(beta0);
                    Tamp = expit(amp);
                    Tphase=exp(phase);
                    Tmu=exp(mu);
                    Tsigma = exp(sigma);
                    Tphi = expit(phi);
                    Tgamma = exp(gamma);
                    Talpha = exp(alpha);
                    Tiota = exp(iota);
                    Trho = expit(rho);
                    Tcohort = expit(cohort);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    from_log_barycentric(&TS_0, &S_0, 4);
                    ")


data %>%
  dplyr::select(time,pop,birthrate) -> covar

data %>%
  dplyr::select(time,deaths) -> infdata

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
       zeronames=c("C","W","D",'contact'),
       statenames=c("S","E","I","C","W",'D','contact'),
       paramnames=c('phi',"sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort",
                    "S_0","E_0","I_0","R_0",'mu',
                    'amp','phase','beta0')
  ) -> m1

## all units in years
## inf period of 5 days
## lat period of 8 days

R0 <- 22
beta <- R0 * 73

# if(args %in% c('west','OR.west','south','OR.south')){
#   amp <- 0.05
#   phase <- 1.5
# }else{
#   amp <- 0.15
#   phase <- 0.5
# }

amp <- 0.1
phase <- 1

params <- c('phi'=0.015,"sigma"=365/8,"gamma"=365/5,"alpha"=1,"iota"= 20,
            "rho"=1,"sigmaSE"=0.05,"psi"=0.1,"cohort"=0,
            "S_0"=0.05,"E_0"=5.14e-05,"I_0"=5.14e-05,"R_0"=0.95,'mu'=1/50,
            'amp' = amp, 'phase' = phase,'beta0'=beta)

m1 <- pomp(m1,params = params)

params['beta0']/params['gamma']

nsim <- 9
x <- simulate(m1,nsim=nsim,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,mapping=aes(x=time,y=deaths,group=sim,color=(sim=="data")))+
  geom_line()+
  scale_color_manual(values=c(`TRUE`="blue",`FALSE`="red"))+
  guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)+
  scale_y_sqrt()+
  theme_bw()+theme(strip.text=element_blank()) -> pl

pl
#require(doSNOW)
require(foreach)
require(doParallel)
#registerDoParallel(cores = 2)

parms.true <- params


if(inSLURM){

  require(doMPI)

  cl <-startMPIcluster()

  nsamples <- 400
  registerDoMPI(cl)

}else{

  require(doParallel)
  cl <- makeCluster(4,'PSOCK')
  registerDoParallel(cl)
  nsamples <- 4

}

est.pars <- c('beta0','amp','phase','cohort',
              'psi','iota','sigmaSE','phi',
              'S_0','E_0','I_0','R_0')

if(!inSLURM){
  stop()
}

#stew(file="local_search.rda",{
w1 <- getDoParWorkers()

mf1<- foreach(i=1:nsamples,
              .packages='pomp',
              .combine=rbind,
              .errorhandling="remove",
              .options.multicore=list(set.seed=TRUE)
) %dorng% {

  tic <- Sys.time()

  guess <- coef(m1)
  guess[est.pars] <- guess[est.pars] * runif(length(est.pars),0.75,1.25)

  guess['amp'] <- runif(1,0,0.4)
  guess['phase'] <- runif(1,0,4)

  #guess <- apply(box,1,function(x)runif(1,x[1],x[2]))
  mf <- mif2(m1,
             #start=c(guess,c(params['alpha'],params['mu'],params['sigma'],params['gamma']),params['rho']),
             start=guess,
             Np=4000,
             Nmif=20,
             cooling.type="hyperbolic",
             cooling.fraction.50=0.5,
             transform=TRUE,
             verbose=TRUE,
             rw.sd=rw.sd(
               beta0=0.02,amp=0.02, phase=0.02,
               phi=0.02, psi=0.02,
               iota=0.02,sigmaSE=0.02,
               S_0=ivp(0.02),E_0=ivp(0.02),I_0=ivp(0.02),R_0=ivp(0.02)
             )
  )

  mf <- continue(mf,
                 Np = 3000,
                 Nmif = 40,
                 cooling.fraction=0.1,
                 verbose=T,
                 rw.sd=rw.sd(
                   #phi = 0.02,psi=0.02,
                   beta0 = 0.01, amp=0.01,phase=0.01
                   #S_0=ivp(0.02),E_0=ivp(0.02),I_0=ivp(0.02),R_0=ivp(0.02)
                 )
  )
  #mf <- continue(mf,Nmif=20,cooling.fraction=0.2)

  pf <- replicate(10, pfilter(mf, Np = 2000))
  ll <- sapply(pf,logLik)
  ll <- logmeanexp(ll, se = TRUE)
  nfail <- sapply(pf,getElement,"nfail")

  toc <- Sys.time()
  etime <- toc-tic
  units(etime) <- "hours"

  data.frame(as.list(coef(mf)),
             loglik = ll[1],
             loglik.se = ll[2],
             nfail.min = min(nfail),
             nfail.max = max(nfail),
             etime = as.numeric(etime))

}
print('finished MIF')

mf1 <- mf1[order(mf1$loglik),]

save(m1,mf1,file=paste0(args,'_monthly_sine.RData'))


