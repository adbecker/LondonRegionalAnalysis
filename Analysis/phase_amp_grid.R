
require(deSolve)
require(animation)

phase_amp_grid <- function(stochastic,grid.size,R0,
                           phase.min,phase.max,amp.min,amp.max){



  S_0 <- 0.05
  I_0 <- 5e-5

  gamma <- 73
  sigma <- 45.6
  N1 <- 1000000
  N2 <- 1*N1

  inSLURM <- c(Sys.getenv('SLURM_JOB_ID') != "")

  seir_step <- function(x, params) {

    S <- E <- I <- R <- N <-  rep(NA,length(times_vector))

    S[1] = x[1]
    E[1] = x[2]
    I[1] = x[3]
    R[1] = x[4]

    N <- S[1]+E[1]+I[1]+R[1]

    v =params['v']
    phase = params['phase']
    mu = params["mu"]
    beta1 = params["beta1"]
    amp = params["amp"]
    sigma = params["sigma"]
    gamma = params["gamma"]

    seas <- beta1 * (1 + amp * cos(2 * pi * times_vector - phase))
    seas <- beta1 * (1 + amp * sin(2 * pi * times_vector + phase))

    for(t in 2:length(S)){
      if(stochastic == T){
        dw <- runif(1,0.9,1.1)
      } else{
        dw <- 1
      }
      #dw <- 1

      births <- v*N[t-1]
      foi <- seas[t]*S[t-1]*I[t-1] / N[t-1] * dw
      inf <- sigma*E[t-1]
      rec <- gamma*I[t-1]

      S[t] = S[t-1] + delta.t*(births - foi - mu*S[t-1])
      E[t] = E[t-1] + delta.t*(foi - inf - mu*E[t-1])
      I[t] = I[t-1] + delta.t*(inf - rec - mu*I[t-1])
      R[t] = R[t-1] + delta.t*(rec -  mu*R[t-1])


      N[t] <- S[t] + E[t] + I[t] + R[t]

    }

    return(data.frame(cbind('time'=times_vector,S,E,I,R)))
  }

  '%!in%' <- function(x,y)!('%in%'(x,y))


  time.step <- 1/52/5
  delta.t <- time.step
  duration <- 100

  times_vector <- seq(from=0, to=duration, by=time.step)

  specfun <- function(ts){

    spec <- spectrum(ts,plot=F)
    options <- 1/spec$freq * time.step
    power <- spec$spec

    #plot(options,power,type='l',xlim=c(0,5))

    bi.ind.ll <-which.min(abs(options - 1.75))
    bi.ind.ul <-which.min(abs(options - 2.25))

    an.ind.ll <-which.min(abs(options - 0.75))
    an.ind.ul <-which.min(abs(options - 1.25))

    yr2 <- power[bi.ind.ll:bi.ind.ul]
    yr1 <- power[an.ind.ll:an.ind.ul]

    pr <- sum(yr2)/sum(yr1)

    freq <- options[which.max(power)]
    return(c(freq,pr))
  }


  if(inSLURM){

    require(doMPI)

    cl <-startMPIcluster()

    registerDoMPI(cl)

  }else{

    require(doParallel)
    cl <- makeCluster(4,'PSOCK')
    registerDoParallel(cl)

  }

  phase.options <- seq(phase.min,phase.max,length=grid.size)
  amp.options <- seq(amp.min,amp.max,length=grid.size)

  start <- expand.grid(phase = phase.options, amp = amp.options)

  nsamples <- nrow(start)
  #stew(file="local_search.rda",{
  w1 <- getDoParWorkers()
  d <- foreach(it=1:nsamples,
               .packages=c('deSolve','ggplot2'),
               .errorhandling = 'remove',
               .combine=rbind,
               .options.multicore=list(set.seed=TRUE)
  ) %dorng% {

    phase <- start[it,1]
    amp <- start[it,2]

    beta1 <- R0 * ( gamma)

    beta2 <- beta1
    eps <- 0
    beta.mix <- beta1 * eps

    v <- 1/50
    mu <- 1/50

    parms_vector <- c(mu=mu, v=v, beta1=beta1, beta2=beta2, beta.mix=beta.mix, gamma=gamma, sigma=sigma, amp=amp,phase=phase )

    parms_vector['beta1'] <- beta1
    parms_vector['amp'] <- amp

    S10 <- S_0*N1
    I10 <- I_0*N1
    E10 <- unname(parms_vector['gamma']*I10 / parms_vector['sigma'])
    R10 <- N1 - S10 - I10 - E10

    S20 <- S_0*N2
    I20 <- I_0*N2
    E20 <- unname(parms_vector['gamma']*I20 / parms_vector['sigma'])
    R20 <- N2 - S20 - I20 - E20



    SIR.output <- seir_step(x=c(S=S10, E=E10, I=I10, R=R10),
                            params=parms_vector)

    I1 <- SIR.output$I
    I2 <- SIR.output$I2

    year <- NA

    results <- c(specfun(I1),year)

    res <- data.frame(rbind(c(phase,amp,as.numeric(results))))
    names(res) <-   c('phase',"amp",'freq','pr','Year')
    res
  }

  d$roundfreq <- round(d$freq)
  d$roundfreq[which(d$roundfreq > 3)] <- '3+ '

  d$roundfreq <- as.factor((d$roundfreq))

  d$year.copy <- NA

  file.name <- paste0('IC_stoch_',R0,'.RData')
  #save(d,file=file.name)

  return(d)
}
