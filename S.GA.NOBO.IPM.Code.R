###############################################################################
## Code for running an integrated population model estimating population dynamics
##    and demographic rates for a population of northern bobwhites (Colinus
##    virginianus) in southern Georgia, USA, 1998 - 2022.
##  Lewis, W. B., J. A. Rectenwald, D. C. Sisson, and J. A. Martin
##############################################################################

require(parallel)

load("S.GA.NOBO.IPM.data.gzip")

str(NOBO.IPM.data)

# See metadata on GitHub page for description of gzip contents

# Putting model and initial values into function for running in parallel
NOBO_MCMC_code <- function(seed, mod.data, mod.constants, monitors){
  
  require(nimble)
  require(truncnorm)
  
  NOBO.model.code <- nimbleCode({
    
    # Priors
    ##############################################################################
    
    # Initial population sizes for April population
    for(r in 1:4){
      S.n.f[r] ~ dunif(0,1000)
    }
    S.nb[1,1,1] <- round(S.n.f[1])
    S.nb[1,2,1] <- round(S.n.f[2])
    S.nb[2,1,1] <- round(S.n.f[3])
    S.nb[2,2,1] <- round(S.n.f[4])
    
    # Informative prior for covey size and availability, both vary by year
    for (q in 1:nyears){
      covey.size[q] ~ T(dnorm(csize.mean[q], sd=csize.sd[q]), csize.min, csize.max)
      availability.logit[q] ~ T(dnorm(avail.mean[q], sd=avail.sd[q]), avail.min, avail.max)
      availability[q] <- ilogit(availability.logit[q])
    }
    
    # Informative prior for detection of coveys. Doesn't vary by year
    pdet.logit ~ T(dnorm(pdet_mean, sd=pdet_sd), pdet_min, Inf)
    pdet <- ilogit(pdet.logit)
    
    for (d in 1:nyears){ # Looping over years
      
      # Juvenile survival
      for(m in 1:n.breeding.months){
        phi.j.daily.logit[m,d] ~ dnorm(phi.j.prior.mean[m], sd=phi.j.prior.sd)
        phi.j.daily[m,d] <- ilogit(phi.j.daily.logit[m,d])
        phi.j[m,d] <- pow(phi.j.daily[m,d],30)
      }
      J.surv[d,1:n.breeding.months] <- c(pow(phi.j[1,d],4),
                                         pow(phi.j[2,d],3),
                                         pow(phi.j[3,d],2),
                                         phi.j[4,d])
      
      # Breeding survival
      phi.b.logit[d] ~ dnorm(phi.b.hyper.logit, sd=phi.b.hyper.logit.sd) # Yearly values come from distribution around hyperparameter
      phi.b[d] <- ilogit(phi.b.logit[d] + DD.phi.b*N.Apr.Tot[d]) # Density-dependence
      
      for (e in 1:2){ # Looping over sexes (1=male, 2=female)
        
        # Non-breeding survival
        for(f in 1:2){ # Looping over ages (1=adult, 2=subadult)
          phi.nb.logit[f,e,d] ~ dnorm(phi.nb.hyper.logit[f,e], sd=phi.nb.hyper.logit.sd[f,e]) # Yearly values come from distribution around hyperparameter
          phi.nb[f,e,d] <- ilogit(phi.nb.logit[f,e,d] + DD.phi.nb*(N.Oct[1,1,d] + N.Oct[1,2,d] + N.Oct[2,1,d] + N.Oct[2,2,d])) # Density-dependence
        }
        
        # Productivity
        for(g in 1:n.breeding.months){
          productivity.log[g,e,d] ~ dnorm(prod.year.hyper.log[g,e], sd=prod.year.hyper.log.sd[g,e]) # Yearly values come from distribution around hyperparameter
          productivity[g,e,d] <- exp(productivity.log[g,e,d] + DD.prod*(N.b[g+2,1,d]+N.b[g+2,2,d]))
        }
      }
      for(g in 1:n.breeding.months){
        constraint.m.prod[g,d] ~ dconstraint(productivity[g,1,d] <= productivity[g,2,d]) # Male productivity has to be <= female productivity
      }
    }
    
    phi.b.hyper.logit <- logit(Survival.b.hyper)
    Survival.b.hyper ~ dunif(0, 1)
    phi.b.hyper.logit.sd ~ dunif(0, 100)
    # Enforcing negative density dependence
    DD.prod ~ dnorm(0, sd=100)
    constraint.DD.prod ~ dconstraint(DD.prod <= 0)
    DD.phi.nb ~ dnorm(0, sd=100)
    constraint.DD.phi.nb ~ dconstraint(DD.phi.nb <= 0)
    DD.phi.b ~ dnorm(0, sd=100)
    constraint.DD.phi.b ~ dconstraint(DD.phi.b <= 0)
    
    for(o in 1:2){ # Looping over sexes (1=male, 2=female)
      for(h in 1:2){ # Looping over ages (1=adults, 2=subadults)
        phi.nb.hyper.logit[h,o] <- logit(Survival.nb.hyper[h,o])
        Survival.nb.hyper[h,o] ~ dunif(0, 1)
        phi.nb.hyper.logit.sd[h,o] ~ dunif(0, 100)
      }
      
      for(u in 1:n.breeding.months){ # Looping over months of chick production (June-September)
        prod.year.hyper.log[u,o] <- log(Prod.hyper[u,o])
        Prod.hyper[u,o] ~ dgamma(1,1)
        prod.year.hyper.log.sd[u,o] ~ dunif(0, 100)
      }
    }
    
    ##############################################################################
    # State Process
    ##############################################################################
    
    
    
    # First Year
    ##############################################################################
    
    for(k in 1:n.breeding.months){
      Chicks.sex[k,1,1] ~ dbin(0.5,Chicks[k,1,1]+Chicks[k,2,1]) # Male chicks
      Chicks.sex[k,2,1] <- Chicks[k,1,1] + Chicks[k,2,1] - Chicks.sex[k,1,1] # Female chicks
    }
    
    # Ratios of age/sex (Ratio) and montly cohorts (MC). Relating count and trap data
    Ratio.Nov[1,1:3] <- c(N.Nov[1,1,1] / (N.Nov[1,1,1] + N.Nov[1,2,1] + N.Nov[2,1,1] + N.Nov[2,2,1]),
                          N.Nov[1,2,1] / (N.Nov[1,1,1] + N.Nov[1,2,1] + N.Nov[2,1,1] + N.Nov[2,2,1]),
                          (N.Nov[2,1,1] + N.Nov[2,2,1]) / (N.Nov[1,1,1] + N.Nov[1,2,1] + N.Nov[2,1,1] + N.Nov[2,2,1]))
    MC.Nov[1,1:4] <- c((Recruit.Nov[1,1,1] + Recruit.Nov[1,2,1]) / (N.Nov[2,1,1] + N.Nov[2,2,1]),
                       (Recruit.Nov[2,1,1] + Recruit.Nov[2,2,1]) / (N.Nov[2,1,1] + N.Nov[2,2,1]),
                       (Recruit.Nov[3,1,1] + Recruit.Nov[3,2,1]) / (N.Nov[2,1,1] + N.Nov[2,2,1]),
                       (Recruit.Nov[4,1,1] + Recruit.Nov[4,2,1]) / (N.Nov[2,1,1] + N.Nov[2,2,1]))
    # Can't separate June and July cohorts reliably in December
    MC.Dec[1,1:3] <- c( (Recruit.Dec[1,1,1] + Recruit.Dec[1,2,1] + Recruit.Dec[2,1,1] + Recruit.Dec[2,2,1]) / (N.Dec[2,1,1] + N.Dec[2,2,1]),
                        (Recruit.Dec[3,1,1] + Recruit.Dec[3,2,1]) / (N.Dec[2,1,1] + N.Dec[2,2,1]),
                        (Recruit.Dec[4,1,1] + Recruit.Dec[4,2,1]) / (N.Dec[2,1,1] + N.Dec[2,2,1]))
    
    # Total abundance and density
    N.Nov.Tot[1] <- N.Nov[1,1,1] + N.Nov[1,2,1] + N.Nov[2,1,1] + N.Nov[2,2,1]
    D.Nov.Tot[1] <- N.Nov.Tot[1] / Area
    N.Apr.Tot[1] <- N.b[1,1,1] + N.b[1,2,1]
    D.Apr.Tot[1] <- N.Apr.Tot[1] / Area
    
    for(p in 1:2){ # Looping over sexes (1=male, 2=female)
      N.b[1,p,1] <- S.nb[1,p,1] + S.nb[2,p,1] # Adults alive at start of April
      
      # Adults surviving to subsequent months of the breeding season
      for(m in 2:n.b.months){
        N.b[m,p,1] ~ dbin(phi.b[1], N.b[m-1,p,1])
      }
      
      # Productivity model. Assuming 50/50 age ratio of chicks. Assuming juvenile survival differs by month of hatching, but doesn't differ by sex
      for (k in 1:n.breeding.months){
        Chicks[k,p,1] ~ dpois(productivity[k,p,1] * N.b[k+2,p,1]) # Chicks produced
        Recruit.Oct[k,p,1] ~ dbin(J.surv[1,k], Chicks.sex[k,p,1]) # Subadult recruits to start of October
        Recruit.Nov[k,p,1] ~ dbin(phi.nb[2,p,1], Recruit.Oct[k,p,1]) # Subadult recruits to start of November
        Recruit.Dec[k,p,1] ~ dbin(phi.nb[2,p,1], Recruit.Nov[k,p,1]) # Subadult recruits to start of December
      }
      
      # October abundance
      N.Oct[1,p,1] ~ dbin(phi.b[1], N.b[n.b.months,p,1])
      N.Oct[2,p,1] <- Recruit.Oct[1,p,1] + Recruit.Oct[2,p,1] + Recruit.Oct[3,p,1] + Recruit.Oct[4,p,1]
      
      # November abundance
      N.Nov[1,p,1] ~ dbin(phi.nb[1,p,1], N.Oct[1,p,1])
      N.Nov[2,p,1] <- Recruit.Nov[1,p,1] + Recruit.Nov[2,p,1] + Recruit.Nov[3,p,1] + Recruit.Nov[4,p,1]
      
      # December abundance
      N.Dec[1,p,1] ~ dbin(phi.nb[1,p,1], N.Nov[1,p,1])
      N.Dec[2,p,1] <- Recruit.Dec[1,p,1] + Recruit.Dec[2,p,1] + Recruit.Dec[3,p,1] + Recruit.Dec[4,p,1]
    }
    
    
    
    # Subsequent Years
    ##############################################################################
    for(t in 2:nyears){
      
      for(k in 1:n.breeding.months){
        Chicks.sex[k,1,t] ~ dbin(0.5,Chicks[k,1,t]+Chicks[k,2,t]) # Male chicks
        Chicks.sex[k,2,t] <- Chicks[k,1,t] + Chicks[k,2,t] - Chicks.sex[k,1,t] # Female chicks
      }
      
      # Ratios of age/sex (Ratio) and montly cohorts (MC). Relating count and trap data
      Ratio.Nov[t,1:3] <- c(N.Nov[1,1,t] / (N.Nov[1,1,t] + N.Nov[1,2,t] + N.Nov[2,1,t] + N.Nov[2,2,t]),
                            N.Nov[1,2,t] / (N.Nov[1,1,t] + N.Nov[1,2,t] + N.Nov[2,1,t] + N.Nov[2,2,t]),
                            (N.Nov[2,1,t] + N.Nov[2,2,t]) / (N.Nov[1,1,t] + N.Nov[1,2,t] + N.Nov[2,1,t] + N.Nov[2,2,t]))
      MC.Nov[t,1:4] <- c((Recruit.Nov[1,1,t] + Recruit.Nov[1,2,t]) / (N.Nov[2,1,t] + N.Nov[2,2,t]),
                         (Recruit.Nov[2,1,t] + Recruit.Nov[2,2,t]) / (N.Nov[2,1,t] + N.Nov[2,2,t]),
                         (Recruit.Nov[3,1,t] + Recruit.Nov[3,2,t]) / (N.Nov[2,1,t] + N.Nov[2,2,t]),
                         (Recruit.Nov[4,1,t] + Recruit.Nov[4,2,t]) / (N.Nov[2,1,t] + N.Nov[2,2,t]))
      # Can't separate June and July cohorts reliably in December
      MC.Dec[t,1:3] <- c( (Recruit.Dec[1,1,t] + Recruit.Dec[1,2,t] + Recruit.Dec[2,1,t] + Recruit.Dec[2,2,t]) / (N.Dec[2,1,t] + N.Dec[2,2,t]),
                          (Recruit.Dec[3,1,t] + Recruit.Dec[3,2,t]) / (N.Dec[2,1,t] + N.Dec[2,2,t]),
                          (Recruit.Dec[4,1,t] + Recruit.Dec[4,2,t]) / (N.Dec[2,1,t] + N.Dec[2,2,t]))
      
      # Total abundance, density, and finite population growth rate
      N.Nov.Tot[t] <- N.Nov[1,1,t] + N.Nov[1,2,t] + N.Nov[2,1,t] + N.Nov[2,2,t]
      D.Nov.Tot[t] <- N.Nov.Tot[t] / Area
      lambda.Nov[t-1] <- N.Nov.Tot[t] / N.Nov.Tot[t-1]
      N.Apr.Tot[t] <- N.b[1,1,t] + N.b[1,2,t]
      D.Apr.Tot[t] <- N.Apr.Tot[t] / Area
      lambda.Apr[t-1] <- N.Apr.Tot[t] / N.Apr.Tot[t-1]
      
      for(p in 1:2){ # Looping over sexes (1=male, 2=female)
        
        # Number of birds surviving from December to April
        for(d in 1:2){ # Looping over ages (1=adults, 2=subadults)
          S.nb[d,p,t] ~ dbin(pow(phi.nb[d,p,t-1],4), N.Dec[d,p,t-1])
        }
        
        N.b[1,p,t] <- S.nb[1,p,t] + S.nb[2,p,t] # Adults alive at start of April
        
        # Adults surviving to subsequent months of the breeding season
        for(m in 2:n.b.months){
          N.b[m,p,t] ~ dbin(phi.b[t], N.b[m-1,p,t])
        }
        
        # Productivity model. Assuming 50/50 age ratio of chicks. Assuming juvenile survival differs by month of hatching, but doesn't differ by sex
        for (k in 1:n.breeding.months){
          Chicks[k,p,t] ~ dpois(productivity[k,p,t] * N.b[k+2,p,t]) # Chicks produced
          Recruit.Oct[k,p,t] ~ dbin(J.surv[t,k], Chicks.sex[k,p,t]) # Subadult recruits to start of October
          Recruit.Nov[k,p,t] ~ dbin(phi.nb[2,p,t], Recruit.Oct[k,p,t]) # Subadult recruits to start of November
          Recruit.Dec[k,p,t] ~ dbin(phi.nb[2,p,t], Recruit.Nov[k,p,t]) # Subadult recruits to start of December
        }
        
        # October abundance
        N.Oct[1,p,t] ~ dbin(phi.b[t], N.b[n.b.months,p,t])
        N.Oct[2,p,t] <- Recruit.Oct[1,p,t] + Recruit.Oct[2,p,t] + Recruit.Oct[3,p,t] + Recruit.Oct[4,p,t]
        
        # November abundance
        N.Nov[1,p,t] ~ dbin(phi.nb[1,p,t], N.Oct[1,p,t])
        N.Nov[2,p,t] <- Recruit.Nov[1,p,t] + Recruit.Nov[2,p,t] + Recruit.Nov[3,p,t] + Recruit.Nov[4,p,t]
        
        # December abundance
        N.Dec[1,p,t] ~ dbin(phi.nb[1,p,t], N.Nov[1,p,t])
        N.Dec[2,p,t] <- Recruit.Dec[1,p,t] + Recruit.Dec[2,p,t] + Recruit.Dec[3,p,t] + Recruit.Dec[4,p,t]
      }
    }
    
    
    
    ##############################################################################
    # Likelihoods
    ##############################################################################
    
    
    # Trapping 
    for (r in 1:nyears.trap){
      trap.am.af.sub[r,1:3] ~ dmulti(Ratio.Nov[r,1:3], trap.N[r]) # November ratio of adult males, adult females, and subadult recruits
      trap.sub.cohort[r,1:4] ~ dmulti(MC.Nov[r,1:4], trap.cohort.N[r]) # November ratio of recruiting subadults from monthly cohorts
    }
    
    # November harvest
    for(v in 1:harv.ratio.years.n){
      harv.cohort.Nov[v,1:4] ~ dmulti(MC.Nov[harv.ratio.years[v],1:4], harv.cohort.Nov.N[v]) # November ratio of subadults from monthly cohorts
    }
    
    # December harvest
    for(u in 1:harv.ratio.years.Dec.n){
      harv.cohort.Dec[u,1:3] ~ dmulti(MC.Dec[harv.ratio.years.Dec[u],1:3], harv.cohort.Dec.N[u]) # December ratio of subadults from monthly cohorts
    }
    
    # Productivity
    for (s in 1:nyears){
      for (u in 1:n.breeding.months){
        for(y in 1:2){
          chicks.nest.mean[u,y,s] ~ T(dnorm(chicks.nest.act[u,y,s], sd=chicks.nest.sd[u,y,s]), 0, Inf) # Incorporating uncertainty in chick production due to missing data
          chicks.nest.act[u,y,s] ~ dpois(productivity[u,y,s] * N.tracked[u,y,s])
        }
      }
    }
    
    # Fall covey counts
    for (a in 1:nyears){
      count[a] ~ dbinom(availability[a] * pdet, coveys[a]) # Number of coveys detected is a function of number of coveys, availability, and detection probability
      coveys[a] ~ dpois(Avail[a] / covey.size[a]) # Number of coveys is a function of number available divided by average covey size
      Avail[a] ~ dbinom(effort[a],N.Nov.Tot[a]) # The number of bobwhite available to be sampled is conditional on total population size and degree of area surveyed
    }
    
    # Telemetry
    # Placing survival estimates in 3-D array, converting from monthly to biweekly survival
    CH.phi[1,1:nyears,1] <- pow(phi.nb[2,1,1:nyears],0.5)
    CH.phi[1,1:nyears,2] <- pow(phi.nb[2,2,1:nyears],0.5)
    CH.phi[2,1:nyears,1] <- pow(phi.nb[1,1,1:nyears],0.5)
    CH.phi[2,1:nyears,2] <- pow(phi.nb[1,2,1:nyears],0.5)
    CH.phi[3,1:nyears,1] <- pow(phi.b[1:nyears],0.5)
    CH.phi[3,1:nyears,2] <- pow(phi.b[1:nyears],0.5)
    
    for (b in 1:nCH){
      for (c in 1:trackperiods[b]){
        CH.state[b,trackdates[b,c]] ~ dbern(CH.phi[CH.age[b,trackdates[b,c]-1], CH.year[trackdates[b,c]-1], CH.sex[b]])
      }
    }
  })
  
  
  # Setting initial values
  # Seed ensures that initial values are acceptable: 
  #   Initial November population size 1 - 5000
  #   Initial April population size <= 4000
  #   Initial capture and harvest ratios never reach 0
  #   Initial covey number never less than observed covey count
  set.seed(1539015)
  N.Nov.Tot.init <- N.Apr.Tot.init <- Avail.init <- coveys.init <- rep(NA, times=mod.constants$nyears)
  phi.b.init=runif(mod.constants$nyears,0.8,1)
  phi.nb.init=array(runif(mod.constants$nyears*2*2,0.8,1), dim=c(2,2,mod.constants$nyears))
  J.surv.init <- matrix(NA, nrow=mod.constants$nyears, ncol=4)
  availability.init <- rep(NA, times=mod.constants$nyears)
  for(i in 1:mod.constants$nyears){
    J.surv.init[i,1] <- (plogis(runif(1,mod.constants$phi.j.prior.mean[1]-2*mod.constants$phi.j.prior.sd,mod.constants$phi.j.prior.mean[1]+2*mod.constants$phi.j.prior.sd))^30)^4
    J.surv.init[i,2] <- (plogis(runif(1,mod.constants$phi.j.prior.mean[2]-2*mod.constants$phi.j.prior.sd,mod.constants$phi.j.prior.mean[2]+2*mod.constants$phi.j.prior.sd))^30)^3
    J.surv.init[i,3] <- (plogis(runif(1,mod.constants$phi.j.prior.mean[3]-2*mod.constants$phi.j.prior.sd,mod.constants$phi.j.prior.mean[3]+2*mod.constants$phi.j.prior.sd))^30)^2
    J.surv.init[i,4] <- plogis(runif(1,mod.constants$phi.j.prior.mean[4]-2*mod.constants$phi.j.prior.sd,mod.constants$phi.j.prior.mean[4]+2*mod.constants$phi.j.prior.sd))^30
    availability.init[i] <- plogis(rtruncnorm(1,mod.constants$avail.min,mod.constants$avail.max,mod.constants$avail.mean[i],mod.constants$avail.sd[i]))
  }
  productivity.init <- array(NA, dim=c(mod.constants$n.breeding.months,2,mod.constants$nyears))
  productivity.init[,2,] <- runif(mod.constants$nyears*mod.constants$n.breeding.months, 0, mod.constants$n.breeding.months)
  for(i in 1:mod.constants$n.breeding.months){
    for(j in 1:mod.constants$nyears){
      productivity.init[i,1,j] <- runif(1,0,productivity.init[i,2,j])
    }
  }
  covey.size.init <- runif(mod.constants$nyears,mod.constants$csize.min,mod.constants$csize.max)
  S.n.f.init <- runif(4,200,600)
  S.nb.init <- N.Oct.init <- N.Nov.init <- N.Dec.init <- array(NA, dim=c(2, 2, mod.constants$nyears))
  S.nb.init[1,1,1] <- round(S.n.f.init[1],0)
  S.nb.init[1,2,1] <- round(S.n.f.init[2],0)
  S.nb.init[2,1,1] <- round(S.n.f.init[3],0)
  S.nb.init[2,2,1] <- round(S.n.f.init[4],0)
  N.b.init <- array(0, dim=c(mod.constants$n.b.months, 2, mod.constants$nyears))
  Chicks.init <- Chicks.sex.init <- Recruit.Oct.init <- Recruit.Nov.init <- Recruit.Dec.init <- array(0, dim=c(mod.constants$n.breeding.months, 2, mod.constants$nyears))
  Ratio.Nov.init <- matrix(NA, nrow=mod.constants$nyears, ncol=3)
  MC.Nov.init <- matrix(NA, nrow=mod.constants$nyears, ncol=4)
  MC.Dec.init <- matrix(NA, nrow=mod.constants$nyears, ncol=3)
  for(p in 1:2){
    N.b.init[1,p,1] <- S.nb.init[1,p,1] + S.nb.init[2,p,1]
    for(m in 2:nrow(N.b.init)){
      N.b.init[m,p,1] <- rbinom(1,N.b.init[m-1,p,1],phi.b.init[1])
    }
    for(k in 1:nrow(Chicks.init)){
      Chicks.init[k,p,1] <- rpois(1,productivity.init[k,p,1] * N.b.init[k+2,p,1])
    }
  }
  for(k in 1:nrow(Chicks.sex.init)){
    Chicks.sex.init[k,1,1] <- rbinom(1,Chicks.init[k,1,1]+Chicks.init[k,2,1],0.5)
    Chicks.sex.init[k,2,1] <- Chicks.init[k,1,1]+Chicks.init[k,2,1]-Chicks.sex.init[k,1,1]
  }
  for(p in 1:2){
    for(k in 1:nrow(Recruit.Oct.init)){
      Recruit.Oct.init[k,p,1] <- rbinom(1,Chicks.sex.init[k,p,1],J.surv.init[1,k])
      Recruit.Nov.init[k,p,1] <- rbinom(1,Recruit.Oct.init[k,p,1],phi.nb.init[2,p,1])
      Recruit.Dec.init[k,p,1] <- rbinom(1,Recruit.Nov.init[k,p,1],phi.nb.init[2,p,1])
    }
    N.Oct.init[1,p,1] <- rbinom(1,N.b.init[nrow(N.b.init),p,1],phi.b.init[1])
    N.Oct.init[2,p,1] <- sum(Recruit.Oct.init[,p,1])
    N.Nov.init[1,p,1] <- rbinom(1,N.Oct.init[1,p,1],phi.nb.init[1,p,1])
    N.Nov.init[2,p,1] <- sum(Recruit.Nov.init[,p,1])
    N.Dec.init[1,p,1] <- rbinom(1,N.Nov.init[1,p,1],phi.nb.init[1,p,1])
    N.Dec.init[2,p,1] <- sum(Recruit.Dec.init[,p,1])
  }
  Ratio.Nov.init[1,1:3] <- c(N.Nov.init[1,1,1]/sum(N.Nov.init[,,1]),
                             N.Nov.init[1,2,1]/sum(N.Nov.init[,,1]),
                             sum(N.Nov.init[2,,1])/sum(N.Nov.init[,,1]))
  MC.Nov.init[1,1:4] <- c(sum(Recruit.Nov.init[1,,1])/sum(N.Nov.init[2,,1]),
                          sum(Recruit.Nov.init[2,,1])/sum(N.Nov.init[2,,1]),
                          sum(Recruit.Nov.init[3,,1])/sum(N.Nov.init[2,,1]),
                          sum(Recruit.Nov.init[4,,1])/sum(N.Nov.init[2,,1]))
  MC.Dec.init[1,1:3] <- c((sum(Recruit.Dec.init[1,,1])+sum(Recruit.Dec.init[2,,1]))/sum(N.Dec.init[2,,1]),
                          sum(Recruit.Dec.init[3,,1])/sum(N.Dec.init[2,,1]),
                          sum(Recruit.Dec.init[4,,1])/sum(N.Dec.init[2,,1]))
  N.Nov.Tot.init[1] <- sum(N.Nov.init[,,1])
  N.Apr.Tot.init[1] <- sum(N.b.init[1,,1])
  Avail.init[1] <- rbinom(1, N.Nov.Tot.init[1], mod.constants$effort[1])
  coveys.init[1] <- rpois(1, Avail.init[1] / covey.size.init[1])
  for(t in 2:mod.constants$nyears){
    for(p in 1:2){
      for(d in 1:2){
        S.nb.init[d,p,t] <- rbinom(1,N.Dec.init[d,p,t-1],pow(phi.nb.init[d,p,t-1],4))
      }
      N.b.init[1,p,t] <- S.nb.init[1,p,t] + S.nb.init[2,p,t]
      for(m in 2:nrow(N.b.init)){
        N.b.init[m,p,t] <- rbinom(1,N.b.init[m-1,p,t],phi.b.init[t])
      }
      for(k in 1:nrow(Chicks.init)){
        Chicks.init[k,p,t] <- rpois(1,productivity.init[k,p,t] * N.b.init[k+2,p,t])
      }
    }
    for(k in 1:nrow(Chicks.sex.init)){
      Chicks.sex.init[k,1,t] <- rbinom(1,Chicks.init[k,1,t]+Chicks.init[k,2,t],0.5)
      Chicks.sex.init[k,2,t] <- Chicks.init[k,1,t]+Chicks.init[k,2,t]-Chicks.sex.init[k,1,t]
    }
    for(p in 1:2){
      for(k in 1:nrow(Recruit.Oct.init)){
        Recruit.Oct.init[k,p,t] <- rbinom(1,Chicks.sex.init[k,p,t],J.surv.init[t,k])
        Recruit.Nov.init[k,p,t] <- rbinom(1,Recruit.Oct.init[k,p,t],phi.nb.init[2,p,t])
        Recruit.Dec.init[k,p,t] <- rbinom(1,Recruit.Nov.init[k,p,t],phi.nb.init[2,p,t])
      }
      N.Oct.init[1,p,t] <- rbinom(1,N.b.init[nrow(N.b.init),p,t],phi.b.init[t])
      N.Oct.init[2,p,t] <- sum(Recruit.Oct.init[,p,t])
      N.Nov.init[1,p,t] <- rbinom(1,N.Oct.init[1,p,t],phi.nb.init[1,p,t])
      N.Nov.init[2,p,t] <- sum(Recruit.Nov.init[,p,t])
      N.Dec.init[1,p,t] <- rbinom(1,N.Nov.init[1,p,t],phi.nb.init[1,p,t])
      N.Dec.init[2,p,t] <- sum(Recruit.Dec.init[,p,t])
    }
    Ratio.Nov.init[t,1:3] <- c(N.Nov.init[1,1,t]/sum(N.Nov.init[,,t]),
                               N.Nov.init[1,2,t]/sum(N.Nov.init[,,t]),
                               sum(N.Nov.init[2,,t])/sum(N.Nov.init[,,t]))
    MC.Nov.init[t,1:4] <- c(sum(Recruit.Nov.init[1,,t])/sum(N.Nov.init[2,,t]),
                            sum(Recruit.Nov.init[2,,t])/sum(N.Nov.init[2,,t]),
                            sum(Recruit.Nov.init[3,,t])/sum(N.Nov.init[2,,t]),
                            sum(Recruit.Nov.init[4,,t])/sum(N.Nov.init[2,,t]))
    MC.Dec.init[t,1:3] <- c((sum(Recruit.Dec.init[1,,t])+sum(Recruit.Dec.init[2,,t]))/sum(N.Dec.init[2,,t]),
                            sum(Recruit.Dec.init[3,,t])/sum(N.Dec.init[2,,t]),
                            sum(Recruit.Dec.init[4,,t])/sum(N.Dec.init[2,,t]))
    N.Nov.Tot.init[t] <- sum(N.Nov.init[,,t])
    N.Apr.Tot.init[t] <- sum(N.b.init[1,,t])
    Avail.init[t] <- rbinom(1, N.Nov.Tot.init[t], mod.constants$effort[t])
    coveys.init[t] <- rpois(1, Avail.init[t] / covey.size.init[t])
  }
  inits.fun <- function() list(S.n.f=S.n.f.init,
                               covey.size=runif(mod.constants$nyears,mod.constants$csize.min,mod.constants$csize.max),
                               availability.logit=qlogis(availability.init),
                               pdet.logit=rtruncnorm(1,mod.constants$pdet_min,10,mod.constants$pdet_mean,mod.constants$pdet_sd),
                               phi.j.daily.logit=matrix(qlogis(runif(mod.constants$nyears*mod.constants$n.breeding.months,0.5,1)),ncol=mod.constants$nyears, byrow=T),
                               phi.b.logit=qlogis(runif(mod.constants$nyears,0.5,1)),
                               phi.nb.logit=array(qlogis(runif(mod.constants$nyears*2*2,0.5,1)), dim=c(2,2,mod.constants$nyears)),
                               productivity.log=array(runif(mod.constants$nyears*2*mod.constants$n.breeding.months,0,1) + rep(rep(c(-1,1),each=mod.constants$n.breeding.months),times=mod.constants$nyears),dim=c(mod.constants$n.breeding.months,2,mod.constants$nyears)), # Ensures lower for males
                               Survival.b.hyper=runif(1,0,1),
                               phi.b.hyper.logit.sd=runif(1,0,1),
                               DD.prod=runif(1,-0.0001,0),
                               DD.phi.nb=runif(1,-0.0001,0),
                               DD.phi.b=runif(1,-0.0001,0),
                               Survival.nb.hyper=matrix(runif(4,0,1),nrow=2),
                               phi.nb.hyper.logit.sd=matrix(runif(4,0,1),nrow=2),
                               Prod.hyper=matrix(runif(2*mod.constants$n.breeding.months,0.00001,3),ncol=2),
                               prod.year.hyper.log.sd=matrix(runif(2*mod.constants$n.breeding.months,0,1),ncol=2),
                               Chicks.sex=Chicks.sex.init,
                               N.b=N.b.init,
                               Chicks=Chicks.init,
                               Recruit.Oct=Recruit.Oct.init,
                               Recruit.Nov=Recruit.Nov.init,
                               Recruit.Dec=Recruit.Dec.init,
                               N.Oct=N.Oct.init,
                               N.Nov=N.Nov.init,
                               N.Dec=N.Dec.init,
                               S.nb=S.nb.init,
                               chicks.nest.act=ceiling(mod.data$chicks.nest.mean),
                               coveys=coveys.init,
                               Avail=Avail.init)
  
  NOBO.model <- nimbleModel(code = NOBO.model.code,
                         data = mod.data,
                         constants = mod.constants,
                         inits = inits.fun())
  NOBO.mcmc.out <- nimbleMCMC(model = NOBO.model,
                         niter = 120000, nchains = 1, nburnin = 30000, thin=25,
                         samplesAsCodaMCMC=TRUE, monitor=monitors, setSeed = seed)
  return(NOBO.mcmc.out)
}


mod.data <- list(trap.am.af.sub=NOBO.IPM.data$trap.am.af.sub,
                 trap.sub.cohort=NOBO.IPM.data$trap.sub.cohort,
                 harv.cohort.Nov=NOBO.IPM.data$harv.cohort.Nov,
                 harv.cohort.Dec=NOBO.IPM.data$harv.cohort.Dec,
                 chicks.nest.mean=NOBO.IPM.data$chicks.nest.mean,
                 count=NOBO.IPM.data$count,
                 CH.state=NOBO.IPM.data$CH.state,
                 constraint.m.prod=matrix(1, nrow=4, ncol=NOBO.IPM.data$nyears),
                 constraint.DD.prod=1,
                 constraint.DD.phi.b=1,
                 constraint.DD.phi.nb=1)

mod.constants <- list(N.tracked=NOBO.IPM.data$N.tracked,
                      effort=NOBO.IPM.data$effort,
                      trap.N=NOBO.IPM.data$trap.N,
                      trap.cohort.N=NOBO.IPM.data$trap.cohort.N,
                      harv.cohort.Nov.N=NOBO.IPM.data$harv.cohort.Nov.N,
                      harv.cohort.Dec.N=NOBO.IPM.data$harv.cohort.Dec.N,
                      harv.ratio.years=NOBO.IPM.data$harv.ratio.years,
                      harv.ratio.years.n=NOBO.IPM.data$harv.ratio.years.n,
                      harv.ratio.years.Dec=NOBO.IPM.data$harv.ratio.years.Dec,
                      harv.ratio.years.Dec.n=NOBO.IPM.data$harv.ratio.years.Dec.n,
                      CH.age=NOBO.IPM.data$CH.age,
                      CH.year=NOBO.IPM.data$CH.year,
                      CH.sex=NOBO.IPM.data$CH.sex,
                      nyears=NOBO.IPM.data$nyears,
                      nyears.trap=NOBO.IPM.data$nyears.trap,
                      n.b.months=NOBO.IPM.data$n.b.months, # April-September
                      n.breeding.months=NOBO.IPM.data$n.breeding.months, # June-September
                      nCH=NOBO.IPM.data$nCH,
                      trackdates=NOBO.IPM.data$trackdates,
                      trackperiods=NOBO.IPM.data$trackperiods,
                      csize.mean=NOBO.IPM.data$csize.mean,
                      csize.sd=NOBO.IPM.data$csize.sd,
                      csize.min=NOBO.IPM.data$csize.min,
                      csize.max=NOBO.IPM.data$csize.max,
                      Area=NOBO.IPM.data$Area,
                      phi.j.prior.mean=NOBO.IPM.data$phi.j.prior.mean,
                      phi.j.prior.sd=NOBO.IPM.data$phi.j.prior.sd,
                      chicks.nest.sd=NOBO.IPM.data$chicks.nest.sd,
                      avail.mean=NOBO.IPM.data$avail.mean,
                      avail.sd=NOBO.IPM.data$avail.sd,
                      avail.min=NOBO.IPM.data$avail.min,
                      avail.max=NOBO.IPM.data$avail.max,
                      pdet_mean=NOBO.IPM.data$pdet_mean,
                      pdet_sd=NOBO.IPM.data$pdet_sd,
                      pdet_min=NOBO.IPM.data$pdet_min)

monitors <- c("S.nb","covey.size","availability","pdet","phi.j.daily","phi.b.logit","DD.phi.b",
              "phi.b","phi.nb.logit","DD.phi.nb","phi.nb","productivity.log","DD.prod",
              "productivity","Survival.b.hyper","phi.b.hyper.logit.sd","Survival.nb.hyper",
              "phi.nb.hyper.logit.sd","Prod.hyper","prod.year.hyper.log.sd","Chicks.sex",
              "Ratio.Nov","MC.Nov","MC.Dec","N.Nov.Tot","D.Nov.Tot","N.Apr.Tot","D.Apr.Tot",
              "N.b","Chicks","Recruit.Oct","Recruit.Nov","Recruit.Dec","N.Oct","N.Nov",
              "N.Dec","lambda.Nov","lambda.Apr","chicks.nest.act","coveys","Avail")

this_cluster <- makeCluster(3) # Creating cluster

NOBO_IPM <- parLapply(cl = this_cluster, X = 1:3, 
                      fun = NOBO_MCMC_code, 
                      mod.data = mod.data,
                      mod.constants = mod.constants,
                      monitors = monitors)

stopCluster(this_cluster)
