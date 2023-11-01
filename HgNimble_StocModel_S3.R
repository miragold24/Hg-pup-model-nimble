#####################
### NIMBLE Hg Sim
#####################

## load libraries, input data and functions
# library(nimble)
# library(coda)
# library(MCMCvis)
#getwd() 
#source("HgNimble_input_S3stoc.R")
#source("HgSim_S3.R")


stoccode <- nimbleCode({
    #priors
    Nborn ~ dpois(500)
    #Nborn ~ dgamma(shape=10000,rate=20) #(pups simulated as discrete&continuous)
    # Nborn ~ dgamma(shape=50000, rate=100) #still mean of 500, but more diffuse
    #start ~ dunif(1,55) # start is 32 in simulation
    #mu.bday ~ dnorm(45,sd=sd.bday)
    #sd.bday ~ dunif(9,10)
    # #constraints
    #one ~ dconstraint(start < mu.bday)
    

# likelihood

  #realBorn[1:ndays] <- rbinom(n = ndays, size = Nborn, prob = Pb_sum[1:ndays])
  #cumsum_born[1:ndays]<-sum(realBorn[1:ndays])
    
    
    ###
    ### derived quantities for log scale calculations 
    ### 
    lCV <- sd.bday/(mu.bday-start)
    lsd <- sqrt(log(lCV^2+1))
    lmu <- log(mu.bday-start)-(0.5*(lsd^2))
    
    
  # Pb transition
  for(t in 2:ndays){
    Pb[t] <- plnorm((t-start), lmu, lsd) - plnorm((t-start-1), lmu, lsd)
    
  }
  Pb_sum[1:ndays] <- Pb[1:ndays]/sum(Pb[1:ndays])
  

  # for (i in 1:ndays) {
  #   cumsum_born[i] <- sum(realBorn[1:i])
  # }
  # account for probabilities of birth across season and 
  # transitions to moult and leave stages 

      ###
      ###  Pm probabilities
      ###

      for (t in 1:amax){
        Pm[t] <- pnorm(t, mu.mday, sd.mday) - pnorm(t-1, mu.mday, sd.mday)
      }

      Pm_sum[1:amax] <- Pm[1:amax]/sum(Pm[1:amax])
      
      ###
      ###  Pl probabilities 
      ###
      
      for (t in 1:amax){
        Pl[t] <- pnorm(t, mu.lday, sd.lday) - pnorm(t-1, mu.lday, sd.lday)
      }
      
      Pl_sum[1:amax] <- Pl[1:amax]/sum(Pl[1:amax])
    
      ###
      ## stochastic process
      ###  
      
     whitePups[1:ndays] <- NimSimS3_W(realBorn[1:ndays],Nborn,Pb_sum[1:amax],Pm_sum[1:amax],Pl_sum[1:amax])
     moultPups[1:ndays] <- NimSimS3_M(realBorn[1:ndays],Nborn,Pb_sum[1:amax],Pm_sum[1:amax],Pl_sum[1:amax])
    #================================================================
    ### observation model:
    # all white and moulted counted on each survey day
    # account for imperfect detection and imperfect classification error
    # #================================================================
     
     #### FUNCTION TO ADD SUM OF TWO BINOMIALS WHICH IS ADDED WITH PMF(B1)+PMF(B2)
     # ============================================================================
     
     for (obs in 1:length(obsdays)){
       wObs[obs] ~ dbinom(POWhite, whitePups[obsdays[obs]])
       mObs[obs] ~ dbinom(POMoult, moultPups[obsdays[obs]])
       
     }

})


## supply data, constants and inits
## supply data, constants and inits
data_list=list(wObs=sim_S3$whiteObs,mObs=sim_S3$moultObs)

constants_list=list(ndays=ndays,amax=amax,obsdays=obsdays,
                    start=start,
                    mu.bday=mu.bday,
                    sd.bday=sd.bday, #params to be estimated
                    mu.mday=mu.mday,sd.mday=sd.mday,mu.lday=mu.lday,sd.lday=sd.lday,
                    POWhite=POWhite,POMoult=POMoult,PCMoult=PCMoult)#,PCWhite=PCWhite)
#Pwhite=POWhite*PCWhite*(POMoult*(1-PCMoult)),
#Pmoult=POMoult*PCMoult)

inits_list <- list(Nborn=runif(1,500,510), 
                   # mm=rep(0,length(obsdays)),mw=rep(0,length(obsdays)),ww=rep(0,
                   # length(obsdays)),wm=rep(0,length(obsdays)),
                   # whiteCount=rep(0,length(obsdays)),moultCount=rep(0,length(obsdays)),
                   Pb=rep(0,ndays),Pm=rep(0,amax),Pl=rep(0,amax)
                   #start = runif(1,31,33),
                   #mu.bday=runif(1,44,46),
                   #sd.bday= 10,
                   #one=1
                    )
################
model <- nimbleModel(code=stoccode,inits=inits_list,data=data_list,
                     constants=constants_list,
                     dimensions=list(born=ndays,
                                     #wmObs=length(obsdays),
                                     #whitePups=ndays,moultPups=ndays,
                                     # moultCount=length(obsdays),
                                     Pb_sum=ndays,Pm_sum=ndays,Pl_sum=ndays))

nc <- 2   ;   ni <- 11000   ;   nb <- 1000   ;   nt <- 0
# Parameters monitored

params = c('Nborn','whitePups','moultPups'),
          # 'whiteCount','moultCount')
           #'start','mu.bday', 'sd.bday')

# #build MCMC object
HgMCMC <- buildMCMC(model,monitors=params)
runMCMC(HgMCMC,niter=1)# can run this now uncompiled, but will be slow. 
                                #set =1 bc with one calculation inits are used
                                #and something wrong with initial values 

## MCMC configuration: what nodes samplers operate on, any parameters needed, etc
#HgMCMCconfiguration <- configureMCMC(model)
# HgMCMCconfiguration$printSamplers()
# # #compile the model object first for the MCMC to compile 
Cmodel<- compileNimble(model,showCompilerOutput = F)
CHgMCMC <- compileNimble(HgMCMC,project=Cmodel)
# CHgMCMC$run(10000) 
# samples<-as.matrix(CHgMCMC$mvSamples)
# plot(samples[ , 'Nborn'], type = 'l', xlab = 'iteration',  
#      ylab = expression(Nborn))
# plot(samples[ , 'mu.bday'], type = 'l', xlab = 'iteration',  
#      ylab = expression(mu.bday))
# plot(samples[ , 'start'], type = 'l', xlab = 'iteration',  
#      ylab = expression(start))
## continue the same run of the MCMC for 10000 more iterations: 
# CHgMCMC$run(10000,reset=FALSE)
# samples_update<- as.matrix(CHgMCMC$mvSamples)
# plot(samples_update[ , 'Nborn'], type = 'l', xlab = 'iteration',  
#  ylab = expression(Nborn))
# plot(samples_update[ , 'mu.bday'], type = 'l', xlab = 'iteration',  
#      ylab = expression(mu.bday))
# plot(samples_update[ , 'start'], type = 'l', xlab = 'iteration',  
#      ylab = expression(start))

#
# now runMCMC on compiled model object/plot results: 
# 
# mcmc.samples <- runMCMC(CHgMCMC,nchains=2,nburnin=1000,niter=10000,samplesAsCodaMCMC = TRUE,summary=TRUE)
# mcmc.samples<-as.mcmc(runMCMC_samples)


####################
mcmc.out <-nimbleMCMC(model=model,nchains =nc,niter = ni, nburnin = nb,
                     monitors = params, samplesAsCodaMCMC = TRUE,summary=TRUE)
# ,WAIC=TRUE


# #mcmc.out$summary
# # ### check convergence and look at plots 
#  mcmc<-mcmc.out$samples
#  coda::gelman.diag(mcmc[,"Nborn"])
# # # ## look at trace plots
#  plot(mcmc[,"Nborn"])
# # 
# MCMCtrace(object = mcmc,
#           pdf = FALSE,
#           ind = TRUE, # separate density lines per chain
#           Rhat=F,
#           n.eff=F,
#           params = "start")

## troubleshooting 

#model$plot()
model$check()
#model$checkBasics()
model$initializeInfo()

