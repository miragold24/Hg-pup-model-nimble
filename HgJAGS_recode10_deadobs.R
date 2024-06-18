# JAGS script - with logistic moult curve, logistic survival, and dead observations


# Load libraries -------------------------
library(jagsUI)
library(coda)


# Load data -----------------------

#source('HgSim_kts.R')  
# source('ndx2.R')

# JAGS code -------------------------

jagsfile <- "Hgpup.txt" 

cat("

model{

#=================================
# set the probability of being white as a function of age
#===================================

pmoult.int <- -12   # ~ dnorm(5,1)               # probability of being moult at birth
pmoult.slope <- 0.6 # ~ dunif(-0.5,0)          # daily decrease in the probability of being moult (increase in probability of moulting)
for(a in 1:ndays){   # loop past the last possible age
  logit(pmoult[a]) <- pmoult.int + pmoult.slope*(a-1)    # deterministic for now
}

#=================================
# set the probability of survival as a logit linear function of age
#===================================

  # surv.int.real ~ dbeta(175,1)                           # survival rate at birth
surv.int <-  4    # log(surv.int.real/(1-surv.int.real))       # logit scale intercept for survival
surv.slope <- 0.1 # ~ dunif(0,.5)                               # daily increase in survival (logit scale) (rate of change in surv) 
for(a in 1:ndays){  # loop past the last possible age
  logit(surv[a]) <- surv.int + surv.slope*(a-1)    # deterministic for now
}

#=================================
# set the daily probability of carcass survival as a logit linear function of days since death
#===================================

pcarcsurv.int <- 4  # ~ dnorm(5,1)               # probability of surviving as a carcass at day of death
pcarcsurv.slope <- 0 #~ dunif(-0.5,0)        # set to 0 for now- constant with age of carcass  # daily decrease in the probability of being observable since time of death
for(deadday in 1:ndays){  # loop past the last possible age
  logit(pcarcsurv[deadday]) <- pcarcsurv.int + pcarcsurv.slope*(deadday-1)   # deterministic for now
}

#=======================================
# set the scheduling distributions 
#=======================================

   # BIRTH DISTRIBUTION
mu.bday ~ dunif(10,70) #dnorm(45,0.01)T(20,60)      ## Prior for mean of birth distribution (1 would mean born just before day 0 [t=1], 10 would mean born just before day 9 [t=10])
sd.bday ~ dunif(1,20)                 ## Prior for standard deviation of birth distribution
prec.bday <- pow(sd.bday,-2)
nconst_Pb <- 1 - ( pnorm(0,mu.bday, prec.bday) + (1 - pnorm((max.day.birth),mu.bday, prec.bday) ) )  # constant for truncating the births distribution
Pb[1] <- 0  # can't be born after 0 time has passed in the season
Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, prec.bday) - pnorm(0:(max.day.birth-1), mu.bday, prec.bday))/nconst_Pb   ## Birth schedule

   # TIME ON COLONY
   # define the probability mass of remaining on the colony for a certain number of days (having left the colony between the previous day and the current day)
mu.days.on.colony ~ dunif(min.age, max.age)       ## Prior for mean of the number of days on the colony
sd.days.on.colony <- 7  #~ dunif(1, 20)        ## Prior for standard deviation of days on colony
prec.days.on.colony <- pow(sd.days.on.colony,-2)  
nconst_birth.to.leave <- 1 - ( pnorm(min.age,mu.days.on.colony, prec.days.on.colony) + (1 - pnorm(max.age,mu.days.on.colony, prec.days.on.colony) ) )
for(day in 1:(min.age+1)){
  pmass.birth.to.leave[day] <- 0   # can not leave before minimum number of days
  cprob.remain[day] <- 1
}
for(day in (min.age+2):(max.age+1)){   # 1 is day of birth, 2 is 1 day after birth, so on
  pmass.birth.to.leave[day] <- (pnorm(day-1, mu.days.on.colony, prec.days.on.colony) -  pnorm(day-2, mu.days.on.colony, prec.days.on.colony))/nconst_birth.to.leave      # define probability of leaving on any given age in days
  cprob.remain[day] <- 1-sum(pmass.birth.to.leave[1:day])  # cumulative probability of remaining on the colony a specific number of days after birth
}
for(day in (max.age+2):(ndays+1)){   
  pmass.birth.to.leave[day] <- 0  # no probability of leaving after specifically this length of time
  cprob.remain[day] <- 0  # no probability of remaining after this length of time
}

    # SURVIVAL OF CARCASSES WITH TIME SINCE DEATH

for(deadday in 1:(ndays+1)){    
  carcsurv.component[1,deadday] <- 0  # can't die at cohort age 0 (zero days after birth) 
}

for(day in 2:(ndays+1)){   # loop through days of cohort, with 1 representing time of birth
  carcsurv.component[day,1] <- 0  # can't die at day 0- haven't had an opportunity to die yet 
  carcsurv.component[day,2] <- cprob.remain[2]*1*(1-surv[1]) * prod(pcarcsurv[1:(day-1)])  # to be surviving carcass, having died in the first day of life, you need to remain on colony,be born, and die between age 0 and 1, and be a surviving carcass at the time of death
  for(deadday in 3:day){  # loop through all possible days of death prior to the current cohort age 
    carcsurv.component[day,deadday] <- cprob.remain[deadday] * prod(surv[1:(deadday-2)]) * (1-surv[(deadday-1)]) * prod(pcarcsurv[(deadday-1):(day-1)])   # to be observable in the dead state, you need to remain on colony,survive to the previous day, and die in the current time step, and the carcass has to survive to remain observable
  }
  for(rest in (day+1):(ndays+1)){
    carcsurv.component[day,rest] <- 0
  }
} 

for(day in 1:(ndays+1)){
  carcsurv[day] <- sum(carcsurv.component[day,1:ndays])  # sum across all the ways of being dead and observable as dead 
}

#=======================
# place cohorts in the proper places over time
#=======================

nborn ~ dunif(100,50000)       ## prior for expected number of births (total across season) ### eventually: function of colony, year and latitude 
for(c in 1:(max.day.birth+1)){
  cohort.size[c] ~ dpois(nborn*Pb[c])
}

for(c in 1:(max.day.birth+1)){    # loop through cohorts
      # SET ALL PRE-BIRTH ABUNDANCES TO ZERO
  for(a in 1:(c-1)){    # loop through all days prior to birth
    for(s in 1:5){      # loop through all stages
      stageN_c[c,a,s] <- 0      # set to zero prior to birth. the second dimension is converting to absolute time (vs cohort time)
    }
  }
     # SET ABUNDANCES WHILE POTENTIALLY ON COLONY
  stageN_c[c,c,1] <- cohort.size[c] * 1 * 1 * 1     # probability of being white after 0 days
  stageN_c[c,c,2] <- cohort.size[c] * 1 * 1 * 0   # moulted after 0 days
  stageN_c[c,c,3] <- cohort.size[c] * 0 * 1       # successfully leaving the colony alive after 0 days 
  stageN_c[c,c,4] <- cohort.size[c] * carcsurv[1]    # dead and surviving carcass after 0 days
  for(a in 2:(ndays-c) ){ 
    stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:(a-1)]) * (1-pmoult[a])     # live, on colony, white (not moulted)
    stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:(a-1)])*(pmoult[a])   # live, on colony, moulted
    stageN_c[c,c+a-1,3] <- cohort.size[c] * (1-cprob.remain[a]) * prod(surv[1:min(max.age,(a-1))])   # live and successfully left the colony
    stageN_c[c,c+a-1,4] <- cohort.size[c] * carcsurv[a]  # dead, surviving on colony as dead
  }
  
  for(a in 1:(ndays-c)){
    stageN_c[c,c+a-1,5] <- cohort.size[c]-(stageN_c[c,c+a-1,1]+stageN_c[c,c+a-1,2]+stageN_c[c,c+a-1,3]+stageN_c[c,c+a-1,4])  # dead and decayed or otherwise removed 
  }
}
  
  # compute the expected number in each stage for each day of the season (final state process step)

for(t in 1:(ndays-1)){ # loop through days of the season
  stageN[t,1] <- sum(stageN_c[1:(max.day.birth+1),t,1])  # tot white across all cohorts
  stageN[t,2] <- sum(stageN_c[1:(max.day.birth+1),t,2])  # tot grey across all cohorts
  stageN[t,3] <- sum(stageN_c[1:(max.day.birth+1),t,3])  # tot leave across all cohorts
  stageN[t,4] <- sum(stageN_c[1:(max.day.birth+1),t,4])  # tot dead and survived carcasses across all cohorts
  stageN[t,5] <- sum(stageN_c[1:(max.day.birth+1),t,5])  # tot dead and gone (not survived carcass) across all cohorts
}

# ========================
# Observation process
# ========================

# set prior on detection and misclassification processes
Pobs <- 0.95 #~ dunif(0.9,0.99)
Pobs_dead <- 0.95  #~ dunif(0.9,0.99)    # set deterministic for now
PCMoult <- 0.84

# normal approximation to the binomial distribution

# for(obs in 1:(length(obsdays))) {
#   white.mu[obs] <- stageN[obsdays[obs],1] * Pobs
#   white.sd[obs] <- pow(stageN[obsdays[obs],1] * Pobs,.5)*(1-Pobs)
#   white.prec[obs] <- pow(white.sd[obs],-2)
#   moult.mu[obs] <- stageN[obsdays[obs],2] * Pobs
#   moult.sd[obs] <- pow(stageN[obsdays[obs],2] * Pobs,.5)*(1-Pobs)
#   moult.prec[obs] <- pow(moult.sd[obs],-2)
#   dead.mu[obs] <- stageN[obsdays[obs],4] * Pobs_dead
#   dead.sd[obs] <- pow(stageN[obsdays[obs],4] * Pobs_dead,.5)*(1-Pobs_dead)
#   dead.prec[obs] <- pow(dead.sd[obs],-2)
#   wObs[obs] ~ dnorm(white.mu[obs],white.prec[obs])
#   mObs[obs] ~ dnorm(moult.mu[obs],moult.prec[obs])
#   dObs[obs] ~ dnorm(dead.mu[obs],dead.prec[obs])
# }

# add misclassification and use binomial process: 

WW_MW <- stageN[,1]+stageN[,2]*(1-PCMoult)
MM <- stageN[,2]*PCMoult  #add misclassification to mObs

for(obs in 1:(length(obsdays))) {
 wObs[obs] ~ dbinom(Pobs, trunc(WW_MW[obsdays[obs]]))   # data node
 mObs[obs] ~ dbinom(Pobs, trunc(MM[obsdays[obs]]))   # data node
 dObs[obs] ~ dbinom(Pobs_dead, trunc(stageN[obsdays[obs],4]))   # data node
}

#no misclassification for moulteds
# for(obs in 1:(length(obsdays))) {
#  wObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],1]))   # data node
#  mObs[obs] ~ dbinom(Pobs, trunc(stageN[obsdays[obs],2]))   # data node
#  dObs[obs] ~ dbinom(Pobs_dead, trunc(stageN[obsdays[obs],4]))   # data node
# }



}
",file=jagsfile)


# data objects for jags -----------------

#Parameters monitored
params <- c( "mu.bday", "sd.bday","nborn",
             "mu.days.on.colony","sd.days.on.colony",
             # "surv.int","surv.slope","pwhite.int","pwhite.slope","pcarcsurv.int", #"pcarcsurv.slope",
             "pmass.birth.to.leave","cprob.remain","Pb",
             "stageN",
             # "Pobs","Pobs_dead",
             "cohort.size")
             

# Define data
data <- list(
  wObs = wObs,
  mObs = mObs,
  dObs = dObs,
  # amax=amax,
  ndays=ndays,
  max.age = max.age,
  min.age = min.age,
  max.day.birth = max.day.birth,
  obsdays=obsdays
)

inits <- function(){
  list(
    mu.bday = 40,
    sd.bday = 20,
    mu.days.on.colony = 30, 
    # sd.days.on.colony= 2,  
    nborn = 5000
    # surv.int = 3,
    # surv.slope=0.03,
    # pmoult.int=-12,
    # pmoult.slope = 0.5,
    # pcarcsurv.int= 4,
    # Pobs = 0.95,
    # Pobs_dead = 0.95
   # cohort.size  = 10000*Pb 
  )
}
# inits()
#######
# run jags

DEBUG=T

if(DEBUG){
  
  nc <- 1   ;   ni <- 5000   ;   nb <- 1000   ;   nt <- 1 ; na <-1000
  
  Hgpup <- jagsUI::jags(data,inits,params,jagsfile,n.chains=nc,n.adapt=na,               
                          n.iter=ni,n.burnin=nb,n.thin=nt, parallel = F) 
}else{
  
  nc <- 3   ;   ni <- 50000   ;   nb <- 20000   ;   nt <- 10 ; na <-2000
  
  system.time(
    Hgpup <- jagsUI::jags(data,inits,params,jagsfile,n.chains=nc,n.adapt=na,               
                          n.iter=ni,n.burnin=nb,n.thin=nt, parallel = T, n.cores = 3) 
  )
}

##########################
# save(Hgpup,file = sprintf("Nborn_Hgpup_%s.RData",Sys.Date()))

# samples <- Hgpup$samples
# 
# library(coda)
# 
# 
# sims <- Hgpup$sims.list
# 
# sims$Pb[1,]
# sims$born[1,]
# 
# plot(sims$Pb[1,])
# plot(sims$Ps4[1,])
# 
# sims$nborn


jagsUI::traceplot(Hgpup,"nborn")

jagsUI::traceplot(Hgpup,"mu.bday")
jagsUI::traceplot(Hgpup,"sd.bday")
jagsUI::traceplot(Hgpup,"mu.days.on.colony")
jagsUI::traceplot(Hgpup,"sd.days.on.colony")


# jagsUI::traceplot(Hgpup,"stageN[60,3]")
# jagsUI::traceplot(Hgpup,"leave[60]")
# jagsUI::traceplot(Hgpup,"surv.int")
# jagsUI::traceplot(Hgpup,"surv.slope")

# jagsUI::traceplot(Hgpup,"pwhite.int")
# jagsUI::traceplot(Hgpup,"pwhite.slope")

# jagsUI::traceplot(Hgpup,"pcarcsurv.int")
# jagsUI::traceplot(Hgpup,"Pobs")
# jagsUI::traceplot(Hgpup,"Pobs_dead")

thissim=1000

Hgpup$sims.list$stageN[thissim,,3]

plot(Hgpup$sims.list$Pb[thissim,],type="l")
# sum(Hgpup$sims.list$Pb[thissim,])

plot(Hgpup$sims.list$pmass.birth.to.leave[thissim,],type="l")
sum(Hgpup$sims.list$pmass.birth.to.leave[thissim,])

plot(Hgpup$sims.list$stageN[thissim,,1],type="l",ylim=c(0,max(Hgpup$sims.list$stageN[thissim,,])),xlab='day of season',ylab="# of pups") #whitepups present
lines(Hgpup$sims.list$stageN[thissim,,2],col="red") #moulteds present
lines(Hgpup$sims.list$stageN[thissim,,3],col="green") # successfully leaving the colony alive after 0 days 
lines(Hgpup$sims.list$stageN[thissim,,4],col="blue")  #dead and survived carcass
lines(Hgpup$sims.list$stageN[thissim,,5],col="brown") #dead and gone (carcass did not survive)
# Add a legend
legend("topleft", legend=c("White pups present", "Moulteds present", "Successfully leaving the colony alive after 0 days", "Dead and survived carcass", "Dead and gone (carcass did not survive)"), 
       col=c("black", "red", "green", "blue", "brown"), lty=1,cex=(0.6),inset = c(0, -0.4),xpd = TRUE)
apply(Hgpup$sims.list$stageN[thissim,,],1,sum)

Hgpup$sims.list$nborn[thissim]
sum(Hgpup$sims.list$cohort.size[thissim,])  # matches with max

plot(Hgpup$sims.list$stageN[thissim,,3],type="l")

Hgpup$sims.list$pwhite.int[thissim]
Hgpup$sims.list$pwhite.slope[thissim]


#???????????????
# day=50;sim=thissim
# whitefunc <- function(day,sim=thissim){
#   plogis(Hgpup$sims.list$pwhite.int[sim]+Hgpup$sims.list$pwhite.slope[sim]*day)
# }
# whitefunc(100,thissim)
# curve(whitefunc,0,100)

plot(Hgpup$sims.list$stageN[thissim,,3],type="l",col="green")

Hgpup$sims.list$nborn[thissim]

Hgpup$sims.list$surv.int[thissim]
Hgpup$sims.list$surv.slope[thissim] 

plot(Hgpup$sims.list$cprob.remain[thissim,],type="l")



