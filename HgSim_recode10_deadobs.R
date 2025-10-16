

# for testing- alternative data generation
# goal: generate wObs and mObs under the process model

ndays=123
min.age=20
max.age=60
max.day.birth=100
obsdays=c(30,40,50,70,80)#,90,100,110)

mu.bday = 45
sd.bday = 10
nborn = 1000
mu.days.on.colony = 35
sd.days.on.colony = 7
surv.int = 4
surv.slope = 0.1
pmoult.int <- -12            # probability of being moult at birth
pmoult.slope <- 0.6
pcarcsurv.int = 4
pcarcsurv.slope=0
Pobs = 0.95
Pobs_dead = 0.95
PCMoult=0.84

Pb <- numeric(max.day.birth)
Pb[1] <- 0  # can't enter after 0 time has passed in the season
Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))   ## Birth schedule
plot(Pb,type="l")

pmass.birth.to.leave <- numeric(max.age+1)
cprob.remain <- numeric(max.age+1)
for(day in 1:(min.age+1)){
  pmass.birth.to.leave[day] <- 0   # can not leave before minimum number of days
  cprob.remain[day] <- 1
}
for(day in (min.age+2):(max.age+1)){   # 1 is day of birth, 2 is 1 day after birth, so on
  pmass.birth.to.leave[day] <- (pnorm(day-1, mu.days.on.colony, sd.days.on.colony) -  pnorm(day-2, mu.days.on.colony, sd.days.on.colony))      # define probability of leaving on any given age in days
  cprob.remain[day] <- 1-sum(pmass.birth.to.leave[1:day])  # cumulative probability of remaining in white a specific number of days
}
pmass.birth.to.leave <- c(pmass.birth.to.leave,rep(0,(ndays+1)-(max.age+1)))
cprob.remain <- c(cprob.remain,rep(0,(ndays+1)-(max.age+1)))

cohort.size <- rmultinom(1,nborn,Pb)[,1]
sum(cohort.size)
plot(cohort.size)

surv <- plogis(surv.int+surv.slope*(0:(ndays+1))) #cumulative daily probability of survival (low in beginning, high as you get older)
plot(surv,type="l")

pmoult <- plogis(pmoult.int+pmoult.slope*(0:(ndays+1)))
plot(pmoult,type="l",xlim=c(0,40))

pcarcsurv <- plogis(pcarcsurv.int + pcarcsurv.slope*(0:ndays)) # # set to 0 for now- constant with age of carcass 
plot(pcarcsurv,type="l",xlim=c(0,40))
# observable dead state

carcsurv.component <- array(0,dim=c((ndays+1),(ndays+1)))
for(deadday in 1:(ndays+1)){
  carcsurv.component[1,deadday] <- 0  # can't die at cohort age 0 (zero days after birth) 
}
day=2
for(day in 2:(ndays+1)){   # loop through days of cohort, with 1 representing time of birth
  carcsurv.component[day,1] <- 0  # can't die at day 0- haven't had an opportunity to die yet 
  carcsurv.component[day,2] <- cprob.remain[2]*1*(1-surv[1]) * prod(pcarcsurv[1:(day-1)])  # to be observable dead, having died in the first day of life, you need to remain on colony,survive to age 0, and die between age 0 and 1, and be observable at the time of death
  for(deadday in 3:day){  # loop through all possible days of death prior to the current cohort age 
    carcsurv.component[day,deadday] <- cprob.remain[deadday] * prod(surv[1:(deadday-2)]) * (1-surv[(deadday-1)]) * prod(pcarcsurv[(deadday-1):(day-1)])   # to be observable in the dead state, you need to remain on colony,survive to the previous day, and die in the current time step, and the carcass has to survive to remain observable
  }
} 

carcsurv <- numeric((ndays+1))
for(day in 1:(ndays+1)){
  carcsurv[day] <- sum(carcsurv.component[day,1:ndays])  # sum across all the ways of being dead and observable as dead 
}

stageN_c <- array(0,dim=c(max.day.birth+1,ndays,5))
c=1
for(c in 1:(max.day.birth+1)){    # loop through cohorts
  # SET ABUNDANCES WHILE POTENTIALLY ON COLONY
  stageN_c[c,c,1] <- cohort.size[c] * 1 * 1 * 1     # probability of being white after a days
  stageN_c[c,c,2] <- cohort.size[c] * 1 * 1 * 0   # molted
  stageN_c[c,c,3] <- cohort.size[c] * 0 * 1       # successfully leaving the colony alive after 0 days 
  stageN_c[c,c,4] <- cohort.size[c] * carcsurv[1]    # dead and on colony and observable dead after 0 days
  # stageN_c[c,c,3] <- (nborn*Pb[c])-(stageN_c[c,c,1]+stageN_c[c,c,2])  # dead or left
  for(a in 2:(ndays-c) ){ 
    stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:(a-1)]) * (1-pmoult[a])     # live, on colony, white
    stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:(a-1)])*pmoult[a]  # live, on colony, moulted
    stageN_c[c,c+a-1,3] <- cohort.size[c] * (1-cprob.remain[a]) * prod(surv[1:min(max.age,(a-1))])   # live and successfully left the colony
    stageN_c[c,c+a-1,4] <- cohort.size[c] * carcsurv[a]  # dead, observable as dead
  }
  for(a in 1:(ndays-c)){
    stageN_c[c,c+a-1,5] <- cohort.size[c]-(stageN_c[c,c+a-1,1]+stageN_c[c,c+a-1,2]+stageN_c[c,c+a-1,3]+stageN_c[c,c+a-1,4])  # unobs dead state
  }
}

plot(stageN_c[50,,3],type="l")

stageN <- array(0,dim=c((ndays-1),5))
for(t in 1:(ndays-1)){ # loop through days of the season
  stageN[t,1] <- sum(stageN_c[1:(max.day.birth+1),t,1])  # tot white across all cohorts
  stageN[t,2] <- sum(stageN_c[1:(max.day.birth+1),t,2])  # tot grey across all cohorts
  stageN[t,3] <- sum(stageN_c[1:(max.day.birth+1),t,3])  # tot leave across all cohorts
  stageN[t,4] <- sum(stageN_c[1:(max.day.birth+1),t,4])  # tot dead and observable across all cohorts
  stageN[t,5] <- sum(stageN_c[1:(max.day.birth+1),t,5])  # tot dead and unobservable across all cohorts
}

plot(stageN[,1],type="l")
plot(stageN[,2],type="l")
plot(stageN[,3],type="l")
plot(stageN[,4],type="l")
plot(stageN[,5],type="l")


wObs <- numeric(length(obsdays))
mObs <- numeric(length(obsdays))
dObs <- numeric(length(obsdays))
#no misclassification for moulted as whitecoats
# for (obs in 1:length(obsdays)) {
#   wObs[obs] <- rbinom(1, trunc(stageN[obsdays[obs], 1]), Pobs)
#   mObs[obs] <- rbinom(1, trunc(stageN[obsdays[obs], 2]), Pobs)
#   dObs[obs] <- rbinom(1, trunc(stageN[obsdays[obs], 4]), Pobs_dead)
# }

# add misclassification of moulteds as whitecoats and use binomial process: 
WW_MW <- stageN[,1]+stageN[,2]*(1-PCMoult)
MM <- stageN[,2]*PCMoult  #add misclassification to mObs

for(obs in 1:(length(obsdays))) {
 wObs[obs] <- rbinom(1,trunc(WW_MW[obsdays[obs]]),Pobs )   # data node
 mObs[obs] <- rbinom(1,trunc(MM[obsdays[obs]]),Pobs )   # data node
 dObs[obs] <- rbinom(1,trunc(stageN[obsdays[obs],4]),Pobs_dead )   # data node
}


#make sure pups sum to total nborn
#apply(stageN[,],1,sum)
