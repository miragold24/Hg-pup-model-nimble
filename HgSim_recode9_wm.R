

# for testing- alternative data generation
# goal: generate wObs and mObs under the process model

ndays=123
min.age=20
max.age=60
max.day.birth=100
obsdays=c(30,40,50,60,70)

mu.bday = 45
sd.bday = 10
nborn = 500
mu.days.on.colony = 31.5
sd.days.on.colony = 7
surv.int = 4
surv.slope = 0.1
pmoult.int = -12
pmoult.slope = 0.6
Pobs = 0.95
PCMoult = 0.84

# Pb <- numeric(max.day.birth)
# Pb[1] <- 0  # can't enter after 0 time has passed in the season
# Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))   ## Birth schedule
# plot(Pb,type="l")


Pb <- numeric(max.day.birth)
Pb[1] <- 0  # can't enter after 0 time has passed in the season
nconst_Pb <- 1 - ( pnorm(0,mu.bday, sd.bday) + (1 - pnorm((max.day.birth),mu.bday, sd.bday) ) )  # constant for truncating the births distribution
Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))/nconst_Pb   ## Birth schedule
plot(Pb,type="l")

# pmass.birth.to.leave <- numeric(max.age+1)
# cprob.remain <- numeric(max.age+1)
# for(day in 1:(min.age+1)){
#   pmass.birth.to.leave[day] <- 0   # can not leave before minimum number of days
#   cprob.remain[day] <- 1
# }
# for(day in (min.age+2):(max.age+1)){   # 1 is day of birth, 2 is 1 day after birth, so on
#   pmass.birth.to.leave[day] <- (pnorm(day-1, mu.days.on.colony, sd.days.on.colony) -  pnorm(day-2, mu.days.on.colony, sd.days.on.colony))      # define probability of leaving on any given age in days
#   cprob.remain[day] <- 1-sum(pmass.birth.to.leave[1:day])  # cumulative probability of remaining in pres4 a specific number of days after becoming s4
# }



nconst_birth.to.leave <- 1 - ( pnorm(min.age,mu.days.on.colony, sd.days.on.colony) + (1 - pnorm(max.age,mu.days.on.colony, sd.days.on.colony) ) )
for(day in 1:(min.age+1)){
  pmass.birth.to.leave[day] <- 0   # can not leave before minimum number of days
  cprob.remain[day] <- 1
}
for(day in (min.age+2):(max.age+1)){   # 1 is day of birth, 2 is 1 day after birth, so on
  pmass.birth.to.leave[day] <- (pnorm(day-1, mu.days.on.colony, sd.days.on.colony) -  pnorm(day-2, mu.days.on.colony, sd.days.on.colony))/nconst_birth.to.leave      # define probability of leaving on any given age in days
  cprob.remain[day] <- 1-sum(pmass.birth.to.leave[1:day])  # cumulative probability of remaining on the colony a specific number of days after birth
}


cohort.size <- rmultinom(1,nborn,Pb)[,1]
sum(cohort.size)

surv <- plogis(surv.int+surv.slope*(0:max.age))
plot(surv,type="l")

pmoult <- plogis(pmoult.int+pmoult.slope*(0:max.age))
plot(pmoult,type="l",xlim=c(0,40))


stageN_c <- array(0,dim=c(max.day.birth+1,ndays,3))
c=1
for(c in 1:(max.day.birth+1)){    # loop through cohorts
  # SET ABUNDANCES WHILE POTENTIALLY ON COLONY
  stageN_c[c,c,1] <- cohort.size[c] * 1 * 1 * 1     # probability of being white after a days
  stageN_c[c,c,2] <- cohort.size[c] * 1 * 1 * 0   # grey
  # stageN_c[c,c,3] <- (nborn*Pb[c])-(stageN_c[c,c,1]+stageN_c[c,c,2])  # dead or left
  for(a in 2:min( (ndays-c),(max.age+1) ) ){ 
    stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a]) *(1-pmoult[a])   # probability of being white after a days
    stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a])*  pmoult[a]    #moulted
  }
  for(a in 1:(ndays-c)){
    stageN_c[c,c+a-1,3] <- cohort.size[c]-(stageN_c[c,c+a-1,1]+stageN_c[c,c+a-1,2])  # dead or left
  }
}

plot(stageN_c[50,,3],type="l")

stageN <- array(0,dim=c((ndays-1),3))
for(t in 1:(ndays-1)){ # loop through days of the season
  stageN[t,1] <- sum(stageN_c[1:(max.day.birth+1),t,1])  # tot white across all cohorts
  stageN[t,2] <- sum(stageN_c[1:(max.day.birth+1),t,2])  # tot grey across all cohorts
  stageN[t,3] <- sum(stageN_c[1:(max.day.birth+1),t,3])  # tot leave across all cohorts
}

plot(stageN[,1],type="l")
plot(stageN[,2],type="l")
plot(stageN[,3],type="l")

wObs <- numeric(length(obsdays))
mObs <- numeric(length(obsdays))

# add misclassification of moulteds as whitecoats and use binomial process: 
WW_MW <- stageN[,1]+stageN[,2]*(1-PCMoult)
MM <- stageN[,2]*PCMoult  #add misclassification to mObs

for(obs in 1:(length(obsdays))) {
  wObs[obs] <- rbinom(1,trunc(WW_MW[obsdays[obs]]),Pobs )   # data node
  mObs[obs] <- rbinom(1,trunc(MM[obsdays[obs]]),Pobs )   # data node
}


# for(obs in 1:(length(obsdays))) {
#  wObs[obs] <- rbinom(1,trunc(stageN[obsdays[obs],1]),Pobs )   # data node
#  mObs[obs] <- rbinom(1,trunc(stageN[obsdays[obs],2]),Pobs )   # data node
# }



#plot expected whitecoats and moulteds vs observed
# plot(stageN[,1],type="l",ylab='# of pups',xlab='days of season')
# lines(WW_MW,col='blue')
# lines(stageN[,2],type="l")
# lines(MM,col='red')
# legend("topright", legend=c("Expected white", "Expected moulted", "Observed white", "Observed moulted"), 
#        col=c("black", "blue", "black", "red",cex=(0.4)))#,inset = c(0, -0.4),xpd = TRUE)


