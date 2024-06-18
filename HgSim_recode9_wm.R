
# for testing- alternative data generation
# goal: generate wObs and mObs under the process model

ndays = 123
min.age = 20
max.age = 60
max.day.birth = 100
#obsdays = seq(20,100,20)  # list of different obsdays combinations
mu.bday = 45
sd.bday = 10
skew.bday=5
nborn = 700
mu.days.on.colony = 31.5
sd.days.on.colony = 2
surv.int = 4
surv.slope = 0.1
pmoult.int = -12
pmoult.slope = 0.6
Pobs = 0.95
PCMoult=0.84

generate_obs <- function(
    ndays = 123,
    min.age = 20,
    max.age = 60,
    max.day.birth = 100,
    obsdays = seq(40,90,12),  # list of different obsdays combinations
    mu.bday = 45,
    sd.bday = 10,
    skew.bday=5,
    nborn = 700,
    mu.days.on.colony = 31.5,
    sd.days.on.colony = 2,
    surv.int = 4,
    surv.slope = 0.1,
    pmoult.int = -12,
    pmoult.slope = 0.6,
    Pobs = 0.95,
    PCMoult=0.84
) {


dskewnorm <- function(x,skew.bday,mu.bday,sd.bday) {
  Pb <- rep(0,ndays)  
  delta <- skew.bday/sqrt(1+(skew.bday^2))
  
  omega <- sqrt(sd.bday^2/(1-(2*delta^2)/pi))
  
  xi <- mu.bday - omega*delta*sqrt(2/pi)
  Pb <- (2/omega)*dnorm((x-xi)/omega) * pnorm(skew.bday*((x-xi)/omega))
  
  # Return the probability vector
  return(Pb)
}


# # Integrate over days of season for prob of being born on each day
integratePb <-  function(days,skew.bday,mu.bday,sd.bday) {
#   # Initialize Pb as a parameter of the same length as ndays
  Pb <- rep(0, length(days))
  for (t in 2:length(days)) {
    # Define the start and end of the interval
    start <- t-1
    end <-  t

    # Integrate over the interval and add to the probability for this day
    Pb[t] <- Pb[t] + integrate(dskewnorm,skew.bday=skew.bday,mu.bday=mu.bday,sd.bday=sd.bday, start, end)$value
  }

  return(Pb)
}
# 
# #skew normal birth distribution
pb<-rep(0,(max.day.birth+1))
Pb<-rep(0,(max.day.birth+1))
pb[1:(max.day.birth+1)] <- integratePb(1:(max.day.birth+1), skew.bday, mu.bday, sd.bday)
Pb[1:(max.day.birth+1)] <- pb[1:(max.day.birth+1)]/sum(pb[1:(max.day.birth+1)])



## normal birth distribution
# Pb <- numeric(max.day.birth)
# Pb[1] <- 0  # can't enter after 0 time has passed in the season
# nconst_Pb <- 1 - ( pnorm(0,mu.bday, sd.bday) + (1 - pnorm((max.day.birth),mu.bday, sd.bday) ) )  # constant for truncating the births distribution
# Pb[2:(max.day.birth+1)] <- (pnorm(1:(max.day.birth), mu.bday, sd.bday) - pnorm(0:(max.day.birth-1), mu.bday, sd.bday))/nconst_Pb   ## Birth schedule
# plot(Pb,type="l")



pmass.birth.to.leave <- numeric(max.age+1)
cprob.remain <- numeric(max.age+1)
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
#plot(surv,type="l")

pmoult <- plogis(pmoult.int+pmoult.slope*(0:max.age))
#plot(pmoult,type="l",xlim=c(0,40))


stageN_c <- array(0,dim=c(max.day.birth+1,ndays,3))
c=1
for(c in 1:(max.day.birth+1)){    # loop through cohorts
  # SET ABUNDANCES WHILE POTENTIALLY ON COLONY
  stageN_c[c,c,1] <- cohort.size[c] * 1 * 1 * 1     # probability of being white after a days
  stageN_c[c,c,2] <- cohort.size[c] * 1 * 1 * 0   # moulted
  # stageN_c[c,c,3] <- (nborn*Pb[c])-(stageN_c[c,c,1]+stageN_c[c,c,2])  # dead or left
  for(a in 2:min( (ndays-c),(max.age+1) ) ){ 
    stageN_c[c,c+a-1,1] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a]) *(1-pmoult[a])   # probability of being white after a days
    stageN_c[c,c+a-1,2] <- cohort.size[c] * cprob.remain[a] * prod(surv[1:a])*  pmoult[a]    #moulted
  }
  for(a in 1:(ndays-c)){
    stageN_c[c,c+a-1,3] <- cohort.size[c]-(stageN_c[c,c+a-1,1]+stageN_c[c,c+a-1,2])  # dead or left
  }
}

#plot(stageN_c[50,,3],type="l")

stageN <- array(0,dim=c((ndays-1),3))
for(t in 1:(ndays-1)){ # loop through days of the season
  stageN[t,1] <- sum(stageN_c[1:(max.day.birth+1),t,1])  # tot white across all cohorts
  stageN[t,2] <- sum(stageN_c[1:(max.day.birth+1),t,2])  # tot grey across all cohorts
  stageN[t,3] <- sum(stageN_c[1:(max.day.birth+1),t,3])  # tot leave across all cohorts
}

# plot(stageN[,1],type="l")
# plot(stageN[,2],type="l")
# plot(stageN[,3],type="l")

# results <- list()

# for (obsdays in obsdays_list) {
wObs <- numeric(length(obsdays))
mObs <- numeric(length(obsdays))

# add misclassification of moulteds as whitecoats and use binomial process: 
WW_MW <- stageN[,1]+stageN[,2]*(1-PCMoult)
MM <- stageN[,2]*PCMoult  #add misclassification to mObs
obs=1
for(obs in 1:(length(obsdays))) {
  wObs[obs] <- rbinom(1,trunc(WW_MW[obsdays[obs]]),Pobs )   # data node
  mObs[obs] <- rbinom(1,trunc(MM[obsdays[obs]]),Pobs )   # data node
}

#plot expected whitecoats and moulteds vs observed
# plot(stageN[,1],type="l",ylab='# of pups',xlab='days of season')#,main='normal')
# lines(WW_MW,col='blue')
# lines(stageN[,2],type="l")
# lines(MM,col='red')
# legend("topright", legend=c("Expected white", "Expected moulted", "Observed white", "Observed moulted"),
#        col=c("black", "blue", "black", "red",cex=(0.4)))#,inset = c(0, -0.4),xpd = TRUE)
results <- list(wObs = wObs, mObs = mObs,Pb=Pb)
# }

return(results)
}

# obsdays_combinations <- list(
#   #c(35,45, 55, 65, 75,85)
#   #seq(10, 50, by = 10),  # Observation days starting at 10, going to 50, spaced 10 days apart
#   seq(50, 90, by = 10),
#   seq(60, 100, by = 10)
# )


sim_wm <- generate_obs()

