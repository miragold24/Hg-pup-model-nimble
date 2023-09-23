
###==========================================
###
### Input data & functions for HgNimble_V3.R
###
###==========================================
library(nimble)
input_list = list(ndays = 123, amax=123,N=500,
                  start = 32, mu.bday = 45, sd.bday = 10,
                  mu.s4day = 17, sd.s4day = 5,
                  mu.mlen = 6, sd.mlen = 2,mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.94, POMoultSD = 0.03,
                  PCWhite = 1,
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
                  obsdays = c(35, 40, 50, 65, 70, 80, 95))
thisenv <-environment()
list2env(input_list, thisenv)
#convert params for log scale
lCV <- sd.bday/(mu.bday-start)
lsd <- sqrt(log(lCV^2+1))
lmu <- log(mu.bday-start)-(0.5*(lsd^2))
####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================


whiteSum <- nimbleFunction(
  run = function(born=double(1),PS4_sum=double(1),Pm_sum=double(1),Pl_sum=double(1)){
    # whitePups <- numeric(value = 0, length = 123)
    ndays=length(born)
    whitePups <- rep(0,ndays)
    
    for (d in 2:(ndays)){
      for (i in 1:d){
        age <- d-i
        whitePups[d] <- whitePups[d] + born[i] * (1-sum(PS4_sum[1:age])) *(1 - sum(Pm_sum[1:age])) * (1 - sum(Pl_sum[1:age]))
      }
    }
    returnType(double(1))
    return(round(whitePups))
  }
)
#whiteSumc <- compileNimble(whiteSum)

#white_sum<-whiteSum(born,PS4,Pm,Pl



S4Sum <- nimbleFunction(
  run = function(born=double(1),PS4_sum=double(1),Pm_sum=double(1),Pl_sum=double(1)){
    ndays=length(born)
    S4Pups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        S4Pups[d] <- S4Pups[d] + born[i]  * sum(PS4_sum[1:age])  * (1-sum(Pm_sum[1:age])) * (1-sum(Pl_sum[1:age]))
      }
    }
    returnType(double(1))
    return(round(S4Pups))
  }
)
#
# S4<-S4Sum(born,PS4_sum,Pm_sum,Pl_sum)
#S4Sumc <- compileNimble(S4Sum)
#test<-S4Sumc(born,PS4,Pm,Pl)


moultSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1)){
    ndays=length(born)
    moultPups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        moultPups[d] <- moultPups[d] + born[i] * sum(Pm_sum[1:age]) * (1-sum(Pl_sum[1:age]))
      }
    }
    returnType(double(1))
    return(round(moultPups))
  }
)


