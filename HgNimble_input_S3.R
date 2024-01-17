library(nimble)
###==========================================
###
### Input data & functions for HgNimble_V3.R
### to calculate pups present in each class on a given day
###==========================================
input_list = list(ndays = 123, amax=123,N=500,
                  start = 32, mu.bday = 45, sd.bday = 10,
                  mu.mlen = 6, sd.mlen = 2,
                  mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  mu.s4day=17,sd.s4day=4,
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.95, POMoultSD = 0.03,
                  PCWhite = 1, PCWhiteS4=0.75,
                  PCWhiteS4SD=0.01,
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  obsdays = c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110),
                  seed=2024)
thisenv <-environment()
list2env(input_list, thisenv)

####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================
# convert birth parameters for lognormal dist
lCV <- sd.bday/(mu.bday-start)
lsd <- sqrt(log(lCV^2+1))
lmu <- log(mu.bday-start)-(0.5*(lsd^2))

# calculate probability of being born on each day
Pb <- rep(0, ndays)

for(t in 2:ndays){
  Pb[t] <- plnorm((t-start+1), lmu, lsd) - plnorm((t-start), lmu, lsd)
}


whiteSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1)){
    ndays=length(born)
    obsdays=c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)
    whitePups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age <- d-i
        #whitePups[d] <-whitePups[d] + born[i] * (1 - sum(Pm_sum[1:age]))
        whitePups[d] <- whitePups[d] + born[i] *(1 - sum(Pm_sum[1:age])) #* (1 - sum(Pl_sum[1:age]))
      }
    }
    
    # Constrain whitePups to be greater than wObs 
   # for (obs in 1:length(obsdays)){
   #   if (whitePups[obsdays[obs]] <= whiteObs[obs]) {
   #     whitePups[obsdays[obs]] <- whiteObs[obs] + 1  # Adjust whitePups to meet the constraint
   #   }
   # }

   #}
    returnType(double(1))
    return(whitePups)
  }
)
    
    
    
    S4Sum <- nimbleFunction(
      run = function(born=double(1),PS4_sum=double(1),Pm_sum=double(1)){
        ndays=length(born)
        obsdays=c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)
        S4Pups <- rep(0,ndays)
        for (d in 2:(ndays)){
          for (i in 1:d){
            age = d-i
            S4Pups[d] <- S4Pups[d] + born[i]  * sum(PS4_sum[1:age])  * (1-sum(Pm_sum[1:age])) # * (1-sum(Pl_sum[1:age]))
          }
        }
        returnType(double(1))
        return(round(S4Pups))
  }
)


moultSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1)){
    ndays=length(born)
    obsdays=c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)
    moultPups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        moultPups[d] <- moultPups[d] + born[i] * sum(Pm_sum[1:age]) * (1-sum(Pl_sum[1:age]))
      }
    }
    
    #Constrain whitePups to be greater than wObs
    # for (obs in 1:length(obsdays)){
    #   if (moultPups[obsdays[obs]] <= moultObs[obs]) {
    #     moultPups[obsdays[obs]] <- moultObs[obs] + 1  # Adjust whitePups to meet the constraint
    #   }
    # }
    
    returnType(double(1))
    return(moultPups)
  }
)
 #moultSumc <-compileNimble(moultSum)
# 
#   wCount <- nimbleFunction(
#     run = function(whiteObs=double(1),ww=double(1),mw=double(1)){
#     # whitePups <- numeric(value = 0, length = 123)
#     obsdays=c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115)
#     whiteCount <- rep(0,length(obsdays))
# #     # Constrain whitePups to be greater than wObs ? 
#     for (obs in 1:length(obsdays)){
#     whiteCount[obs] = ww[obs]+mw[obs]
#     if (whiteCount[obs] <= whiteObs[obs]) {
#     whiteCount[obs] <- whiteObs[obs] + 1  # Adjust whiteCount to meet the constraint
#     }
#   
# }
#     returnType(double(1))
#     return(whiteCount)
#   }
# )
# # w<-wCount(born,whiteCount)
# # w
# 
#   mCount <- nimbleFunction(
#     run = function(moultObs=double(1),mm=double(1),wm=double(1)){
#       # whitePups <- numeric(value = 0, length = 123)
#       obsdays=c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115)
#       moultCount <- rep(0,length(obsdays))
#       #     # Constrain whitePups to be greater than wObs ? 
#       for (obs in 1:length(obsdays)){
#         moultCount[obs] = mm[obs]+wm[obs]
#         if (moultCount[obs] <= moultObs[obs]) {
#           moultCount[obs] <- moultObs[obs] + 1  # Adjust moultCount to meet the constraint
#         }
#         
#       }
#       returnType(double(1))
#       return(moultCount)
#     }
#   )
  
