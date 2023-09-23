
###==========================================
###
### Input data & functions for HgNimble_V3.R
### to calculate pups present in each class on a given day
###==========================================
input_list = list(ndays = 123, amax=123,N=500,
                  start = 32, mu.bday = 45, sd.bday = 10,
                  mu.mlen = 6, sd.mlen = 2,mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.94, POMoultSD = 0.03,
                  PCWhite = 1, 
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  obsdays = c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115),
                 seed=2024)
thisenv <-environment()
list2env(input_list, thisenv)

####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================

whiteSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1),whiteObs=double(1)){
    # whitePups <- numeric(value = 0, length = 123)
    ndays=length(born)
    obsdays=c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115)
    whitePups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age <- d-i
        whitePups[d] <- whitePups[d] + born[i] *(1 - sum(Pm_sum[1:age])) * (1 - sum(Pl_sum[1:age]))
      }
    }
    
    # Constrain whitePups to be greater than wObs ? 
   # for (obs in 1:length(obsdays)){
   #   if (whitePups[obsdays[obs]] <= whiteObs[obs]) {
   #     whitePups[obsdays[obs]] <- whiteObs[obs] + 1  # Adjust whitePups to meet the constraint
   #   }
   # }
    returnType(double(1))
    return(whitePups)
  }
)
# whiteSumc <- compileNimble(whiteSum)
# born=model$born
# Pm_sum=model$Pm_sum
# Pl_sum=model$Pl_sum
#white<-whiteSum(born,Pm_sum,Pl_sum)

moultSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1),moultObs=double(1)){
    ndays=length(born)
    obsdays=c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115)
    moultPups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        moultPups[d] <- moultPups[d] + born[i] * sum(Pm_sum[1:age]) * (1-sum(Pl_sum[1:age]))
      }
    }
    
    # Constrain whitePups to be greater than wObs
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
# wCount <- nimbleFunction(
#   run = function(born=double(1),whiteCount=doulbe(1)){
#     # whitePups <- numeric(value = 0, length = 123)
#     ndays=length(born)
#     born=sim_S3$realBorn
#     obsdays=c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115)
#     whiteCount <- rep(0,ndays)
#     # Constrain whitePups to be greater than wObs ? 
#     for (obs in 1:length(obsdays)){
#       if (whiteCount[obsdays[obs]] <= sim_S3$whiteCount[obs]) {
#         whiteCount[obsdays[obs]] <- sim_S3$whiteCount[obs] + 1  # Adjust whitePups to meet the constraint
#       }
#     }
#     returnType(double(1))
#     return(whiteCount)
#   }
# )
# w<-wCount(born,whiteCount)
# w
