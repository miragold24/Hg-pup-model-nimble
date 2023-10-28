library(nimble)
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
                  POMoult = 0.95, POMoultSD = 0.03,
                  PCWhite = 1, 
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  skew.bday = 5, 
                  obsdays = c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110),
                 seed=2024)
thisenv <-environment()
list2env(input_list, thisenv)

####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================



#### will brms load into nimble function?
PbSkew<-function(mu.bday,sd.bday,skew.bday){ 
  library(brms)
  skew.bday=5
  Pb <- rep(0, ndays)
  
  for(t in 2:ndays){
    Pb[t] <- pskew_normal(q = t, mu = mu.bday, sigma = sd.bday, alpha = skew.bday) - 
      pskew_normal(q = t-1, mu = mu.bday, sigma = sd.bday, alpha = skew.bday)
  }
  
  Pb <- Pb/sum(Pb)
  return(Pb)
}

PbSkew(mu.bday,sd.bday,skew.bday)
NimSimPbSkew <- nimbleRcall(
  prototype=function(mu.bday=double(1),sd.bday=double(1),skew.bday=double(1)
  ){},
  
  returnType=double(1),
  Rfun = 'PbSkew'
  
)
# NimSimPbSkew(mu.bday,sd.bday,skew.bday)
# 
# # using Pb skewnormal to check nborn estimation in nimble (no other priors will be used bc calculated here
# #not in model)
# PbSkew<-function(Pb){  
#   Pb=simseals$Pb
#   Pb <- Pb/sum(Pb)
#   return(Pb)
# }
# 
# PbSkew()
# NimSimPbSkew <- nimbleRcall(
#   prototype=function(Pb=double(1)
#   ){},
#   
#   returnType=double(1),
#   Rfun = 'PbSkew'
#   
# )
#NimSimPbSkew()



whiteSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1),whiteObs=double(1)){
    ndays=length(born)
    obsdays=c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)
    whitePups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age <- d-i
        whitePups[d] <- whitePups[d] + born[i] *(1 - sum(Pm_sum[1:age])) * (1 - sum(Pl_sum[1:age]))
      }
    }
    
    # Constrain whitePups to be greater than wObs 
   for (obs in 1:length(obsdays)){
     if (whitePups[obsdays[obs]] <= whiteObs[obs]) {
       whitePups[obsdays[obs]] <- whiteObs[obs] + 1  # Adjust whitePups to meet the constraint
     }
   }

   #}
    returnType(double(1))
    return(whitePups)
  }
)



moultSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1),moultObs=double(1)){
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
    for (obs in 1:length(obsdays)){
      if (moultPups[obsdays[obs]] <= moultObs[obs]) {
        moultPups[obsdays[obs]] <- moultObs[obs] + 1  # Adjust whitePups to meet the constraint
      }
    }
    
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
  
  
  
  ### properly calculate PWhite and PMoult for Observation process 
  
  
  # pWhite<- nimbleFunction(
  #   run = function(born=double(1),Pb_sum=double(1),Pm_sum=double(1),Pl_sum=double(1)){
  #     ndays=length(born)
  #     for (d in 2:(ndays)){
  #       for (i in 1:d){
  #         age <- d-i
  #         PWhite[d] <- Pb_sum[i]*(1 - sum(Pm_sum[1:age])) * (1 - sum(Pl_sum[1:age])*
  #                       POWhite*PCWhite)+(sum(Pm_sum[1:age]))*(1 - sum(Pl_sum[1:age]))*
  #                       POMoult*(1-PCMoult)
  #       }
  #     }
  #     returnType(double(1))
  #     return(PWhite)
  #   }
  # )
  # born=model$born
  # Pb_sum=model$Pb_sum
  # Pm_sum=model$Pm_sum
  # Pl_sum=model$Pl_sum
  # pwhite<-pWhite(born,Pm_sum,Pl_sum)
  # 
  
  
#   pMoult<- nimbleFunction(
#     run = function(born=double(1),Pb_sum=double(1),Pm_sum=double(1),Pl_sum=double(1)){
#       ndays=length(born)
#       for (d in 2:(ndays)){
#         for (i in 1:d){
#           age <- d-i
#           PMoult[d] <- Pb_sum[i]*(sum(Pm_sum[1:age])) * (1 - sum(Pl_sum[1:age]))*
#                         POMoult*PCMoult+(1-sum(Pm_sum[1:age]))*
#                         (1 - sum(Pl_sum[1:age]))*POWhite*(1-PCWhite)
#         }   
#         
#         ###
#         ###  Probability of moulting across season
#         ##
#         pmoult_check<-Pb_sum*(Pm_sum) *(1 - Pl_sum)*
#           POMoult*PCMoult+(1-Pm_sum)*
#           (1 - Pl_sum)*POWhite*(1-PCWhite)
#       }
#       returnType(double(1))
#       return(PMoult)
#     }
#   )
#   born=model$born
#   Pb_sum=model$Pb_sum
#   Pm_sum=model$Pm_sum
#   Pl_sum=model$Pl_sum
#   pmoult<-pMoult(born,Pb_sum,Pm_sum,Pl_sum)
# lines(pmoult)
#   