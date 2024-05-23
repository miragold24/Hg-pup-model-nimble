
###==========================================
###
### Input data & functions for HgNimble_Model_S4.R
###
###==========================================

library(nimble)
input_list = list(ndays = 123, amax=123,N=500,
                  start = 32, mu.bday = 45, sd.bday = 10,
                  mu.s4day = 17, sd.s4day = 5,
                  mu.mlen = 6, sd.mlen = 2,mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  skew.bday=5, x=1:123,
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.94, POMoultSD = 0.03,
                  PCWhite = 1,PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  obsdays = c(40, 50, 60,70,80,90))#c(35,45,55,65,75,85,95,105,115))#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)) #c(35, 60, 75,90,100,110))#30:123)##
thisenv <-environment()
list2env(input_list, thisenv)


# convert paramaters for skew normal
# delta <- skew.bday/sqrt(1+(skew.bday^2)) #location
# omega<-sqrt(sd.bday^2/(1-(2*delta^2)/pi)) #scaling factor
# xi <- mu.bday - omega*delta*sqrt(2/pi) #skewness
# # #x<- 1:123
####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================

whiteSum <- nimbleFunction(
  run = function(born=double(1),PS4_sum=double(1),Pm_sum=double(1)){
    # whitePups <- numeric(value = 0, length = 123)
    ndays=length(born)
    obsdays=  c(40, 50, 60,70,80,90)#c(35,45,55,65,75,85,95,105,115)#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110) #c(35, 60, 75,90,100,110)#30:123
    whitePups <- rep(0,ndays)
    
    for (d in 2:(ndays)){
      for (i in 1:d){
        age <- d-i
        whitePups[d] <- whitePups[d] + born[i] * (1-sum(PS4_sum[1:age])) *(1 - sum(Pm_sum[1:age]))
      }
    }
    returnType(double(1))
    return(whitePups)
  }
)
#test compile
#whiteSumc <- compileNimble(whiteSum)





S4Sum <- nimbleFunction(
  run = function(born=double(1),PS4_sum=double(1),Pm_sum=double(1)){
    ndays=length(born)
    obsdays=   c(40, 50, 60,70,80,90)#c(35,45,55,65,75,85,95,105,115)#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)#c(35, 60, 75,90,100,110)#30:123
    S4Pups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        S4Pups[d] <- S4Pups[d] + born[i]  * sum(PS4_sum[1:age])  * (1-sum(Pm_sum[1:age])) # * (1-sum(Pl_sum[1:age]))
      }
    }
    returnType(double(1))
    return(S4Pups)
  }
)
#test compile
#S4Sumc <- compileNimble(S4Sum)



moultSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1)){ #,moultObs=double(1)
    ndays=length(born)
    obsdays=  c(40, 50, 60,70,80,90)#c(35,45,55,65,75,85,95,105,115)#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)#c(35,75,90,110)#c(35, 60, 75,90,100,110)#30:123
    moultPups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        moultPups[d] <- moultPups[d] + born[i] * sum(Pm_sum[1:age])* (1-sum(Pl_sum[1:age])) 
      }
    }
    
    #Constrain moultPups to be greater than mObs
    # for (obs in 1:length(obsdays)){
    #   if (moultPups[obsdays[obs]] <= moultObs[obs]) {
    #     moultPups[obsdays[obs]] <- moultObs[obs] + 1  
    #   }
    # }

    returnType(double(1))
    return(moultPups)
  }
)
#test compile
#moultSumc <- compileNimble(moultSum)

#====================================================
###
### SKEW NORMAL functions
###
#====================================================
# PdfSkewNorm <- nimbleFunction(
#   run = function(x=double(1),skew.bday=integer(),mu.bday=integer(),sd.bday=integer()){
#   ndays = 123
#    # convert paramaters for skew normal
#   delta <- skew.bday/sqrt(1+(skew.bday^2)) #location
#   omega<-sqrt(sd.bday^2/(1-(2*delta^2)/pi)) #scaling factor
#   xi <- mu.bday - omega*delta*sqrt(2/pi) #skewness
#    (2/omega)*dnorm((x-xi)/omega) *
#     pnorm(skew.bday*((x-xi)/omega))
# 
#   }
# )

# dskewnorm <- function(x,skew.bday,mu.bday,sd.bday) {
#   Pb <- rep(0,ndays)  
#   #x=1:123
#   delta <- skew.bday/sqrt(1+(skew.bday^2))
#     
#     omega <- sqrt(sd.bday^2/(1-(2*delta^2)/pi))
#     
#     xi <- mu.bday - omega*delta*sqrt(2/pi)
#     Pb <- (2/omega)*dnorm((x-xi)/omega) * pnorm(skew.bday*((x-xi)/omega))
#     
#     # Return the probability vector
#     return(Pb)
#   }
# 
# 
# dskew<-nimbleRcall(function(x=double(0),skew.bday=integer(),mu.bday=integer(),sd.bday=integer()){},Rfun ="dskewnorm",
#                    returnType=double(1))
# 
# test<-dskew(x,skew.bday,mu.bday,sd.bday)
# 
# ### 1) integrate over days of season for prob of being born on each day
# integratePb <-  function(ndays,skew.bday,mu.bday,sd.bday) {
#     # Initialize Pb as a parameter of the same length as ndays
#     Pb <- rep(0, ndays)
#     
#     for (t in 2:ndays) {
#       #for (i in 1:int) {
#         # Define the start and end of the interval
#         start <- t-1   
#         end <-  t
#         
#         # Integrate over the interval and add to the probability for this day
#         Pb[t] <- Pb[t] + integrate(dskew,skew.bday=skew.bday,mu.bday=mu.bday,sd.bday=sd.bday, start, end)$value
#     }
#     
#     # Return the probability vector
#     return(Pb)
# }
# 
#  # ndays <- 123
#  #result <- integratePb(ndays,skew.bday,mu.bday,sd.bday)
# #plot(result)
# intPb<-nimbleRcall(function(ndays=double(0),skew.bday=integer(),mu.bday=integer(),sd.bday=integer()){},Rfun ="integratePb",
#                    returnType=double(1))
# ##compileNimble(intPb)




### 2) (EJ METHOD)
# # now need to make this a cdf and make sure it matches
# slice <-  function(ndays,nslice,skew.bday,mu.bday,sd.bday) {
#   ndays <- 123
#   nslice <- 1000
#   Pb <- rep(0, ndays)
#   for (i in 1:ndays){
#     for (j in 1:nslice){
#       dayslice <- i + (j/nslice)
#       Pb[i] <- Pb[i] + dskewnorm(dayslice, skew.bday, mu.bday, sd.bday)
#     }
#   }
#   return(Pb)
# }
# 
# 
# #Pb<-pborn/sum(pborn)
# slicePb<-nimbleRcall(function(ndays=double(0),nslice=double(0),skew.bday=integer(),mu.bday=integer(),sd.bday=integer()){},Rfun ="slice",
#                      returnType=double(1))
# #compileNimble(slicePb)


