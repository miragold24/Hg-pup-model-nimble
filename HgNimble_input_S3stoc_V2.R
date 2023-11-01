
###==========================================
###
### Input data & functions for HgNimble_StocModel_V3.R
###
###==========================================
library(nimble)
input_list = list(ndays = 123, amax=123,
                  start = 32, mu.bday = 45, sd.bday = 10,
                  mu.s4day = 17, sd.s4day = 5,
                  mu.mlen = 6, sd.mlen = 2,mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.94, POMoultSD = 0.03,
                  PCWhite = 1, 
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,seed=2024,
                  obsdays =  c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)
                 )
thisenv <-environment()
list2env(input_list, thisenv)
####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================

simSealsS3_W <- function(realBorn,Nborn,Pb_sum,Pm_sum,Pl_sum){  
  
  #########
  # initialize matrices
  ndays=123
  amax=123
  Nborn=500
  set.seed(2024)
  
  realBorn <- 0
  while(sum(realBorn) != Nborn){
    realBorn <- rbinom(n = ndays, size = Nborn, prob = Pb_sum)
  }
  
  leavers <- matrix(rep (0, amax*amax), nrow = amax)
  realMoult <- matrix(rep(0, ndays*amax), nrow = ndays, ncol = amax)
  realLeave <- matrix(nrow = ndays, ncol = amax)
  
  for (i in 1:ndays){
    
    while(sum(realMoult[i,]) != sum(realBorn[i])){
      realMoult[i, ] <- rbinom(n = amax, size = realBorn[i], prob = Pm_sum)
    }
    
    for (j in 1:amax){
      while(sum(leavers[j,]) != realMoult[i, j]){
        leavers[j, ] <- c(rep(0, j), rbinom(size = realMoult[i, j], prob = Pl_sum, n = amax)[j:(amax-1)])
      } # end while
    } # end for j
    
    realLeave[i,] <- colSums(leavers)
    
  } # end for i
  
  # now count number of individuals per day
  
  # initialize vectors
  
  tmin <- min(which(realBorn!=0))
  
  moultvec <- rep(0, ndays) # number that have moulted by each day
  leavevec <- rep(0, ndays) # number that have left by each day
  
  for (i in tmin:ndays){
    
    mv <- 0
    lv <- 0
    
    rowmax <- i - tmin + 1
    
    for (j in tmin:i){
      
      for (k in 1:rowmax){
        
        mv <- sum(mv, realMoult[j, k])
        lv <- sum(lv, realLeave[j, k])
        
      } # end for k
      
      rowmax <- rowmax - 1
    }
    
    moultvec[i] <- mv
    leavevec[i] <- lv
    
  } # end for i
  
  whitePups <- cumsum(realBorn) - moultvec
  moultPups <- moultvec - leavevec
  
  # Pw=POWhite*PCWhite*POMoult*(1-PCMoult)
  # Pm=POMoult*PCMoult
  
  #list2env(sim.out, globalenv())
  return(whitePups)
  
}

NimSimS3_W <- nimbleRcall(
  prototype=function(realBorn=double(1),Nborn=double(0),
                     Pb_sum=double(1),Pm_sum=double(1),Pl_sum=double(1)
  ){},
  
  returnType=double(1),
  Rfun = 'simSealsS3_W'
  
)

#nimsimS3W<-NimSimS3_W(realBorn=sim_S3$realBorn,Nborn=500,
                      #Pb_sum=sim_S3$Pb,Pm_sum=sim_S3$Pm,Pl_sum=sim_S3$Pl)




simSealsS3_M <- function(realBorn,Nborn,Pb_sum,Pm_sum,Pl_sum){  
  
  #########
  # initialize matrices
  ndays=123
  amax=123
  Nborn=500
  set.seed(2024)

  realBorn <- 0
  while(sum(realBorn) != Nborn){
    realBorn <- rbinom(n = ndays, size = Nborn, prob = Pb_sum)
  }
  
  leavers <- matrix(rep (0, amax*amax), nrow = amax)
  realMoult <- matrix(rep(0, ndays*amax), nrow = ndays, ncol = amax)
  realLeave <- matrix(nrow = ndays, ncol = amax)
  
  for (i in 1:ndays){
    
    while(sum(realMoult[i,]) != sum(realBorn[i])){
      realMoult[i, ] <- rbinom(n = amax, size = realBorn[i], prob = Pm_sum)
    }
    
    for (j in 1:amax){
      while(sum(leavers[j,]) != realMoult[i, j]){
        leavers[j, ] <- c(rep(0, j), rbinom(size = realMoult[i, j], prob = Pl_sum, n = amax)[j:(amax-1)])
      } # end while
    } # end for j
    
    realLeave[i,] <- colSums(leavers)
    
  } # end for i
  
  # now count number of individuals per day
  
  # initialize vectors
  
  tmin <- min(which(realBorn!=0))
  
  moultvec <- rep(0, ndays) # number that have moulted by each day
  leavevec <- rep(0, ndays) # number that have left by each day
  
  for (i in tmin:ndays){
    
    mv <- 0
    lv <- 0
    
    rowmax <- i - tmin + 1
    
    for (j in tmin:i){
      
      for (k in 1:rowmax){
        
        mv <- sum(mv, realMoult[j, k])
        lv <- sum(lv, realLeave[j, k])
        
      } # end for k
      
      rowmax <- rowmax - 1
    }
    
    moultvec[i] <- mv
    leavevec[i] <- lv
    
  } # end for i
  
  whitePups <- cumsum(realBorn) - moultvec
  moultPups <- moultvec - leavevec
  
  # Pw=POWhite*PCWhite*POMoult*(1-PCMoult)
  # Pm=POMoult*PCMoult
  
  #list2env(sim.out, globalenv())
  return(moultPups)
  
}

NimSimS3_M <- nimbleRcall(
  prototype=function(realBorn=double(1),Nborn=double(0),
                     Pb_sum=double(1),Pm_sum=double(1),Pl_sum=double(1)
  ){},
  
  returnType=double(1),
  Rfun = 'simSealsS3_M'
  
)


# PmfSumBin<- function(dw, dm, pcw, pcm) {
#   # Calculate pmf for all possible values of yw
#   # Calculate the second part, and complete the equation, for each
#   # value of yw from 0 to (dw+wm)
#   pyw <- numeric(dw + dm + 1)
#   for(yw in 0:(dw + dm)){
#     i <- 0:(min(dw, yw))
#     pdww <- dbinom(i, dw, pcw)
#     pdmw <- dbinom(yw - i, dm, (1 - pcm))
#     pyw[yw+1] <- sum(pdww * pdmw)
#   }
#   # Check it sums to 1!
#   sum(pyw)
#   return(pyw)
# }
# dw <- 50
# dm <- 90
# pcw <- 0.9
# pcm <- 0.8
# 
# sum_binom<-PmfSumBin(dw,dm,pcw,pcm)
# sum(sum_binom)

# 
# PmfSumBin <- nimbleFunction(
#   setup = function(dw = integer(), dm = integer(), pcw = double(), pcm = double()) {
#     # Define the parameters
#     dw <- dw
#     dm <- dm
#     tpcw <- pcw
#     pcm <- pcm
#   },
#   run = function() {
#     # Calculate pmf for all possible values of yw
#     # Calculate the second part, and complete the equation, for each
#     # value of yw from 0 to (dw+dm)
#     pyw <- numeric(dw + dm + 1)
#     for(yw in 0:(dw + tdm)){
#       i <- 0:(min(dw, yw))
#       pdww <- dbinom(i, dw, this$pcw)
#       pdmw <- dbinom(yw - i, dm, (1 - pcm))
#       pyw[yw+1] <- sum(pdww * pdmw)
#     }
#     # Check it sums to 1!
#     returnType(double())
#     return(pyw)
#   }
# )
# Bin<- NimbleComPmfSumBin()
