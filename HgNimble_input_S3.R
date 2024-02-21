library(nimble)
###==========================================
###
### Input data & functions for HgNimble_Model_V3.R
### Steps:
### 1) set obsdays
### 2) source HgSim_S3.R
### 3) run HgNimble_Model_S3.R
### to calculate pups present in each class on a given day

###==========================================
input_list = list(ndays = 123, amax=123,N=500,
                  start = 32, mu.bday = 45, sd.bday = 10,
                  skew.bday=5,
                  mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  mu.s4day=17,sd.s4day=4,
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.95, POMoultSD = 0.03,
                  PCWhite = 1, PCWhiteS4=0.75,
                  PCWhiteS4SD=0.01,
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  obsdays =  c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110), #c(35, 40, 50, 65, 70, 80, 95),
                  seed=2024)
thisenv <-environment()
list2env(input_list, thisenv)

####==============================================================
###
### Define NIMBLE functions 
###
###===============================================================

SkewNorm<- nimbleFunction(
  run=function(Pb_sum=double(1)){
    Pb <- rep(0, ndays)
    
    for(t in 2:ndays){
      Pb[t] <- pskew_normal(q = t, mu = mu.bday, sigma = sd.bday, alpha = skew.bday) -
        pskew_normal(q = t-1, mu = mu.bday, sigma = sd.bday, alpha = skew.bday)
    }
    
    Pb_sum <- Pb/sum(Pb)
    
  }
)




whiteSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1)){
    ndays=length(born)
    obsdays=c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110) #c(35, 40, 50, 65, 70, 80, 95)#
    whitePups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age <- d-i
        whitePups[d] <- whitePups[d] + born[i] *(1 - sum(Pm_sum[1:age])) # * (1 - sum(Pl_sum[1:age]))
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
    


moultSum <- nimbleFunction(
  run = function(born=double(1),Pm_sum=double(1),Pl_sum=double(1)){
    ndays=length(born)
    obsdays= c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110) #c(35, 40, 50, 65, 70, 80, 95)#
    moultPups <- rep(0,ndays)
    for (d in 2:(ndays)){
      for (i in 1:d){
        age = d-i
        moultPups[d] <- moultPups[d] + born[i] * sum(Pm_sum[1:age]) * (1-sum(Pl_sum[1:age]))
      }
    }
    
    #Constrain moultPups to be greater than mObs
    # for (obs in 1:length(obsdays)){
    #   if (moultPups[obsdays[obs]] <= moultObs[obs]) {
    #     moultPups[obsdays[obs]] <- moultObs[obs] + 1  # Adjust moultPups to meet the constraint
    #   }
    # }
    # 
    returnType(double(1))
    return(moultPups)
  }
)
 #moultSumc <-compileNimble(moultSum)

dskewnorm <- function(x,skew.bday,mu.bday,sd.bday) {
  delta <- skew.bday/sqrt(1+(skew.bday^2))
  
  omega <- sqrt(sd.bday^2/(1-(2*delta^2)/pi))
  
  xi <- mu.bday - omega*delta*sqrt(2/pi)
  (2/omega)*dnorm((x-xi)/omega) * pnorm(skew.bday*((x-xi)/omega))
}


integratePb <-  function(ndays,int,skew.bday,mu.bday,sd.bday) {
  # Initialize Pb as a parameter of the same length as ndays
  Pb <- numeric(ndays)
  
  for (day in 1:ndays) {
    # Initialize the probability for this day
    Pb_day <- 0
    
    for (i in 1:int) {
      # Define the start and end of the interval
      start <- (day - 1) + (i - 1) / int
      end <- start + 1 / int
      
      # Integrate over the interval and add to the probability for this day
      Pb_day <- Pb_day + integrate(dskewnorm,skew.bday=skew.bday,mu.bday=mu.bday,sd.bday=sd.bday, start, end)$value
    }
    
    # Store the probability for this day
    Pb[day] <- Pb_day
  }
  
  # Return the probability vector
  return(Pb)
}

intPb<-nimbleRcall(function(ndays=double(0),int=double(0),skew.bday=integer(),mu.bday=integer(),sd.bday=integer()){},Rfun ="integratePb",
                   returnType=double(1))