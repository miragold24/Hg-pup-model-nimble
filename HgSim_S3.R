#==============================================================================
### 
### Simulation of grey seal pupswith only 3 processes (no S4 before moult): 
###  whitecoat, moult, leave 
###
#===============================================================================
source("leaveFunction.R")
simSealsS3 <- function(N = 500, 
   ndays = 123, amax = 123, mu.bday = 45, sd.bday = 10,
   birthcurve = "skewnormal",
   skew.bday = 5, 
   mu.mday = 23, sd.mday = 5,
   mu.lday = 31.5, sd.lday = 7,
   POWhite = 0.98, POWhiteSD = 0.01,
   POMoult = 0.95, POMoultSD = 0.03,
   PCWhite = 1, 
   PCMoult = 0.82, PCMoultSD = 0.09, 
   #PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
   #test all days as observed/surveyed, then reduce frequency of obsdays:
   obsdays = c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110) ,
   #                #c(35, 40, 50, 65, 70, 80, 95)
   #popproc: biological process can be generated as: 
   #         "detcont"=deterministic and continuous
   #         "stocdis"= stochastic and discrete
   popproc = "detcont", 
   # obsproc: observation process can be generated as: 
   #          "exact"= obs generated with no variation, 
   #          "normappx"= obs are normal approximated using binomial variance
   #          "binomial" = obs are binomially distributed 
   obsproc = "binom",seed=2024) {  
  
  #============================================================================
  ### BIOLOGICAL PROCESS
  ### Simulate births over season and transition to moulted and leaving colony: 
  ###  Pups born with white coats (stage 1-3) and begin to moult (stage 4) at  
  ###  ~ 17 days old. Pups complete the moult ~23 days and leave the colony for
  ###  sea ~31.5 days. Mortality is ~10% (98% of which occurs in stages 1-4)
  ###  
  ### Biological process error that may affect pup production estimates:
  ###      a) mortality not currently accounted for in model
  ###      b) 31.5 days likely too low for time to leave colony for sea
  ###      c) stage 4 not currently accounted for in model
  ###      *TO DO: -add mortality parameter
  ###              -add stage between whitecoat and moulted: 'moulting'(S4)
  ###              -sensitivity analysis: impact of # and timing of surveys
  #============================================================================
  library(lubridate)

  # initialize matrices

  set.seed(seed)
  
  if(birthcurve == "lognormal"){
    
    # convert birth parameters for lognormal dist
    lCV <- sd.bday/(mu.bday-start)
    lsd <- sqrt(log(lCV^2+1))
    lmu <- log(mu.bday-start)-(0.5*(lsd^2))
    
    # calculate probability of being born on each day
    Pb <- rep(0, ndays)
    
    for(t in 2:ndays){
      Pb[t] <- plnorm((t-start+1), lmu, lsd) - plnorm((t-start), lmu, lsd)
    }
    
  }
  
   if(birthcurve == "skewnormal"){
     # function to calculate prob of being born each day as CDF of skewnormal:
     
     dskewnorm <- function(x,skew.bday,mu.bday,sd.bday) {
       Pb <- rep(0,ndays)  
       #x=1:123
       delta <- skew.bday/sqrt(1+(skew.bday^2))
       
       omega <- sqrt(sd.bday^2/(1-(2*delta^2)/pi))
       
       xi <- mu.bday - omega*delta*sqrt(2/pi)
       Pb <- (2/omega)*dnorm((x-xi)/omega) * pnorm(skew.bday*((x-xi)/omega))
       
       # Return the probability vector
       return(Pb)
     }
     
     # Integrate over days of season for prob of being born on each day
     integratePb <-  function(ndays,skew.bday,mu.bday,sd.bday) {
       # Initialize Pb as a parameter of the same length as ndays
       Pb <- rep(0, ndays)
       
       for (t in 2:ndays) {
         # Define the start and end of the interval
         start <- t-1   
         end <-  t
         
         # Integrate over the interval and add to the probability for this day
         Pb[t] <- Pb[t] + integrate(dskewnorm,skew.bday=skew.bday,mu.bday=mu.bday,sd.bday=sd.bday, start, end)$value
       }
       
       return(Pb)
     }
     
     Pb <- integratePb(ndays,skew.bday,mu.bday,sd.bday)
     # make sure area under curve sums to 1
     Pb <- Pb/sum(Pb)
   } # end if

  ###
  ### Calculate probability of moulting as a function of age
  ###
  Pm <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax - 1), mu.mday, sd.mday)
  #Ensure probabilities sum to 1
  Pm <- Pm / sum(Pm)

  
  # *Use PmPl2.PlGivenm (prob of leaving given moulted) in leaveFunction.R
  Pl <- pnorm(1:amax, mu.lday, sd.lday) - pnorm(0:(amax - 1), mu.lday, sd.lday)
  #Ensure probabilities sum to 1
  Pl <- Pl / sum(Pl)


  ###
  ### simulate perfect birth, moult, leave process
  ###
  if (popproc == "detcont"){
    N <- 500
      realBorn <- Pb*N 
      
      whitePups <- rep(0, ndays)
      moultPups <- rep(0, ndays)

      for (d in 2:(ndays)){
        for (i in 1:d){
          age = d-i
          #calculate whitePups
          whitePups[d]<- whitePups[d] + realBorn[i] * (1 - sum(Pm[1:age])) 
          # Calculate moultPups
          moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) * 
            (1-sum(PmPl2$PlGivenm[1:age]))


        } # end for i

      } # end for d
      
    } # end if detcont

  
  
  ###
  ### simulate random birth, moult, leave process
  ###
  if (popproc == "stocdis"){
    N=500
    realBorn <- 0
    
    #The t function transposes this matrix, so each row of realBorn represents a draw.
    realBorn <- as.numeric(t(rmultinom(1, size = N, prob = Pb)))

    
    realMoult <- matrix(rep(0, ndays * amax), nrow = ndays, ncol = amax)
    cumMoult <- numeric(ndays)# stores cumulative sum of moulted pups
    for (i in 2:ndays) { 
      for (j in 1:amax) {
        realMoult[i, ] <- as.numeric(t(rmultinom(1, size = realBorn[i], prob = Pm)))
        
        
        # Extract moulting data for current and previous days
        M <- realMoult[1:i, 1:j]
      }
      # Flip the upper triangle to get daily and cumulative moulting
      flip.M <- M[i:1, ]
      cumMoult[i] <- sum(M) - sum(flip.M[upper.tri(flip.M)]) #sum of fliped M is whitePups present each day
    }
    
    
    #
    # calculate pups that leave conditional on having moulted
    #
    realLeave = array(data = 0, dim = c(ndays,amax,amax))
    for (dob in 1:nrow(realMoult)) {
      for (age in 1:ncol(realMoult)) {
        realLeave[dob,age,] <- as.vector(rmultinom(1, size = realMoult[dob,age], prob = Pl))
      }
    }
    
    realLeave2d = apply(realLeave, c(1,3), sum)
    dayLeave = matrix(0, nrow = ndays, ncol = amax)
    
    for(dob in 1:nrow(realLeave2d)){
      for(dol in 1:ncol(realLeave2d)){
        day_leave = min(dol + dob - 1,123)
        dayLeave[dob,day_leave] = dayLeave[dob,day_leave] + realLeave2d[dob,dol]
      }
    }
  
    cumLeave<-cumsum(apply(dayLeave, 2, sum))

    whitePups<- cumsum(realBorn) - cumMoult
    moultPups <- cumMoult - cumLeave

}

  #============================================================================
  ### OBSERVATION PROCESS
  ### Simulate observations of whitecoat and moulted pups
  ### account for imperfect detection and imperfect classification error
  ### Data:this simulation is based on data collected from aerial surveys. From
  ###      2012, digital images stitched together to display an entire colony 
  ###      have been used to count pups. Sources of obs error that may affect
  ###      pup production estimates:
  ###      a) the possibility of missing pups during counting (POWhite/POMoult)
  ###      b) misclassifying pups (whitecoat and moulted) (PCWhite/PCMoult)
  ###       -currently whitecoats cannot be misclassified as moulted (PCWhite=1)
  ###      *TO DO: add stage between whitecoat and moulted: 'moulting'(S4)
  ###              
  #============================================================================
  
  
   # if (obsproc == "exact"){
   # 
   #  whiteObs <- whitePups*POWhite*PCWhite + moultPups*POMoult*(1-PCMoult)
   #  moultObs <- moultPups*POMoult*PCMoult
   # }
     
  if (obsproc == "exact"){
    whiteCount <- whitePups*POWhite
    moultCount <- moultPups*POMoult

    # white classified as white
    WW <- whiteCount*PCWhite
    # #whitecoats classified as moulted
    WM <- whiteCount-WW
    # # moulted classified as moulted
    MM <- moultCount*PCMoult
    # moulted classified as white
    MW <- moultCount*(1-PCMoult)

   }

  if (obsproc == "normappx"){
    sdwhite <- sqrt(whitePups * POWhite * (1 - POWhite) + 
                      moultPups * (1 - PCMoult) * POMoult * 
                      (1 - (1 - PCMoult) * POMoult))
    sdmoult <- sqrt(moultPups * POMoult * PCMoult * (1 - POMoult * PCMoult))

    whiteCount <- rnorm(n = length(whitePups), mean=whitePups * POWhite + 
                          moultPups * POMoult * (1 - PCMoult), sd=sdwhite)
    moultCount <- rnorm(n = length(moultPups), mean= moultPups * POMoult * 
                          PCMoult, sd=sdmoult)

   
    # white classified as white:
    WW <- round(whiteCount)
    # No whitecoats misclassified as moulted
    WM <- 0
    # moulted classified as moulted
    MM <- round(moultCount)
    # This is already acounted for in moultCount calculation (**clunky!**)
    MW <- 0

     }

  # if (obsproc == "binom"){
  #   whiteCount <- rbinom(n = length(whitePups), size = whitePups, POWhite)
  #   moultCount <- rbinom(n = length(moultPups), size = moultPups, POMoult)
  #   # add in imperfect classification
  #   # white classified as white
  #   WW <- rbinom(n = length(whiteCount), size = whiteCount, PCWhite)
  #   # white classified as moulted
  #   WM <- whiteCount - WW
  #   # moulted classified as moulted
  #   MM <- rbinom(n = length(moultCount), size = moultCount, PCMoult)
  #   # moulted classified as white:
  #   # - correctly observed moulteds - correcty classified moulteds
  #   MW <- moultCount - MM  
  # 
  #  }
  
  if (obsproc == "binom"){
  whiteCount <- rbinom(n = length(whitePups), size = round(whitePups), POWhite)
  moultCount <- rbinom(n = length(moultPups), size = round(moultPups), POMoult)
  WW <- whiteCount
  # No whitecoats misclassified as moulted
  WM <- 0
  # moulted classified as moulted
  MM <- moultCount
  # This is already acounted for in moultCount calculation (**clunky!**)
  MW <- 0
  
  }
  #simulated observation data of white and moulted:
  whiteObs <- WW + MW
  moultObs <- MM + WM

  wObs <- whiteObs[obsdays]
  mObs <- moultObs[obsdays]

  #==========================================================================


  sim.out <- list(realBorn = realBorn, Pb = Pb,Pm = Pm, Pl=PmPl2$PlGivenm,
                   wObs = wObs, mObs = mObs, whiteObs=whiteObs,moultObs=moultObs,
                 whitePups=whitePups,moultPups=moultPups)
  list2env(sim.out, globalenv())
  return(sim.out)

}

sim_S3<-simSealsS3()


