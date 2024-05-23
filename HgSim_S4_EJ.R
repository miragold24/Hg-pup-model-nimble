### ===================
### HG pup simulation with 4 stages: whitecoat, moulting (S4) moulted, leave 
### ===================
source("leaveFunction.R")
input_list = list(ndays = 123, amax=123,N=500, mu.bday = 45, sd.bday = 10,
                  mu.s4day = 17, sd.s4day = 5,
                  mu.mlen = 6, sd.mlen = 2,mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  skew.bday=5, 
                  POWhite = 0.98, POWhiteSD = 0.01,
                  POMoult = 0.94, POMoultSD = 0.03,
                  PCWhite = 1,PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
                  PCMoult = 0.82, PCMoultSD = 0.09,
                  obsdays = 1:123) #c(35,45,55,65,75,85,95,105,115))#c(50, 60,70,80,90))#c(35, 50, 75,90,105,120))#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110)) #
thisenv <-environment()
list2env(input_list, thisenv)


simSealsS4 <- function(N = 500, 
   ndays = 123, amax = 123,
   start = 32, mu.bday = 45, sd.bday = 10, 
   birthcurve = "skewnormal",
   skew.bday = 5,
   mu.mday = 23, sd.mday = 5,
   mu.lday = 31.5, sd.lday = 7,
   mu.s4day=17,sd.s4day=5,
   POWhite = 0.98, POWhiteSD = 0.01,
   POMoult = 0.95, POMoultSD = 0.03,
   PCWhite = 1, 
   PCMoult = 0.82, PCMoultSD = 0.09, 
   PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
   obsdays =  1:123, #c(35,45,55,65,75,85,95,105,115),#c(50, 60,70,80,90),#c(35, 50, 75,90,105,120),#30:123, ,#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110),#
   popproc = "stocdis", #detcont = deterministic and continuous, option for 
   # "stocdis" - stochastic and discrete
   obsproc = 'binom',seed=2024) {  

   library(lubridate)

  #########
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
  } #end if

  #===========================
  if(birthcurve == "skewnormal"){
    
    #below checks if simulated matching model can return values
    # skew.bday=0.3335911
    # mu.bday=Cmodel$mu.bday
    # sd.bday=Cmodel$sd.bday
    dskewnorm <- function(x,skew.bday,mu.bday,sd.bday) {
      Pb <- rep(0,ndays)  
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
  

  #=============================================================================
  ###
  ### TRANSITION PROBABILITIES
  ###
  ### calculate the probability of transitioning into moulting (PS4) or complete 
  ### moult (Pm) on any given day by subtracting the CDF on current day from 
  ### previous day to get PS4 and Pm for current day: 
  #=============================================================================
  PS4 <- rep(0, amax)
  
  ###  PS4 as function of age
  PS4<- pnorm(1:amax, mu.s4day, sd.s4day) - pnorm(0:(amax - 1), mu.s4day, sd.s4day)
  
  PS4 <- PS4/sum(PS4)
  
  
  Pm<- rep(0, amax)
  # #Calculate probability of moulting as a function of age
  Pm <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax - 1), mu.mday, sd.mday)
  #Ensure probabilities sum to 1
  Pm <- Pm / sum(Pm)

  #Calculate probability of leaving as a function of age
  Pl <- pnorm(1:amax, mu.lday, sd.lday) - pnorm(0:(amax - 1), mu.lday, sd.lday)
  #Ensure probabilities sum to 1
  Pl <- Pl / sum(Pl)

  # simulate perfect birth, moult, leave process
  if (popproc == "detcont"){
    #N <- 500
      realBorn <- Pb*N 
      
      whitePups <- rep(0, ndays)
      S4Pups <- rep(0, ndays)
      moultPups <- rep(0, ndays)
      #leavePups <- rep(0, ndays)
      for (d in 2:(ndays)){
        for (i in 1:d){
          age = d-i
          whitePups[d]<- whitePups[d] + realBorn[i] * (1-sum(PS4[1:age])) *(1 - sum(Pm[1:age])) 
          S4Pups[d] <- S4Pups[d] + realBorn[i]  * sum(PS4[1:age])  * (1-sum(Pm[1:age])) 
          moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) * (1-sum(PmPl2$PlGivenm[1:age]))
         # leavePups[d] <- leavePups[d] + moultPups[d] * sum(PmPl2$PlGivenm[1:age])
        } # end for i

      } # end for d
      
    } # end if detcont

  
  
  # simulate randomized birth, moult, leave process
  
  
  if (popproc == "stocdis"){
    N=500
    realBorn <- 0
    
    #The t function transposes this matrix, so each row of realBorn represents a draw.
    realBorn <- as.numeric(t(rmultinom(1, size = N, prob = Pb)))
    realS4 <- matrix(rep(0, amax*ndays), nrow = ndays)
    realMoult <- matrix(rep(0, amax*ndays), nrow = ndays)
    realLeave <- matrix(rep(0, amax*ndays), nrow = ndays)
    
    for (i in 1:ndays){
      
      while(sum(realS4[i,]) != sum(realBorn[i])){
        realS4[i,] <- rbinom(n = amax, size = realBorn[i], prob = PS4)
      }
      while(sum(realMoult[i,]) != sum(realS4[i])){
        realMoult[i, ] <- rbinom(n = amax, size = realS4[i], prob = Pm)
      }
      while(sum(realLeave[i, ]) != sum(realMoult[i])){
        realLeave[i, ] <- rbinom(n = amax, size = realMoult[i], prob = Pl)
      }
    } # end for i
    
    
    # now count number of individuals per day
    
    # initialize vectors
    
    tmin <- min(which(realBorn!=0))
    
    s4vec <- rep(0, ndays)
    moultvec <- rep(0, ndays) # number that have moulted by each day
    leavevec <- rep(0, ndays) # number that have left by each day
    
    # this takes the upper left triangular part of the moult and leave matrices
    
    for (i in tmin:ndays){
      
      s4v <- 0
      mv <- 0
      lv <- 0
      
      rowmax <- i - tmin + 1
      
      for (j in tmin:i){
        
        for (k in 1:rowmax){
          s4v <- sum(s4v, realS4[j, k])
          mv <- sum(mv, realMoult[j, k])
          lv <- sum(lv, realLeave[j, k])
          
        } # end for k
        
        rowmax <- rowmax - 1
      }
      s4vec[i] <- s4v
      moultvec[i] <- mv
      leavevec[i] <- lv
      
    }
    
    # determine who's there on each day
    
    whitePups <- cumsum(realBorn) - s4vec
    S4Pups <- s4vec-moultvec
    moultPups <- moultvec - leavevec
    
    moultPups[which(moultPups<0)] <- 0
    allPups <- cumsum(realBorn) - leavevec
    
    if (any(allPups != (whitePups + S4Pups + moultPups))){
      warning("Pups don't add up. This may be because some pups left before moulting.")
    }
    
  } # end if
  
  
  if (obsproc == "exact"){
    whiteObs <- whitePups*POWhite + moultPups*POMoult*(1-PCMoult)+
                S4Pups*POWhite*PCWhiteS4
    moultObs <- moultPups*POMoult*PCMoult +S4Pups*POWhite*(1-PCWhiteS4)

  WW <- whiteObs
  MW <-0
  MM <- moultObs
  WM <-0
  whiteCount <- WW
  moultCount <-  MM
  
  }

  if (obsproc == "pois"){
     whiteM<-rep(0,ndays)
     moultM<-rep(0,ndays)
     whiteCount<-rep(0,ndays)
     moultCount<-rep(0,ndays)
    for (i in 1:length(obsdays)){
      whiteM[i] <- whitePups[obsdays[i]] * POWhite + S4Pups[obsdays[i]] * POWhite*PCWhiteS4 + moultPups[obsdays[i]] * POMoult*(1-PCMoult)
      moultM[i] <- S4Pups[obsdays[i]] * POWhite*(1-PCWhiteS4) + moultPups[obsdays[i]] * POMoult*PCMoult

      whiteCount[i] <- rpois(n = length(whiteM[obsdays[i]]), whiteM[i])
      moultCount[i] <- rpois(n = length(moultM[obsdays[i]]), moultM[i])
    }
    
  }
  
  
  if (obsproc == "normappx"){
    whiteM <-rep(0,ndays)
    moultM <- rep(0,ndays)
    whiteV<-rep(0,ndays)
    moultV<-rep(0,ndays)
    whiteS<-rep(0,ndays)
    moultS<-rep(0,ndays)
    PCStge4=PCWhiteS4
    for(i in 1:length(obsdays)) {

      moultS[i] <- max(sqrt(moultV[i]), 0.00000001)
      whiteM[i] <- whitePups[obsdays[i]] * POWhite               + S4Pups[obsdays[i]] * POWhite*PCStge4                       + moultPups[obsdays[i]] * POMoult*(1-PCMoult)
      whiteV[i] <- whitePups[obsdays[i]] * POWhite * (1-POWhite) + S4Pups[obsdays[i]] * POWhite*PCStge4 * (1-POWhite*PCStge4) + moultPups[obsdays[i]] * POMoult*(1-PCMoult) * (1-POMoult*(1-PCMoult))
      whiteS[i] <- max(sqrt(whiteV[i]), 0.00000001)
      moultM[i] <- S4Pups[obsdays[i]] * POWhite*(1-PCStge4)                         + moultPups[obsdays[i]] * POMoult*PCMoult
      moultV[i] <- S4Pups[obsdays[i]] * POWhite*(1-PCStge4)*(1-POWhite*(1-PCStge4)) + moultPups[obsdays[i]] * POMoult*PCMoult * (1-POMoult*PCMoult)
      moultS[i] <- max(sqrt(moultV[i]), 0.00000001)
      
      whiteCount[i] <- rnorm(n = length(whiteM[obsdays[i]]), mean=whiteM[obsdays[i]], sd=whiteS[i])
      moultCount[i] <- rnorm(n = length(moultM[obsdays[i]]), mean= moultM[obsdays[i]],sd=moultS[i])
      
    }

     ## CLASSIFICATION -> POBS
    
    # whiteV <- whitePups * POWhite * (1-POWhite) + 
    #   S4Pups * POWhite*PCStge4 * (1-POWhite*PCStge4) + 
    #   moultPups* POMoult*(1-PCMoult) * (1-POMoult*(1-PCMoult))
    # whiteS <- max(sqrt(whiteV), 0.00000001)
    # 
    # moultV <- S4Pups * POWhite*(1-PCStge4)*(1-POWhite*(1-PCStge4)) + 
    #   moultPups * POMoult*PCMoult * (1-POMoult*PCMoult)
    # 
    # whiteM <- whitePups * POWhite  + S4Pups * POWhite*PCStge4   + moultPups[obsdays[i]] * POMoult*(1-PCMoult)
    # whiteCount <- rnorm(n = length(whitePups), mean=whitePups * POWhite +
    #                       moultPups * (1 - PCMoult)+
    #                       S4Pups*POWhite*PCWhiteS4, sd=whiteS)
    # 
    # whiteCount<- whiteCount*POWhite
    # 
    # moultCount <- rnorm(n = length(moultPups), mean= moultPups * PCMoult +
    #                       S4Pups*(1-PCWhiteS4), sd=sdmoult)
    # 
    # moultCount <- moultCount*POMoult


  }

   if (obsproc == "binom"){

     #White detected as white stage 1-3 and S4
     WW1_3 <- rbinom(n = length(whitePups), size = round(whitePups), 1)
     WWS4  <- rbinom(n = length(S4Pups), size = round(S4Pups), PCWhiteS4)
     #All whitecoats correctly classified as whitecoats:
     WW <- WW1_3 + WWS4
     #moulteds incorrectly classified as whitecoats
     MW <- rbinom(n = length(moultPups), size = round(moultPups), 1-PCMoult)
     # Add detection probability (which is a consequence of how an observer 
     # classifies the animal)
     whiteCount <-  rbinom(n = length(WW+MW), size = round(WW+MW), POWhite)
     
     #Moulted classified as moulted
     MM <- rbinom(n = length(moultPups), size = round(moultPups), PCMoult)
     # S4 incorrectly detected as moulted
     WMS4 <- rbinom(n = length(S4Pups), size = round(S4Pups), 1-PCWhiteS4)
     
     #add detection probability for all pups classified as moulted:
     moultCount<-rbinom(n = length(MM+WMS4), size = round(MM+WMS4), POMoult)
    

   }

  # simulated observation data of white and moulted:
  whiteObs <- whiteCount
  moultObs <- moultCount
  # whiteObs <- WW + MW
  # moultObs <- MM + WM
  
  
  wObs <- whiteObs[obsdays]
  mObs <- moultObs[obsdays]
 
  sim.out <- list( realBorn = realBorn, Pb = Pb,PS4=PS4,Pm = Pm,
                  whiteObs = whiteObs, moultObs = moultObs,whiteCount=whiteCount,moultCount=moultCount,
                  wObs = wObs, mObs = mObs, white1_4=whitePups+S4Pups,
                 whitePups=whitePups,S4Pups=S4Pups,moultPups=moultPups)
  list2env(sim.out, globalenv())
  return(sim.out)

}

sim_S4<-simSealsS4()

####################################
# visualize actual pups on colony and counted pups: 

# Define colorblind-friendly palette
# cbPalette <- c("#E69F00", "#cc8400", "#56B4E9","#3c7da3")
# # Plot the lines for whitePups and S4Pups
# plot(whitePups + S4Pups, typ = 'l', col = cbPalette[1], ylab = '# pups', xlab = 'day of season', main = 'Present vs Observed Pups')
# lines(moultPups, col = cbPalette[3])
# 
# # Add points for whiteCount and moultCount on observation days
# #obsdays <- c(35, 50, 75, 90, 105, 120)
# points(x = obsdays, y = whiteCount[obsdays], col = "#cc8400", pch = 19)
# points(x = obsdays, y = moultCount[obsdays], col = "#3c7da3", pch = 19)
# 
# # Create legend with points
# legend("topright", legend = c("White present", "White observed", "Moult present", "Moult observed"), col = cbPalette, lty = c(1,NA,1,NA), pch = c(NA, 19, NA, 19))
# 
