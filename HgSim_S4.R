### ===================
### HG sim with only 3 stages: whitecoat, moult, leave 
### ===================
source("leaveFunction.R")
simSealsS4 <- function(N = 500, 
   ndays = 123, amax = 123, 
   start = 32, mu.bday = 45, sd.bday = 10, 
   #mu.s4day = 17, sd.s4day = 5,
   mu.mlen = 6, sd.mlen = 2,
   mu.mday = 23, sd.mday = 5,
   mu.lday = 31.5, sd.lday = 7,
   mu.s4day=17,sd.s4day=5,
   POWhite = 0.98, POWhiteSD = 0.01,
   POMoult = 0.95, POMoultSD = 0.03,
   PCWhite = 1, 
   PCMoult = 0.82, PCMoultSD = 0.09, 
   PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
   obsdays = 33:123,#c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110),
   popproc = "detcont", #detcont = deterministic and continuous, option for "stocdis" - stochastic and discrete
   obsproc = "exact",seed=2024) {  

   library(lubridate)

  #########
  # initialize matrices

  set.seed(seed)

  # convert birth parameters for lognormal dist
  lCV <- sd.bday/(mu.bday-start)
  lsd <- sqrt(log(lCV^2+1))
  lmu <- log(mu.bday-start)-(0.5*(lsd^2))
  
  # calculate probability of being born on each day
  Pb <- rep(0, ndays)
  
  for(t in 2:ndays){
    Pb[t] <- plnorm((t-start+1), lmu, lsd) - plnorm((t-start), lmu, lsd)
  }
  Pb <- Pb / sum(Pb)
  
  PS4 <- rep(0, amax)
  
  ###  PS4 as function of age
  PS4<- pnorm(1:amax, mu.s4day, sd.s4day) - pnorm(0:(amax - 1), mu.s4day, sd.s4day)
  
  PS4 <- PS4/sum(PS4)
  
  # #Calculate probability of moulting as a function of age
  Pm <- pnorm(1:amax, mu.mday, sd.mday) - pnorm(0:(amax - 1), mu.mday, sd.mday)
  #Ensure probabilities sum to 1
  Pm <- Pm / sum(Pm)

  # #Calculate probability of leaving as a function of age
  # Pl <- pnorm(1:amax, mu.lday, sd.lday) - pnorm(0:(amax - 1), mu.lday, sd.lday)
  # #Ensure probabilities sum to 1
  # Pl <- Pl / sum(Pl)

  # simulate perfect birth, moult, leave process
  if (popproc == "detcont"){
    N <- 500
      realBorn <- Pb*N 
      
      whitePups <- rep(0, ndays)
      S4Pups <- rep(0, ndays)
      moultPups <- rep(0, ndays)

     # wleavers <- rep(0,ndays)
      #mleavers <- rep(0,ndays)
     # #totPups <- rep(0, ndays)
      for (d in 2:(ndays)){
        for (i in 1:d){
          age = d-i
          #calculate whitePups where no whitePups can leave
          whitePups[d]<- whitePups[d] + realBorn[i] * (1-sum(PS4[1:age])) *(1 - sum(Pm[1:age])) #* (1-sum(Pl[1:age]))
          S4Pups[d] <- S4Pups[d] + realBorn[i]  * sum(PS4[1:age])  * (1-sum(Pm[1:age])) 
          #moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) *(1-sum(Pl[1:age]))
          moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) * (1-sum(PmPl2$PlGivenm[1:age]))


        } # end for i

      } # end for d
      
    } # end if detcont

  
  
  # simulate randomized birth, moult, leave process
  if (popproc == "stocdis"){
    set.seed(seed)

    realBorn <- 0
    while(sum(realBorn) != N){
      realBorn <- rbinom(n = ndays, size = N, prob = Pb)
    }
    
    leavers <- matrix(rep (0, amax*amax), nrow = amax)
    realMoult <- matrix(rep(0, ndays*amax), nrow = ndays, ncol = amax)
    realLeave <- matrix(nrow = ndays, ncol = amax)
    
    for (i in 1:ndays){
      while(sum(realMoult[i,]) != sum(realBorn[i])){
        realMoult[i, ] <- rbinom(n = amax, size = realBorn[i], prob = Pm)
      }

      for (j in 1:amax){
        while(sum(leavers[j,]) != realMoult[i, j]){
          leavers[j, ] <- c(rep(0, j), rbinom(size = realMoult[i, j], prob = Pl, n = amax)[j:(amax-1)])
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
    allPups <- cumsum(realBorn) - leavevec
    
  } # end if stocdisc 
  
  if (obsproc == "exact"){
    
    whiteCount <- whitePups*POWhite 
    S4Count <- S4Pups*POWhite
    moultCount <- moultPups*POMoult
    
    
    # add in imperfect classification
    
    # S4 classified as white
    S4W <- S4Count*PCWhiteS4
    # white classified as white
    WW <- whiteCount + S4W
    # white classified as moulted
    WM <- S4Count - S4W
    # moulted classified as moulted
    MM <- moultCount*PCMoult
    # moulted classified as white
    MW <- moultCount - MM
    

  }
  
  
  if (obsproc == "normappx"){
    set.seed(seed)
    
    sdwhite <- sqrt(whitePups * POWhite * (1 - POWhite) + moultPups * (1 - PCMoult) * POMoult * (1 - (1 - PCMoult) * POMoult)) +
      S4Pups * POWhite * PCWhiteS4 * (1 - POWhite*PCWhiteS4)
    sdmoult <- sqrt(moultPups * POMoult * PCMoult * (1 - POMoult * PCMoult)) +
      S4Pups*POWhite*(1-PCWhiteS4)*(1 - POWhite*(1-PCWhiteS4))
    
    whiteCount <- rnorm(n = length(whitePups), mean=whitePups * POWhite + 
                          moultPups * POMoult * (1 - PCMoult)+ S4Pups*POWhite*PCWhiteS4, sd=sdwhite)
    moultCount <- rnorm(n = length(moultPups), mean= moultPups * POMoult * PCMoult + 
                          S4Pups*POWhite*(1-PCWhiteS4), sd=sdmoult)
      

  # white classified as white:
  WW <- whiteCount
  #No whitecoats misclassified as moulted
  WM <- 0
  # moulted classified as moulted
  MM <- moultCount
  # This is already acounted for in moultCount calculation for rnorm but needed for function to run
  MW <- 0
  
  }
    
  if (obsproc == "binom"){
    set.seed(seed)
    whiteCount <- rbinom(n = length(whitePups), size = round(whitePups), POWhite)
    moultCount <- rbinom(n = length(moultPups), size = round(moultPups), POMoult)
    # add in imperfect classification
    # white classified as white
    WW <- rbinom(n = length(whiteCount), size = whiteCount, PCWhite)
    # white classified as moulted
    WM <- round(whiteCount) - WW
    # moulted classified as moulted
    MM <- rbinom(n = length(moultCount), size = moultCount, PCMoult)
    # moulted classified as white
    MW <- round(moultCount) - MM  #correctly observed moulteds - correcty classified moulteds


  }
  

  # # simulated observation data of white and moulted:
  whiteObs <- WW + MW
  moultObs <- MM + WM

  wObs <- whiteObs[obsdays]
  mObs <- moultObs[obsdays]

   

  sim.out <- list( realBorn = realBorn, Pb = Pb,Pm = Pm, #totPups=totPups,
                  whiteCount=whiteCount,moultCount=moultCount,S4Count=S4Count, #Pw=Pw,Pm=Pm,
                  whiteObs = wObs, moultObs = mObs,#mleavers=mleavers,
                 whitePups=whitePups,S4Pups=S4Pups,moultPups=moultPups)
  list2env(sim.out, globalenv())
  return(sim.out)

}

sim_S4<-simSealsS4()


## check that whiteObs+moultObs= whiteCount+moultCount
#sim_S3$whiteCount[obsdays]+sim_S3$moultCount[obsdays]
#sim_S3$whiteObs+sim_S3$moultObs
####################################3

#        lty = c(1, 1), col = c("darkgrey", "red","blue"), lwd = 2,bty = "n", cex = 0.8)
