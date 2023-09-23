### ===================
### HG sim with only 3 stages: whitecoat, moult, leave 
### ===================

simSealsS3 <- function(N = 500, 
   ndays = 123, amax = 123, 
   start = 32, mu.bday = 45, sd.bday = 10, 
   #mu.s4day = 17, sd.s4day = 5,
   mu.mlen = 6, sd.mlen = 2,
   mu.mday = 23, sd.mday = 5,
   mu.lday = 31.5, sd.lday = 7,
   POWhite = 0.98, POWhiteSD = 0.01,
   POMoult = 0.94, POMoultSD = 0.03,
   PCWhite = 1, 
   PCMoult = 1, PCMoultSD = 0.09, 
   #PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
   realizeObsPars = TRUE,
   obsdays = c(45,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110,115),
   popproc = "detcont", #detcont = deterministic and continuous, option for "stocdis" - stochastic and discrete
   obsproc = "binom",seed=2024) {  

   library(lubridate)

  #########
  # initialize matrices

  set.seed(seed)
  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  
  # convert birth parameters for lognormal dist
  lCV <- sd.bday/(mu.bday-start)
  lsd <- sqrt(log(lCV^2+1))
  lmu <- log(mu.bday-start)-(0.5*(lsd^2))
  
  # calculate probability of being born on each day
  # calculate probability of being born on each day
  Pb <- rep(0, ndays)
  
  for(t in 2:ndays){
    Pb[t] <- plnorm((t-start), lmu, lsd) - plnorm((t-start-1), lmu, lsd)
  }
  
  Pb <- Pb/sum(Pb)
  
  # calculate the probability of ending moult at each day post-S4
  
  Pm <- rep(0, amax)
  
  for (t in 1:amax){
    Pm[t] <- pnorm(t, mu.mday, sd.mday) - pnorm(t-1, mu.mday, sd.mday)
  }

  Pm <- Pm/sum(Pm)
  
  # calculate the probability of leaving at each age
  
  Pl <- rep(0, amax)
  
  for (t in 2:amax){
    Pl[t] <- pnorm(t, mu.lday, sd.lday) - pnorm(t-1, mu.lday, sd.lday)
  }
  
  Pl <- Pl/sum(Pl)

  
  # simulate perfect birth, moult, leave process
  if (popproc == "detcont"){
      
      realBorn <- Pb*N 
      
      whitePups <- rep(0, ndays)
      moultPups <- rep(0, ndays)
      totPups <- rep(0, ndays)
      for (d in 2:(ndays)){
        for (i in 1:d){
          age = d-i
          whitePups[d]<- whitePups[d] + realBorn[i] *(1 - sum(Pm[1:age])) * (1 - sum(Pl[1:age]))
          moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) * (1-sum(Pl[1:age]))
          totPups[d]<- whitePups[d]+moultPups[d]
        } # end for d
      } # end for i
      
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
    moultCount <- moultPups*POMoult
    
    # white classified as white:
    WW <- whiteCount*PCWhite
    #whitecoats classified as moulted is observed whitecoats - WW
    WM <- whiteCount-WW
    # moulted classified as moulted
    MM <- moultCount*PCMoult
    # moulted classified as white
    MW <- moultCount - MM

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
  
  if (obsproc == "binom_OBS"){
  set.seed(seed)
  whiteCount <- rbinom(n = length(whitePups), size = round(whitePups), POWhite)
  moultCount <- rbinom(n = length(moultPups), size = round(moultPups), POMoult)
  # white classified as white
  WW <- whiteCount
  # white classified as moulted
  WM <- 0
  # moulted classified as moulted
  MM <- moultCount
  # moulted classified as white
  MW <- 0
  
  }
  
  # simulated observation data of white and moulted: 
  whiteObs <- WW + MW
  moultObs <- MM + WM

  wObs <- whiteObs[obsdays]
  mObs <- moultObs[obsdays]
  

  sim.out <- list( realBorn = realBorn, Pb = Pb,Pm = Pm, Pl = Pl,totPups=totPups,
                  whiteCount=whiteCount,moultCount=moultCount,
                  whiteObs = round(wObs), moultObs = round(mObs),
                  whitePups = whitePups, moultPups = moultPups)
  #list2env(sim.out, globalenv())
  return(sim.out)

}

sim_S3<-simSealsS3()


