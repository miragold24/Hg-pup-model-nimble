### ===================
### HG pup simulation with 4 stages: whitecoat, moulting (S4) moulted, leave 
### options for various observation processes
### ===================

input_list = list(ndays = 123, amax=123,N=500, mu.bday = 45, sd.bday = 10,
                  mu.s4day = 17, sd.s4day = 5,
                  mu.mlen = 6, sd.mlen = 2,mu.mday=23,sd.mday=5,
                  mu.lday = 31.5, sd.lday = 7,
                  skew.bday=5,
                  POWhite = 0.98, 
                  POMoult = 0.94,
                  PCWhite = 1,
                  PCMoult = 0.82, 
                  obsdays = c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110))#1:123) #65:75)#c(65, 70, 75, 80,85,90, 95,100))#1:123) #c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110))
thisenv <-environment()
list2env(input_list, thisenv)


simSealsS4 <- function(N = 500, 
                       ndays = 123, amax = 123,
                       start = 32, mu.bday = 45, sd.bday = 10, 
                       skew.bday = 5,
                       mu.mday = 23, sd.mday = 5,
                       mu.lday = 31.5, sd.lday = 7,
                       mu.s4day=17,sd.s4day=5,
                       POWhite = 0.98, 
                       POMoult = 0.95, 
                       PCWhite = 1, 
                       PCMoult = 0.82, 
                       PCWhiteS4 = 0.75, 
                       obsdays = c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110),# 1:123, #65:75,#c(65, 70, 75, 80,85,90, 95,100), #1:123, #
                       obsproc = "binom_obs",seed=2024) {  

    set.seed(seed)
    #skew normal birth curve
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
  
   
  
  # generate pups and transitions

    realBorn <- 0
    
    #The t function transposes this matrix, so each row of realBorn represents a draw.
    realBorn <- as.numeric(t(rmultinom(1, size = N, prob = Pb)))
    realS4 <- matrix(rep(0, ndays * amax), nrow = ndays, ncol = amax)
    for (i in 1:ndays){
      
      while(sum(realS4[i,]) != sum(realBorn[i])){
        realS4[i,] <- rbinom(n = amax, size = realBorn[i], prob = PS4)
      }
      
    }
    
    s4vec <- rep(0, ndays)
    tmin <- min(which(realBorn!=0))
    
    #take upper triangular of realS4
    for (i in tmin:ndays){
      
      
      s4v <- 0
      rowmax <- i - tmin + 1
      
      for (j in tmin:i){
        for (k in 1:rowmax){
          
          s4v <- sum(s4v, realS4[j, k])
          
        } # end for k
        
        rowmax <- rowmax - 1
      }
      
      s4vec[i] <- s4v
    } # end for i
    
    cumS4<-s4vec
   
    
    realMoult = array(data = 0, dim = c(ndays,amax,amax))
    for (dob in 1:nrow(realS4)) {
      for (age in 1:ncol(realS4)) {
        realMoult[dob,age,] <- as.vector(rmultinom(1, size = realS4[dob,age], prob = Pm))
      }
    }
    
    realMoult2d = apply(realMoult, c(1,3), sum)
    dayMoult = matrix(0, nrow = ndays, ncol = amax)
    
    for(dob in 1:nrow(realMoult2d)){
      for(dom in 1:ncol(realMoult2d)){
        day_moult = min(dom + dob - 1,123)
        dayMoult[dob,day_moult] = dayMoult[dob,day_moult] + realMoult2d[dob,dom]
      }
    }
    
    cumMoult<-cumsum(apply(dayMoult, 2, sum))
    
    
    #
    # calculate pups that leave conditional on having moulted
    #
    realLeave = array(data = 0, dim = c(ndays,amax,amax))
    for (dob in 1:nrow(realMoult)) {
      for (age in 1:ncol(realMoult)) {
        realLeave[dob,age,] <- as.vector(rmultinom(1, size = realMoult2d[dob,age], prob = Pl)) #PmPl2$PlGivenm))
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
    
    # whitePups<- cumsum(realBorn) - cumS4
    # S4Pups <- cumS4-cumMoult
    # moultPups <- cumMoult - cumLeave
    # leavePups <-colSums(dayLeave)
    # allPups <- cumsum(realBorn) - cumLeave
    
    whitePups <- pmax(cumsum(realBorn) - cumS4, 0)
    S4Pups <- pmax(cumS4 - cumMoult, 0)
    moultPups <- pmax(cumMoult - cumLeave, 0)
    leavePups <- pmax(colSums(dayLeave), 0)
    allPups <- pmax(cumsum(realBorn) - cumLeave, 0)
    white1_4<-whitePups+S4Pups
  
  if (obsproc == "exact"){
    whiteObs <- whitePups*POWhite + moultPups*POMoult*(1-PCMoult)+
      S4Pups*POWhite*PCWhiteS4
    moultObs <- moultPups*POMoult*PCMoult +S4Pups*POWhite*(1-PCWhiteS4)
    
    whiteCount <- whiteObs
    moultCount <- moultObs
    
  }
  
  if (obsproc == "pois"){
    whiteO<-rep(0,ndays)
    moultO<-rep(0,ndays)
    whiteCount<-rep(0,ndays)
    moultCount<-rep(0,ndays)

      # whiteM <- whitePups* POWhite + S4Pups * POWhite*PCWhiteS4 + moultPups * POMoult*(1-PCMoult)
      # moultM <- S4Pups * POWhite*(1-PCWhiteS4) + moultPups * POMoult*PCMoult
      whiteO <- whitePups* POWhite + S4Pups * POWhite
      moultO <- moultPups * POMoult
      
      whiteCount <- rpois(n = length(whiteM), round(whiteM))
      moultCount <- rpois(n = length(moultM), round(moultM))
    }
    
  
  
  if (obsproc == "normappx_DT"){
    # POBS only
    
    # POBS AND CLASSIFICATION ERROR
    whiteCount <-rep(0,ndays)
    moultCount<-rep(0,ndays)
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

      whiteCount[i] <- abs(rnorm(n = length(whiteM[obsdays[i]]), mean=whiteM[obsdays[i]], sd=whiteS[i]))
      moultCount[i] <- abs(rnorm(n = length(moultM[obsdays[i]]), mean= moultM[obsdays[i]],sd=moultS[i]))
    }
  }
    
    if (obsproc == "normappxC_O"){
    ## CLASSIFICATION -> POBS
    sdwhiteC <- sqrt(moultPups * (1 - PCMoult) * (1 - (1 - PCMoult)) +
                      S4Pups * PCWhiteS4 * (1 - PCWhiteS4))

    sdmoultC <- sqrt(moultPups *  PCMoult * (1 - PCMoult) +
                      S4Pups*(1-PCWhiteS4)*(1 - (1-PCWhiteS4)))

    sdwhiteO <- sqrt(whitePups*POWhite * (1-POWhite))
    sdmoultO  <- sqrt(moultPups*POWhite * (1-POWhite))

    whiteC <- rnorm(n = length(whitePups), mean=whitePups +
                          moultPups * (1 - PCMoult)+
                          S4Pups*PCWhiteS4, sd=sdwhiteC)
    moultC <- rnorm(n = length(moultPups), mean= moultPups*PCMoult, sd=sdmoultC)

    whiteCount<- rnorm(n=length(whiteC),mean=whiteC*POWhite,sd=sdwhiteO)
    moultCount <- rnorm(n=length(moultC),mean=moultC*POMoult,sd=sdmoultO)
    
    }
    
    if (obsproc == "normappx_EJ"){
      #PObs=0.95
      whiteDetectSD <- sqrt(whitePups*POWhite* (1-POWhite))
      moultDetectSD <- sqrt(moultPups*POMoult * (1-POMoult))

      whiteDetect <- rnorm(length(whitePups),
                           mean = whitePups*POWhite,
                           sd = whiteDetectSD)

      moultDetect <- rnorm(length(moultPups),
                           mean = moultPups*POMoult,
                           sd = moultDetectSD)

      whiteClassSD <- sqrt(abs(whiteDetect)*PCWhite * (1-PCWhite))
      moultClassSD <- sqrt(abs(moultDetect)*PCMoult * (1-PCMoult))
      WW <- rnorm(length(whiteClassSD),
                  mean = whiteDetect*PCWhite,
                  sd = whiteClassSD)
      WM <- whiteDetect - WW
      MM <- rnorm(length(moultDetect),
                  mean = moultDetect*PCMoult,
                  sd = moultClassSD)
      MW <- moultDetect - MM

      whiteCount <- abs(WW + MW)
      moultCount <- abs(MM + WM)
    }
  
  
  if (obsproc == "binom_obs"){
    
    white1_4<-whitePups+S4Pups
    whiteCount <- rbinom(n=length(white1_4), size =(white1_4), POWhite)
    moultCount <- rbinom(n=length(moultPups), size =(moultPups), POMoult)
    # whiteCount <- rbinom(n=length(whitePups), size =(whitePups), POWhite)
    # whiteS4 <- rbinom(n=length(S4Pups), size =(S4Pups), POWhite)
    # whiteCount <- white+whiteS4
    # moultCount <- rbinom(n=length(moultPups), size =(moultPups), POMoult)
    # # pups counted as white or moulted come from a binomial with size= all pups classified
    # # into that age stage and probability of observing them (which is influenced by their classification)
    # whiteC <- whitePups + moultPups*(1-PCMoult)+ S4Pups*PCWhiteS4
    # moultC <- moultPups*PCMoult +S4Pups*(1-PCWhiteS4)
    # whiteCount <- rbinom(n=length(whiteC), size =(round(whiteC)), POWhite)
    # moultCount <- rbinom(n=length(moultC), size =(round(moultC)), POWhite)
  }
    
    


    if (obsproc == "binom_CO"){
    #White detected as white stage 1-3 and S4
    WWS4  <- rbinom(n = length(S4Pups), size = round(S4Pups), PCWhiteS4)
    #All whitecoats correctly classified as whitecoats:
    WW <- whitePups + WWS4
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
    moultCount <- rbinom(n = length(MM+WMS4), size = round(MM+WMS4), POMoult)

    if (any(whitePups + WWS4 > whitePups + S4Pups)) {
      warning("Warning: The total correctly classified white 1-4 is greater than the total actual white 1-4.")
    }
    
    }
    
    if (obsproc == "binom_OC"){
    whiteDetected <- rbinom(length(whitePups), whitePups, POWhite)
    whiteS4Detected <- rbinom(length(S4Pups),S4Pups,POWhite)
    moultDetected <- rbinom(length(moultPups), moultPups, POMoult)
    
    WW <- whiteDetected+whiteS4Detected
    WM <- rbinom(length(whiteS4Detected), whiteS4Detected, 1-PCWhiteS4)
    MM <- rbinom(length(moultDetected), moultDetected, PCMoult)
    MW <- rbinom(length(moultDetected), moultDetected, 1-PCMoult)
    whiteCount <- WW + MW
    moultCount <- MM + WM
    }
    
    
  # simulated observation data of white and moulted:
  wObs <- whiteCount[obsdays]
  mObs <- moultCount[obsdays]
  
  sim.out <- list( realBorn = realBorn, Pb = Pb,PS4=PS4,Pm = Pm,
                    whiteCount=whiteCount,moultCount=moultCount,
                   wObs = wObs, mObs = mObs, white1_4=whitePups+S4Pups,
                   whitePups=whitePups,S4Pups=S4Pups,moultPups=moultPups,leavePups=leavePups)
  list2env(sim.out, globalenv())
  return(sim.out)
  
}

sim_S4<-simSealsS4()


#check wObs+mObs always less than total available pups present on colony
#any(wObs+mObs>whitePups+S4Pups+moultPups)
####################################
# visualize actual pups on colony and counted pups: 

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


