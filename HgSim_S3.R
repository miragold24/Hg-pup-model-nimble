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
   POMoult = 0.95, POMoultSD = 0.03,
   PCWhite = 1, 
   PCMoult = 0.82, PCMoultSD = 0.09, 
   #PCWhiteS4 = 0.75, PCWhiteS4SD = 0.01,
   obsdays = c(35,55, 60, 65, 70, 75, 80,85,90, 95,100,105,110),
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
    N <- 500
      realBorn <- Pb*N 
      
      whitePups <- rep(0, ndays)
      moultPups <- rep(0, ndays)

     # wleavers <- rep(0,ndays)
      #mleavers <- rep(0,ndays)
     # #totPups <- rep(0, ndays)
      for (d in 2:(ndays)){
        for (i in 1:d){
          age = d-i
          #calculate whitePups where no whitePups can leave
          whitePups[d]<- whitePups[d] + realBorn[i] * (1 - sum(Pm[1:age])) #* (1-sum(Pl[1:age]))
          # Calculate moultPups
          moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) *(1-sum(Pl[1:age]))
          # moultPups[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:age]) * (1-sum(Pl[1:age]))
          # #moultPups1[d] <- moultPups[d] + realBorn[i] * sum(Pm[1:thisage]) * (1-sum(PmPl2$PlGivenm[1:thisage]))


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
  
  if (obsproc == "stoc_binom"){
    set.seed(seed)
    whiteCount <- rbinom(n = length(whitePups), size = whitePups, POWhite)
    moultCount <- rbinom(n = length(moultPups), size = moultPups, POMoult)
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
  # 
  # ### USE THIS FOR OBSERVATIONS WITH OBSERVATION ERROR ONLY (basic model)
  # whiteObs <- whiteCount
  # moultObs <- moultCount
  # simulated observation data of white and moulted:
  whiteObs <- WW + MW
  moultObs <- MM + WM

  wObs <- whiteObs[obsdays]
  mObs <- moultObs[obsdays]

   

  sim.out <- list( realBorn = realBorn, Pb = Pb,Pm = Pm, Pl = Pl, #totPups=totPups,
                  whiteCount=whiteCount,moultCount=moultCount, #Pw=Pw,Pm=Pm,
                  whiteObs = wObs, moultObs = mObs,#mleavers=mleavers,
                 whitePups=whitePups,moultPups=moultPups)
  list2env(sim.out, globalenv())
  return(sim.out)

}

sim_S3<-simSealsS3()


## check that whiteObs+moultObs= whiteCount+moultCount
#sim_S3$whiteCount[obsdays]+sim_S3$moultCount[obsdays]
#sim_S3$whiteObs+sim_S3$moultObs
####################################3
#VISUALIZE 

# #check plots are the same for realBorn (set.seed)
# plot(sim_S3_binom$realBorn)
# lines(sim_S3_binom$realBorn,lty=2,col='green')
# 
# plot(sim_S3_binom$whitePups)
# lines(sim_S3_binom$moultPups,lty=2,col='red')
# 
# plot(sim_S3_exact$whiteCount)
# lines(sim_S3_binom$whiteCount,lty=2,col='red')
# sim_S3_exact$whiteCount[obsdays]
# sim_S3_exact$whitePups[obsdays]
# sim_S3_binom$whiteCount[obsdays]
# 

# 
# ymax <- max(mleavers)
# plot(whitePups, type = "l",ylim=c(0,500), ylab = "Count", 
#      xlab = "Day of Season", main = "Counts based on leave|moulted",
#      col = "darkgrey", lwd = 2)
# #lines(whitePups, col = 'orange', lty = 2, lwd = 2)
# #lines(moultPups1, col = 'purple', lwd = 2)
# lines(moultPups, col = 'red', lty = 1, lwd = 2)
# lines(mleavers, col = 'blue', lty = 1, lwd = 2)
# legend("topleft", legend = c("White pups",  
#                              "Moult pups (leave|moulted)","mleavers"), 
#        lty = c(1, 1), col = c("darkgrey", "red","blue"), lwd = 2,bty = "n", cex = 0.8)
