calculate_PmPl <- function(amax, mu.lday, sd.lday, mu.mday, sd.mday) {
  
  #Function for the density of leavers at some given age in continuous time
  #Basically just a truncated normal density
  #The normalization constant ensures that the total probability of the density 
  # function over the entire range of ages sums up to 1. It is calculated by 
  # subtracting the cumulative probabilities of the left truncation point (0) and 
  # the right truncation point (amax) from 1. 
  # This accounts for the truncation of the normal distribution and ensures that the density function is properly scaled.
  
  dleave <- function(age, amax, mu.lday, sd.lday){
    norm_const <- 1 - pnorm(0, mu.lday, sd.lday) - pnorm(amax, mu.lday, sd.lday, lower.tail = FALSE)
    res <- dnorm(age, mu.lday, sd.lday) / norm_const
    return(res)
  }
  
  #Function for the proportion of seals that are moulted at some age in continuous time
  # Basically just a truncated normal cdf
  pmoulted <- function(age, amax, mu.mday, sd.mday){
    norm_const <- 1 - pnorm(0, mu.mday, sd.mday) - pnorm(amax, mu.mday, sd.mday, lower.tail = FALSE)
    res <- (pnorm(age, mu.mday, sd.mday) - pnorm(0, mu.mday, sd.mday)) / norm_const
    return(res)
  }
  
  #Function for non-normalized density of moulted leavers (i.e., leaving and moulted)
  dmoulted_leave <- function(age, amax, mu.lday, sd.lday, mu.mday, sd.mday){
    res <- dleave(age, amax, mu.lday, sd.lday) * pmoulted(age, amax, mu.mday, sd.mday)
    return(res)
  }
  
  #Calculate the daily prob of leaving and moulted
  #integrate over proportion of pups moulted that will leave on given day given truncated 0-amax 
  # $value accesses the definite integral of dmoulted_leave function for specified interval
  norm_const <- integrate(dmoulted_leave, 0, amax, amax = amax, mu.lday = mu.lday, 
                          sd.lday = sd.lday, mu.mday = mu.mday, sd.mday = sd.mday)$value
  PmPl <- numeric(amax)
  for(age in 1:amax){
    res <- integrate(dmoulted_leave, age - 1, age, amax = amax, mu.lday = mu.lday, 
                     sd.lday = sd.lday, mu.mday = mu.mday, sd.mday = sd.mday)
    PmPl[age] <- res$value  
  }
  
  #Function for non-normalized density of non-moulted leavers 
  # (i.e., leaving and non-moulted)
  dnmoulted_leave <- function(age, amax, mu.lday, sd.lday, mu.mday, sd.mday){
    res <- dleave(age, amax, mu.lday, sd.lday) * (1 - pmoulted(age, amax, mu.mday, sd.mday))
    return(res)
  }
  
  #Calculate the daily prob of non-moulters leaving
  #  #In next line, this is actually just the complement of the last norm_const
  #  # but recalcualte here for robustness
  # norm_const <- integrate(dnmoulted_leave, 0, amax, amax = amax, mu.lday = mu.lday,
  #     sd.lday = sd.lday, mu.mday = mu.mday, sd.mday = sd.mday)$value
  PnmPl <- numeric(amax)
  for(age in 1:amax){
    res <- integrate(dnmoulted_leave, age - 1, age, amax = amax, mu.lday = mu.lday, 
                     sd.lday = sd.lday, mu.mday = mu.mday, sd.mday = sd.mday)
    PnmPl[age] <- res$value  
  }
  
  res <- data.frame(Age = 1:amax - 0.5, 
                    PlAndm = PmPl, 
                    PlAndnm = PnmPl,
                    
                    PlGivenm = PmPl / norm_const, 
                    PlGivennm = PnmPl / (1 - norm_const))
}
# Parameters
amax <- 123
age=1:amax-0.5
mu.lday <- 31.5
sd.lday <- 7
mu.mday <- 23
sd.mday <- 5

PmPl2 <- calculate_PmPl(amax, mu.lday, sd.lday, mu.mday, sd.mday)
