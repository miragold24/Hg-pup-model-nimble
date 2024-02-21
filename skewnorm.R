dskewnorm <- function(x, mu, sigma, alpha){
  2*dnorm(x, mu, omega^2) * pnorm(alpha*((x-mu)/omega))
}

dskewnorm(0, 0, 1, 0)

#homemade pdf
x <- seq(0, 10, by =0.1)
y <- dskewnorm(x, mu = 5, sigma = 1, alpha = 4)
plot(x, y/sum(y))

# probability density function
x <- seq(0, 10, by = 0.1)
#plot(x,dsn(x, xi=5, omega=1, alpha=4, log=FALSE))
s <- dsn(x, xi=5, omega=1, alpha=4, log=FALSE)
lines(x, s/sum(s))

# cumulative density function
x <- seq(0, 10, by = 0.1)
#plot(x,psn(x, xi=5, omega=1, alpha=4, log=FALSE))
p <- psn(x, xi=5, omega=1, alpha=4, log=FALSE)
deltap <- p[2:length(p)]- p[1:(length(p)-1)]
lines(x[2:(length(x))], deltap, col = "red")


plot(x, s/sum(s), type = "l")
lines(x[1:(length(x)-1)], deltap)
lines(x[2:(length(x)-1)], deltap)
lines(x[2:(length(x))], deltap)

#https://www.sciencedirect.com/science/article/pii/S209012321000069X


# comparing lognormal and skewnormal

ndays = 123
start = 32
mu.bday = 45
sd.bday = 10


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
Pb.ln <- rep(0, ndays)

for(t in 2:ndays){
  Pb.ln[t] <- plnorm((t-start), lmu, lsd) - plnorm((t-start-1), lmu, lsd)
}

Pb.ln <- Pb.ln/sum(Pb.ln)

library(sn)

skew.bday <- 5


# using package sn

mu.bday <- 45
sd.bday <- 10
alpha.bday <- 5

Pb.sn <- rep(0, ndays)

for(t in 2:ndays){
  Pb.sn[t] <- psn(x = t, xi = mu.bday, omega = sd.bday, alpha = skew.bday) - 
    psn(x = t-1, xi = mu.bday, omega = sd.bday, alpha = skew.bday)
}

Pb.sn <- Pb.sn/sum(Pb.sn)




Pb.sn <- rep(0, ndays)

for(t in 2:ndays){
  Pb.sn[t] <- pskew_normal(q = t, mu = mu.bday, sigma = sd.bday, alpha = skew.bday) - 
    pskew_normal(q = t-1, mu = mu.bday, sigma = sd.bday, alpha = skew.bday)
}

Pb.sn <- Pb.sn/sum(Pb.sn)


plot(Pb.ln, col = "black", type = "l")
lines(Pb.sn, col = "blue")


# check what the actual mean birthday is

checkBday <- function(Pb.vec){
  
  numbers <- round(10000*Pb.vec)
  bday.vec <- vector()
  for (i in 1:length(Pb.vec)){
    bday.vec <- c(bday.vec, rep(i, numbers[i]))
  }
  return(c(mean(bday.vec), sd(bday.vec)))}

checkBday(Pb.sn)



#### compare sn and dskewnormal given mu, sigma, alpha

library(sn)

mu <- 45
sigma <- 10
# recover normal by testing with alpha =0, no skew
alpha <-5

# convert to params for skew normal

delta <- alpha/sqrt(1+(alpha^2))

omega <- sqrt(sigma^2/(1-(2*delta^2)/pi))

xi <- mu - omega*delta*sqrt(2/pi)

# PDF of skewnormal: 
dskewnorm <- function(x, xi, omega, alpha){
  (2/omega)*dnorm((x-xi)/omega) * 
    pnorm(alpha*((x-xi)/omega))
}
###
### calculate Pb at each 100th of a day interval and sum
###

#visualize against packages

Pb.dskewnorm <- dskewnorm(x = 1:123, xi = xi,  omega = omega, alpha = alpha)
#Pb.dskewnorm <- Pb.dskewnorm/max(Pb.dskewnorm)
plot(Pb.dskewnorm/sum(Pb.dskewnorm), type = "l")
lines(dsn(1:123, xi, omega, alpha)/sum(dsn(1:123, xi, omega, alpha)), col = "blue")
lines(dskew_normal(1:123, mu, sigma, alpha)/sum(dskew_normal(1:123, mu, sigma, alpha)), col = "red")  
lines(dskew_normal(1:123, xi = xi, omega = omega, alpha = alpha)/sum(dskew_normal(1:123, xi = xi, omega = omega, alpha = alpha)), col = "green")

# now need to make this a cdf and make sure it matches

ndays <- 123
nslice <- 1000
pborn <- rep(0, ndays)
for (i in 1:ndays){
  for (j in 1:nslice){
    dayslice <- i + (j/nslice)
    pborn[i] <- pborn[i] + dskewnorm(dayslice, xi, omega, alpha)
  }
}
plot(1:123, pborn/sum(pborn))

plot(1:123,cumsum(pborn/sum(pborn)), type = "l")
lines(psn(1:123, xi, omega, alpha), col = "blue")
lines(pskew_normal(1:123, mu, sigma, alpha), col = "red")  
lines(pskew_normal(1:123, xi = xi, omega = omega, alpha = alpha), col = "green")


###

# check output from TMB

pbtmb <- as.numeric(rep$value[which(names(rep$value) == "pborn1")])

plot(1:ndays, cumsum(pbtmb), type = "l")
lines(psn(1:123, xi, omega, alpha), col = "blue")
lines(pskew_normal(1:123, mu, sigma, alpha), col = "red")  
lines(pskew_normal(1:123, xi = xi, omega = omega, alpha = alpha), col = "green")


#### check bday

mean(rep(1:123, times = round(10000*pbtmb)))
mean(rep(1:123, times = round(dsn(1:123, xi, omega, alpha)*10000)))
#################################################################################################
# function for calculating the skew normal
# convert to params for skew normal

mu <- 45
sigma <- 10
# recover normal by testing with alpha =0, no skew
alpha <-5

delta <- alpha/sqrt(1+(alpha^2))

omega <- sqrt(sigma^2/(1-(2*delta^2)/pi))

xi <- mu - omega*delta*sqrt(2/pi)

# PDF of skewnormal: 
dskewnorm <- function(x, xi, omega, alpha){
  (2/omega)*dnorm((x-xi)/omega) * pnorm(alpha*((x-xi)/omega))
}

dskewnorm <- function(x){
  mu <- 45
  sigma <- 10
  # recover normal by testing with alpha =0, no skew
  alpha <-5
  
  delta <- alpha/sqrt(1+(alpha^2))
  
  omega <- sqrt(sigma^2/(1-(2*delta^2)/pi))
  
  xi <- mu - omega*delta*sqrt(2/pi)
  (2/omega)*dnorm((x-xi)/omega) * pnorm(alpha*((x-xi)/omega))
}


  
# Initialize Pb as a vector of the same length as ndays
Pb <- numeric(ndays)
# Number of intervals per day
int <- 100

for (day in 1:ndays) {
  # Initialize the probability for this day
  Pb_day <- 0
  
  for (i in 1:int) {
    # Define the start and end of the interval
    start <- (day - 1) + (i - 1) / int
    end <- start + 1 / int
    
    # Integrate over the interval and add to the probability for this day
    Pb_day <- Pb_day + integrate(dskewnorm, start, end)$value
  }
  
  # Store the probability for this day
  Pb[day] <- Pb_day
}




sum(Pb)
lines(Pb,col="green")

sum(Pb)


### nimble function?
dskewnorm <- nimbleFunction(
  run = function(x = double(0)) {
    mu <- 45
    sigma <- 10
    alpha <- 5
    
    delta <- alpha/sqrt(1+(alpha^2))
    
    omega <- sqrt(sigma^2/(1-(2*delta^2)/pi))
    
    xi <- mu - omega*delta*sqrt(2/pi)
    (2/omega)*dnorm((x-xi)/omega) * pnorm(alpha*((x-xi)/omega))
    
    Pb <- numeric(ndays)
    # Number of intervals per day
    int <- 100
    
    for (day in 1:ndays) {
      # Initialize the probability for this day
      Pb_day <- 0
      
      for (i in 1:int) {
        # Define the start and end of the interval
        start <- (day - 1) + (i - 1) / int
        end <- start + 1 / int
        
        # Integrate over the interval and add to the probability for this day
        Pb_day <- Pb_day + nimIntegrate(dskewnorm, start, end,param=x)$value
      }
      
      # Store the probability for this day
      Pb[day] <- Pb_day
    }
    
    
  }
)




################
#visualize in R: 
#Convert paramters for skew normal
# mu <- 45
# sigma <- 10
# alpha <- 5 #skew.bday
# alpha1=5
# #x=1:123
# dskewnorm <- function(x, xi, omega, alpha1){
#   #convert paramaters for skew normal
#   delta <- alpha/sqrt(1+(alpha^2)) #location
# 
#   omega <- sqrt(sigma^2/(1-(2*delta^2)/pi)) #scale
# 
#   xi <- mu - omega*delta*sqrt(2/pi) #skewness
# 
#   Pb<- (2/omega)*dnorm((x-xi)/omega) *
#     pnorm(alpha*((x-xi)/omega))
# }
# Pb.dskewnorm <- dskewnorm(x = 1:123, xi = xi,  omega = omega, alpha1 = alpha)
# 
# #lines(1:123,cumsum(Pb.dskewnorm/sum(Pb.dskewnorm)), type = "l",col='green')
# ############################
# 
# SkewNorm <- function(mu.bday, sd.bday, skew.bday, ndays) {
#   delta <- skew.bday/sqrt(1+(skew.bday^2))
#   omega <- sqrt((sd.bday^2)/(1-(2*(delta^2))/3.14159))
#   xi <- mu.bday - omega*delta*sqrt(2/3.14159)
#   
#   pbornslices <- rep(NA, ndays*10)
#   slices <- seq(0, ndays, length.out = ndays*10)
#   start <- seq(0, length(pbornslices), by = 10)
#   # calculate pborn at each 100th of a day interval
#   for (i in 1:(ndays*100)){
#     pbornslices[i] = 2/omega*dnorm(x = (slices[i] - xi)/omega, mean = 0, sd = 1)*pnorm(skew.bday*((slices[i]-xi)/omega))
#   }
#   
#   pborn <- rep(NA, ndays)
#   # add up the 100 slices for each day of the season
#   for(i in 1:ndays){
#     pborn[i] <- sum(pbornslices[start[i]:(start[i]+10)])
#   }
#   
#   Pb_sum <- pborn/sum(pborn)
#   return(Pb_sum)
# }

# # Paramater values:
# mu.bday <- 45
# sd.bday <- 10
# skew.bday <- 5
# ndays <- 123
# 
# result <- SkewNorm(mu.bday, sd.bday, skew.bday, ndays)
# lines(result,col='red')
# print(result)

#SkewNormNim <- nimbleRcall(function(mu.bday=double(),sd.bday=double(),skew.bday=double(),ndays=double()){}, Rfun = 'SkewNormNim',
#returnType = double())
#skew<-compileNimble(SkewNormNim)
############################
# #convert paramaters for skew normal
# delta <- alpha1/sqrt(1+(alpha1^2)) #location
# 
# omega<-sqrt(sigma^2/(1-(2*delta^2)/pi)) #scale
# 
# xi <- mu - omega*delta*sqrt(2/pi) #skewness
# 
# x=1:123
# PbSkew <- nimbleFunction(
#   run = function(x=double(),xi=double(),omega=double(),alpha1=double(),log= integer(default = 0)){
#     Pb=rep(0,123)  #numeric(length(x))
#     
#     delta <- alpha1/sqrt(1+(alpha1^2)) #location
#     
#     omega<-sqrt(sigma^2/(1-(2*delta^2)/pi)) #scale
#     
#     xi <- mu - omega*delta*sqrt(2/pi) #skewness
#     
# Pb <- (2/omega)*dnorm((x-xi)/omega) * 
#     pnorm(alpha1*((x-xi)/omega))
# 
# Pb<-Pb/sum(Pb)
# returnType(double())
# return(Pb)
# 
# })
# 
# Pbskew <- nimbleRcall(function(x=double(),xi=double(),omega=double(),alpha1=double()){}, Rfun = 'dskewnorm',
#                      returnType = double(1))
# registerDistributions(list(PbSkew = list(BUGSdist = "PbSkew(xi,omega,alpha)",range = c(0, 500))))
# 
# 
# Pbskew<-PbSkew(x=1:123,xi,omega,alpha)

# Pbskewnim<- compileNimble(PbSkew)
# # Compile the nimble function
# compiled_PbSkew <- compileNimble(PbSkew)
# 

# plot(cumsum(Pbskew/sum(Pbskew)), type = "l",col='green')
# lines(cumsum(Pb_sum/sum(Pb_sum)), type = "l",col='red')

# Pb.dskewnorm <- dskewnorm(x = 1:123, xi = xi,  omega = omega, alpha1 = alpha)

# Nim manual ch 12 cases in which you’ll need to explicitly call registerDistributions are when you want to do any of the following:
#   
#1) provide alternative parameterizations,
#2) indicate a distribution is discrete, and
#3) provide the range of possible values for a distribution.
#registerDistributions(list(PbSkew = list(BUGSdist = "PbSkew(x,mu,sigma,alpha)"),pqAvail = TRUE,range=c(0,123)))

