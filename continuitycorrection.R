
library(tigerstats)
dbinom(x = 47, size = 100, prob = .5)
## [1] 0.0665905

dnorm(x = 47, mean = 50, sd = 5)
## [1] 0.06664492

pnorm(q = 495.3367, mean = 500, sd = 0.02) - pnorm(q = 500.7207, mean = 500, sd = 0.02)
## [1] 0.06657389

pnormGC(c(495.3367,500.7207), "between", mean = 500, sd = 0.02, graph = TRUE)


probs2 = dbinom(1:500, size=500, prob=.95)
x<-probs2
#binomial distributions
barplot(probs2, names.arg=c(1:100), space=0, xlim=c(85,100), ylim=c(0,0.2))

#normal approximations
lines((c(1:100)-.5), dnorm(c(1:100), 95, sqrt(4.75)), type="h", lwd=2, col="red")
curve(dnorm(c(x+.5), mean=95, sd=sqrt(4.75)), from=85, to=101, xlim = c(85:100), add=T, col="blue")

#adjust the binomial probs and mean
probs2 = dbinom(1:100, size=100, prob=.65)
x <- probs2
barplot(probs2, names.arg=c(1:100), space=0, xlim=c(47,83), ylim=c(0,0.09))

lines((c(1:100)-.5), dnorm(c(1:100), 65, sqrt(22.75)), type="h", lwd=2, col="red")
curve(dnorm(c(x+.5), mean=65, sd=sqrt(22.75)), from=47, to=83, xlim = c(47:83), add=T, col="blue")
################################################################
#So a normal distribution won’t work when the probability p is close to 0 or 1 or when the number of trials n is small.
#p close to 1, larger N
probs2 = dbinom(1:100, size=100, prob=.95)
x<-probs2
barplot(probs2, names.arg=c(1:100), space=0, xlim=c(85,100), ylim=c(0,0.2))

lines((c(1:100)-.5), dnorm(c(1:100), 95, sqrt(4.75)), type="h", lwd=2, col="red")
curve(dnorm(c(x+.5), mean=95, sd=sqrt(4.75)), from=85, to=101, xlim = c(85:100), add=T, col="blue")



#p small, but small N
probs2 = dbinom(1:10, size=10, prob=.3)
x <- probs2
barplot(probs2, names.arg=c(1:10), space=0, xlim=c(0,10), ylim=c(0,0.3))

lines((c(1:10)-.5), dnorm(c(1:10), 3, sqrt(2.1)), type="h", lwd=2, col="red")
curve(dnorm(c(x+.5), mean=3, sd=sqrt(2.1)), from=0, to=10, xlim = c(0:10), add=T, col="blue")
      
#But if we adjust them things will improve. The left plot has a probability closer to 0.5 and in the right plot n⋅p
#is equal to 5. 

#prob=0.6
probs2 = dbinom(1:100, size=100, prob=.65)
x <- probs2
barplot(probs2, names.arg=c(1:100), space=0, xlim=c(47,83), ylim=c(0,0.09))

lines((c(1:100)-.5), dnorm(c(1:100), 65, sqrt(22.75)), type="h", lwd=2, col="red")

curve(dnorm(c(x+.5), mean=65, sd=sqrt(22.75)), from=47, to=83, xlim = c(47:83), add=T, col="blue")

# prob =0.5
probs2 = dbinom(1:10, size=10, prob=.5)
x <- probs2
barplot(probs2, names.arg=c(1:10), space=0, xlim=c(0,10), ylim=c(0,0.3))
#Adds a red line showing a normal approximation to the binomial distribution, using a 
#mean of 5 and a standard deviation of sqrt(2.5). The -.5 adjustment is a continuity correction.
lines((c(1:10)-.5), dnorm(c(1:10), 5, sqrt(2.5)), type="h", lwd=2, col="red")
#Adds a blue curve showing another normal approximation, shifted slightly to the right for better visualization
curve(dnorm(c(x+.5), mean=5, sd=sqrt(2.5)), from=0, to=10, xlim = c(0:10), add=T, col="blue")


### with Nborn 500, POWhite=0.98, POMoult=0.94
probs2 = dbinom(1:500, size=500, prob=.98)
x <- probs2
barplot(probs2, names.arg=c(1:500), space=0, xlim=c(0,500), ylim=c(0,0.3))
#Adds a red line showing a normal approximation to the binomial distribution, using a 
#mean of 500 and a standard deviation of 0.01. The -.5 adjustment is a continuity correction.
lines((c(1:500)-.5), dnorm(c(1:500), 500, 0.01), type="h", lwd=2, col="red")
#Adds a blue curve showing another normal approximation, shifted slightly to the right for better visualization
curve(dnorm(c(x+.5), mean=500, sd=0.01), from=0, to=500, xlim = c(0:500), add=T, col="blue")




x <- sim_S4$whitePups # Value to evaluate the cumulative distribution function at
mean <- mean(sim_S4$whiteCount)  # Mean of the normal distribution
sd <- model$sdwhite  # Standard deviation of the normal distribution

dnorm_cdf_approx <- function(x, mean, sd) {
  pnorm(x + 0.5, mean, sd) - pnorm(x - 0.5, mean, sd)
}


# Assuming you have values for x, mean, and sd
result <- pnorm(x + 0.5, mean, sd) - pnorm(x - 0.5, mean, sd)
print(result)




#######################
#log density function
ld <- nimbleFunction(
  run = function(x = double(), mean = double(), sd = double()) {
    returnType(double())
    return(log(dnorm(x, mean, sd)))
  }
)


#################################
# #defining the discretized pdf as the area under the curve from x-0.5 to x+0.5.
# dnormCDFw <- nimbleFunction(
# 
#   run = function(x = double(), mean = double(), sd = double(), log = integer(default = 0)) {
#     ans <- pnorm(x + 0.5, mean = mean, sd = sd) - pnorm(x -0.5, mean = mean, sd = sd)
#     if(log) return (log(ans))
#     return(ans)
#     returnType(double())
#   }
# )
# 
# 
# registerDistributions(list(dnormCDFw = list(BUGSdist = "dnormCDFw(mean,sd)", range = c(0, 500))))
# 
# dnormCDFm <- nimbleFunction(
# 
#   run = function(x = double(), mean = double(), sd = double(), log = integer(default = 0)) {
#     ans <- pnorm(x + 0.5, mean = mean, sd = sd) - pnorm(x-0.5, mean = mean, sd = sd)
#     if(log) return (log(ans))
#     return(ans)
#     returnType(double())
#   }
# )
# registerDistributions(list(dnormCDFm = list(BUGSdist = "dnormCDFm(mean,sd)", range = c(0, 500))))