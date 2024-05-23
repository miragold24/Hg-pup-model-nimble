### Outputs/plots for MCMC runs / visualizations
###############################################
library(coda)
library(MCMCvis)
library(bayesplot)
###
###visualize priors
####

Nborn.prior<-rgamma(10000,shape=485, scale=1)
quantile(Nborn.prior)
plot(density(Nborn.prior))
mean(Nborn.prior)
Nborn.prior<-rgamma(10000,shape=4950, rate=10)
quantile(Nborn.prior)
mean(Nborn.prior)
lines(density(Nborn.prior))

skew.prior<-runif(10000,4,8)
quantile(skew.prior)
plot(density(skew.prior))


mubday.prior<-runif(10000,33,60)
quantile(mubday.prior)
plot(density(mubday.prior))


sd.prior<-runif(10000,7,12)
quantile(sd.prior)
plot(density(sd.prior))

### shifted beta priors 
# Set parameter values
alpha.m=2
beta.m=2
shift.m=10 # Shift the distribution to start at 33
# The "scale" parameter in a shifted beta distribution
# multiplies all values in the original beta distribution by 
# a constant factor. 
scale.m=40


# sd.bday
# alpha.sd=4  # Shape parameter (alpha > 0)
# beta.sd=4    # Shape parameter (beta > 0)
# shift.sd=6   # Shift the distribution to start at 6
# scale.sd=5  # 'stretch' to 11


alpha.sd <- 2
beta.sd <- 2
shift.sd <- 6
scale.sd <- 5

alpha.m <- 2
beta.m <- 2
shift.m <- 20
scale.m <- 30


alpha.s <- 2
beta.s <- 2
shift.s <- 6
scale.s <- 5


alpha.s=2
beta.s=2
shift.s=3
scale.s=5
# Set parameter values for the Poisson prior
lambda <- 500


# Number of samples
n_samples <- 1000


# Generate random samples from the Poisson distribution for the Nborn prior
Nborn_prior <- rpois(n_samples, lambda)
# Generate random samples from the shifted scale beta distribution for sd.bday
sdbday.prior <- rbeta(n_samples, alpha.sd, beta.sd) * scale.sd + shift.sd
quantile(sdbday.prior)
# Generate random samples for mu.bday and skew.bday
mubday.prior <- rbeta(n_samples, alpha.m, beta.m) * scale.m + shift.m
quantile(mubday.prior)

skew.prior <- rbeta(n_samples, alpha.s, beta.s) * scale.s + shift.s
quantile(skew.prior)
# Plot histograms
par(mfrow=c(1,1))
# hist(shifted_scale_samples, main = "Shifted Scale Beta for sd.bday", xlab = "Values")
# hist(mu_shifted_samples, main = "Shifted Scale Beta for mu.bday", xlab = "Values")
# hist(skew_shifted_samples, main = "Shifted Scale Beta for skew.bday", xlab = "Values")

# Plot densities
# Plot density
plot(density(Nborn_prior), main = "Nborn Prior")
plot(density(sdbday.prior), main = paste("sd.bday prior (shift =", shift.sd, ", scale =", scale.sd, ", alpha =", alpha.sd, ", beta =", beta.sd, ")", sep = ""))
plot(density(mubday.prior), main = paste("mu.bday prior (shift =", shift.m, ", scale =", scale.m, ", alpha =", alpha.m, ", beta =", beta.m, ")", sep = ""))
plot(density(skew.prior), main = paste("skew.bday prior (shift =", shift.s, ", scale =", scale.s, ", alpha =", alpha.s, ", beta =", beta.s, ")", sep = ""))
##mcmc output 
mcmc<-mcmc.out$samples

#Pair plots to check correlation:

mcmc_pairs(mcmc)

#coda::gelman.diag(mcmc[,"sd.bday"])
#look at trace plots
#plot(mcmc[,"Nborn"])
#trace and posterior density of params:
PR=cbind(Nborn.prior,mubday.prior,sdbday.prior,skew.prior)
MCMCtrace(object = samples,
          pdf = FALSE, 
          type='density',
          ind = TRUE, # separate density lines per chain
          Rhat=T,
          n.eff=T,
          exact = TRUE,
          priors = PR,
          #gvals = GV,
          params = "all")
# generating values for each parameter used to simulate data
GV <- c(-10, -5.5, -15)
MCMCtrace(object = samples,
          pdf = FALSE, 
          ind =T, # separate density lines per chain
          Rhat=F,
          n.eff=F,
          # gvals = GV
          params = "all")

#bayesplot::mcmc_dens_overlay(samples, pars = c('nInd','mu.bday','sd.bday','skew.bday'))

#95% credible intervals:
quantile(mcmc$chain1[,'Nborn'],probs=c(2.5,97.5)/100)


# check means:
nInd.samples<-  c(samples$chain1[,'nInd'], 
                   samples$chain2[,'nInd'])
mean(nInd.samples)

 start.samples<-  c(mcmc$chain1[,'start'], 
                     mcmc$chain2[,'start'])
 mean(start.samples)
 
mu.bday.samples<-  c(samples$chain1[,'mu.bday'], 
                     samples$chain2[,'mu.bday'])
 mean(mu.bday.samples)
 
sd.bday.samples<-  c(samples$chain1[,'sd.bday'], 
                     samples$chain2[,'sd.bday'])
mean(sd.bday.samples)


#stochastic model visualizaiton: 
# what is best distribution for births of pups over season?

par(mfrow=c(1,1))
### See if poisson approximates binomial using simulated data 
### (large N and small P is where this works: poisson 'limit' of binomial) 
pois_b<-rpois(123,500*Pb)
binom_b<-rbinom(123,500,Pb)
plot(pois_b,typ='l')
lines(binom_b,col='red')

### compare with multinomial
multi_b<-rmulti(1,500,Pb)
lines(multi_b,col='blue')
# plot(rmulti(1,500,Pb),type="l",ylab="pup count",xlab="day of season")


#note these are not normally distributed so normality cannot be assumed,
# doesn't this mean that normal approximation is not a good choice? 


#visualize actual pups on colony vs counted pups:
plot(whitePups+S4Pups,typ='l',col="purple",ylab='Count',xlab='day of season',main='actual pups present vs counted pups')
lines(whiteCount,col='orange')
lines(moultPups,col='red')
lines(moultCount,col='blue')

legend("topright", legend = c("whitePresent", "whiteCounted", "moultPresent","moultCounted"), col = c("purple", "orange", "red",'blue'), lty = 1)
