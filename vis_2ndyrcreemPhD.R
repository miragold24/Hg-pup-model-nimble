### visualizations for talks


#check simulation: stocdis, norm appx before source
#source("HgSim_s4.R")




#### DETERMINISTIC 

par(mfrow=c(1,1))
plot(sim_S4$realBorn,typ="l",ylab="count",xlab="day of season",main='Daily pups born over season')
max_pup<-which.max(sim_S4$realBorn)
abline(v=max_pup,col="red")

# Plot input to births -Pb

plot(sim_S4$Pb,type="l",xlab="day of season",ylab="birth probability",main="Birth curve")
max_index <- which.max(sim_S4$Pb)
abline(v = max_index, col = "red")

#######
# Define the start date and the number of days
start_date <- as.Date("2024-10-01")
num_days <- 123

# Generate the sequence of dates
date_sequence <- seq(start_date, by = "day", length.out = num_days)

# Calculate the indices for the first, middle, and last dates
first_index <- 1
middle_index <- 39
last_index <- num_days

# Extract the first, middle, and last dates
dates_to_display <- c(date_sequence[first_index], date_sequence[middle_index], date_sequence[last_index])

# Plot input to births - Pb
plot(sim_S4$Pb, type = "l", lwd = 2, xlab = "day of season", 
     ylab = "Birth probability", main = "Birth curve", 
     xlim = c(1, num_days), 
     xaxt = "n")  # Suppress automatic x-axis

# Add custom x-axis labels for the first, middle, and last dates
axis(side = 1, at = c(first_index, middle_index, last_index), labels = format(dates_to_display, "%b %d"))
# Add a vertical line at the maximum value of Pb
max_index <- which.max(sim_S4$Pb)
abline(v = max_index, col = "red",lwd="2")
#mtext(text = format(date_sequence[max_index], "%b %d"), side = 1, line = 1, at = max_index)


#check if this is around 40-45 (40 with skew norm)
# start_date <- as.Date("2024-10-01")
# end_date <- as.Date("2025-11-08")
# 
# # Calculate the difference in days
# num_days <- as.numeric(difftime(end_date, start_date, units = "days"))
# num_days
###############################


plot(sim_S4$PS4,type="l",col='blue',xlab="age",ylab="Probability density",xlim=c(0,60))#,main="Probability distribution of transitions across age classes")
lines(sim_S4$Pm,col='green')
lines(PmPl2$PlGivenm,col='black')
legend("topright", legend = c("PS4", "Pm", "Pl"), col = c("blue", "green", "black"), lty = 1)
#legend("topright", legend = c("stage4", "moulted", "leaving"), col = c("blue", "green", "black"), lty = 1)

###
### DAILY COUNTS PUPS PRESENT ON COL ACROSS AGE CLASSES (LEAVERS NOT INCLUDED)
### 


plot(sim_S4$whitePups,type="l",col='purple',xlab="day of season",ylab="count",main="Pups present on colony across age classes")
lines(sim_S4$S4Pups,col='orange')
lines(sim_S4$moultPups,col='green')
# Find the index of the maximum value

legend("topright", legend = c("whitecoat", "stage4", "moulted"), col = c("purple", "orange", "green"), lty = 1)
###
### DAILY COUNTS PUPS PRESENT ON COL ACROSS AGE CLASSES (LEAVERS INCLUDED)
### 

plot(sim_S4$whitePups,type="l",col='purple',xlab="day of season",ylab="count",main="Pups present on colony across age classes")
lines(sim_S4$S4Pups,col='orange')
lines(sim_S4$moultPups,col='green')
lines(sim_S4$leavePups,col='black',lty='dotted')
# Find the index of the maximum value

legend("topright", legend = c("whitecoat", "stage4", "moulted", 'leavers'), col = c("purple", "orange", "green", 'black'), lty = c(1, 1, 1, 2), cex = 0.8)


## colorblind friendly ggplot version of daily counts pups across age classses (leavers included)
library(RColorBrewer)


date_sequence <- seq(start_date, by = "day", length.out = num_days)

# Calculate the indices for the first, middle, and last dates
first_index <- 1
second_index <- 31
third_index <- 61
fourth_index <- 92
last_index <- num_days

# Extract the first, middle, and last dates
dates_to_display <- c(date_sequence[first_index], date_sequence[second_index],date_sequence[third_index],
                      date_sequence[fourth_index], date_sequence[last_index])
# Define a colorblind-friendly color palette
colorblind_palette <- brewer.pal(4, "Set2")

plot(sim_S4$whitePups, type = "l", col = colorblind_palette[1], lwd = 2,
     xlab = "Day of season", ylab = "Count", main = "Pups present on colony across age classes",
     xaxt = "n")
lines(sim_S4$S4Pups, col = colorblind_palette[2], lwd = 2)
lines(sim_S4$moultPups, col = colorblind_palette[3], lwd = 2)
lines(sim_S4$leavePups, col = colorblind_palette[4], lty = "dotted", lwd = 2)
# Add custom x-axis labels for the first, middle, and last dates
axis(side = 1, at = c(first_index, second_index, third_index, fourth_index, last_index), labels = format(dates_to_display, "%b %d"))

# Add legend with adjusted font size (cex)
legend("topright", 
       legend = c("stage 1-3", "stage 4", "moulted", "leavers"), 
       col = colorblind_palette, 
       lty = c(1, 1, 1, 2), 
       lwd = 2,
       cex = 0.9)

# Find the index of the maximum value (optional - not shown in ggplot)
# max_index <- which.max(sim_S4$whitePups + sim_S4$S4Pups + sim_S4$moultPups + sim_S4$leavePups)

# Find the index of the maximum value (optional - not shown in ggplot)
# max_index <- which.max(sim_S4$whitePups + sim_S4$S4Pups + sim_S4$moultPups + sim_S4$leavePups)

# ### remove white bc in slide, only trying to demonstrate multinomials
# plot(sim_S4$moultPups,type="l",col='green',xlab="day of season",ylab="count")#,main="Daily pup counts across age class over season")
# lines(sim_S4$S4Pups,col='blue')
# #lines(sim_S4$moultPups,col='green')
# lines(sim_S4$leavePups,col='black')
# # Find the index of the maximum value
# 
# legend("topright", legend = c( "m","S4","l"), col = c("green", "blue",'black'), lty = 1)




####################

#### Plot stochastic births 
### NOTE: births are multinomial here (see later for poisson)

plot(sim_S4$realBorn,typ="l",ylab="count",xlab="day of season",main='Daily pups born over season')
max_pup<-which.max(sim_S4$realBorn)
abline(v=max_pup,col="red")

# Plot input to births -Pb

plot(sim_S4$Pb,type="l",xlab="day of season",ylab="probability",main="Probability of birth over season")
max_index <- which.max(sim_S4$Pb)
abline(v = max_index, col = "red")



plot(sim_S4$PS4,type="l",col='blue',xlab="age",ylab="Probability density",xlim=c(0,60))#,main="Probability distribution of transitions across age classes")
lines(sim_S4$Pm,col='green')
lines(PmPl2$PlGivenm,col='black')
legend("topright", legend = c("PS4", "Pm", "Pl"), col = c("blue", "green", "black"), lty = 1)
#legend("topright", legend = c("stage4", "moulted", "leaving"), col = c("blue", "green", "black"), lty = 1)



plot(sim_S4$whitePups,type="l",col='red',xlab="day of season",ylab="count")#,main="Daily pup counts across age class over season")
lines(sim_S4$S4Pups,col='blue')
lines(sim_S4$moultPups,col='green')
lines(sim_S4$leavePups,col='black')
# Find the index of the maximum value

legend("topright", legend = c("white", "stage4", "moulted","leavers"), col = c("red", "blue", "green",'black'), lty = 1)

### remove white bc in slide, only trying to demonstrate multinomials
plot(sim_S4$moultPups,type="l",col='green',xlab="day of season",ylab="count")#,main="Daily pup counts across age class over season")
lines(sim_S4$S4Pups,col='blue')
#lines(sim_S4$moultPups,col='green')
lines(sim_S4$leavePups,col='black')
# Find the index of the maximum value

legend("topright", legend = c( "m","S4","l"), col = c("green", "blue",'black'), lty = 1)

###
### STOCHASTIC VERSION 2: Poisson births  (compare with multinomial)
###


white=0
for(t in 1:ndays) {
  scale=500
  
  white[t] <- rpois(1,scale * Pb[t]) ## Birth process (all pups born into 'white' class)
}
white
plot(white,typ='l')
max_w<-which.max(white)
lines(white,col='blue')
abline(v = max, col = "orange")



realBorn <- 0
#The t function transposes this matrix, so each row of realBorn represents a draw.
realBorn <- as.numeric(t(rmultinom(1, size = N, prob = Pb)))
plot(realBorn,typ='l')
max<-which.max(realBorn)
abline(v = max, col = "orange")

lines(white,col='blue')
abline(v = max_w, col = "orange")
