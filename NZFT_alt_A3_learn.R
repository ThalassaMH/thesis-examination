

library(popbio)
library(mc2d) # Load betaPERT distribution to describe uncertain parameters
library(sjmisc)
library(tidyr)
library(plyr)
library(reshape2) 
library(dplyr)
library(ggplot2)
library(lme4)
#library(Hmisc)
library(DHARMa)
library(MuMIn)
library(lattice)
library(sjPlot)
library(plotrix)
library(truncnorm)
library(gridExtra)


#######################################################################################
## Function to translate mean and sd of survival into parameters of a beta distribution
#######################################################################################

beta.params <- function(x_bar, sd.x){
  u <- x_bar*(1-x_bar)/sd.x^2-1
  alpha <- x_bar*u
  beta <- (1-x_bar)*u
  return(list(alpha = alpha, beta = beta))
}

#####################################################################################
## Load expert elicitation data
#####################################################################################

## excluding expert ten for colonisation round:
biol_elicit_raw<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\Fairy Terns and SDM\\NZFT\\DATA\\FAIRY-TERN-R\\df\\Biol_elicit_1to9.csv", header=T)

biol_elicit_raw[,8]<-list(NULL) #drop last blank column

## generate mean values for input into models

## filter for round 2

biol<-biol_elicit_raw %>%
  filter(round=="two") %>%
  group_by(parameter) %>%
  summarise(low.mean = mean(low, na.rm=T),best.mean = mean(best, na.rm=T), high.mean = mean(high, na.rm=T),
            low.sd = sd(low, na.rm=T), best.sd = sd(best, na.rm=T), high.sd = sd(high, na.rm=T))

biol<-as.data.frame(biol)

biol # data frame with all the values needed for the alternative strategies' models
biol$parameter<-trimws(biol$parameter)

# Function that draws the rpert value for each sim from the elicited parameter dataframe 

rperty<-function(param){
  
  df<-biol[biol$parameter == param,]
  rpert_out<-rpert(1,min = df[,2], mode = df[,3], max = df[,4]) # draw the rpert distribution, low, best, high
  df_out<-data.frame(param, rpert_out) # put the result of the rpert into a df
  return(df_out)
  
}


###################################################################
# FIELD 1 + CAPTIVE 1 (ALTERNATIVE 3) PLUS LEARNING
###################################################################

## In this scenario, Field 1 procedures are implemented plus captive solutions ##

# Captive management and manipulation:

# First clutch harvested from good pairs (50% of all breeding pairs). Harvest between 6-18 days incubation. 
# Eggs lifted, put in incubator & taken to purpose built aviary with full time staff. 
# Hatched and given access to aviary in Kaipara. Move into flight aviary with behavioural training. 
# Release as cohort at 3-9 months from flight aviary in Kaipara. 
# Supplementary food provided at aviary until it is no longer being taken. 
# Annual, for 10 years

## plus in this new scenario, productivity in captive rearing facilities is reduced by 50% for the first 3 years to show a learning process for ex situ management


#### giving us the assumptions that: 

    # juveniles are in captivity for approx 6 months after fledge
    # juveniles are released into the wild for 6 months before becoming immatures

    # "good pairs" are a minimum of 50% of the current breeding females
    #  In the fecundity equation,  multiply the number of females attempting to breed by 0.5 to give the harvest rule
    # After harvest, 52% of females will relay


#########################################################################
# Define the number of years with predictions and the Monte Carlo setting

nyears <- 50  # Number of years (projection time frame)
nsim <- 10000 # Number of replicate populations simulated
nharv <- 11   # Number of years harvesting takes place
nlearn <- 3   # Number of years of learning


# Define empty arrays, vectors and initial stage-specific population sizes
# There will now be 5 life stages: captive managed juveniles and immatures, wild juvs and imms plus one adult stage

Ni <- c(0, 0, 3, 0, 14) # starting population vector
N <- array(0, dim = c(5, nyears + 1, nsim)) # edited for five life stages (two semi-captive stages)
N[,1,] <- Ni 

alive <- matrix(0, nrow = nyears, ncol = nsim) 
r <- matrix(0, nrow = nyears, ncol = nsim)
mean.r.a3 <- lambda.a3 <- persist.a3 <- numeric()   # empty vectors for r, lambda, p(persistence > 2)
sj.val.a3 <- si.val.a3 <- sa.val.a3 <- Fa.val.a3 <- Fac.val.a3 <- Car.val.a3 <- numeric()   # empty vectors for sensitivity



# Define constants or non-elicited, fixed values

X <- 0.47                        #'cost' of management-reduction in prob of fledging from egg, calc from GLMM outputs
infer <-0.679                    # mean proportion of fertile eggs 
infer.e <- 0.156                 # uncertainty expressed as SD of the mean
c <- 1.73                        # mean clutch size (first clutch attempt)
c.e <- 0.45                      # uncertainty expressed as standard deviation of the mean
harvest <- 0.5                   # harvest rule, set to 50% of breeding females
relay <- 0.52
relay <- as.numeric(relay)
learn <- numeric()

# Define temporal variation 

Sjuv.t <- 0.001                 # temporal variability of juvenile survival (would be expressed as SD on logit scale)
Simm.t <- 0.001                 # temporal variability of immature survival (would be expressed as SD on logit scale)
Sad.t <- 0.001                  # temporal variability of immature survival (would be expressed as SD on logit scale)
h.t <- 0.2046                  # temporal variability of hatch prob (from model m2) expressed as untransformed SD
fl.t <- 0.2389                 # temporal variability of fledge prob (from model f3) expressed as untransformed SD


# Define carrying capacity, if using a fixed value

#Car <- 17      # mean 'most likely' mode value from field 1 elicitation



for (s in 1:nsim){                                          # Loop over replicate populations
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  
  p_out<-lapply(biol$parameter,rperty)
  p_out<- do.call(rbind, p_out)
  
  # Generate mean demographic values from elicited data
  
  # elicited survival parameters
  
  sj.sim <- as.numeric(filter(p_out, param == "sjuv.f1") %>% select(rpert_out))
  si.sim <- as.numeric(filter(p_out, param == "simm.f1") %>% select(rpert_out))
  sa.sim <- as.numeric(filter(p_out, param == "sad.f1") %>% select(rpert_out))
  
  sjc.sim <- as.numeric(filter(p_out, param == "sjuv.cap") %>% select(rpert_out)) # captive juvenile survival (annual)
  sjr.sim <- as.numeric(filter(p_out, param == "sjuv.rel") %>% select(rpert_out)) # released juvenile survival (annual)
  sir.sim <- as.numeric(filter(p_out, param == "simm.rel") %>% select(rpert_out)) # released immature survival (annual)
  
  
  # elicited breeding parameters
  
  h.sim <- as.numeric(filter(p_out, param == "h.f1") %>% select(rpert_out))    # wild params
  fl.sim <- as.numeric(filter(p_out, param == "fl.f1") %>% select(rpert_out))  # wild params
  
  hc.sim <- as.numeric(filter(p_out, param == "h.cap") %>% select(rpert_out))    # captive params
  flc.sim <- as.numeric(filter(p_out, param == "fl.cap") %>% select(rpert_out))  # captive params
  
  K.sim <- as.numeric(filter(p_out, param == "K.f1") %>% select(rpert_out))
  M.sim <- as.numeric(filter(p_out, param == "M.f1") %>% select(rpert_out))
  c.sim <- rtruncnorm(1, a=0, b=2, mean=c, sd=c.e)                              # not elicited for field 1
  
  car.sim <- as.numeric(filter(p_out, param=="Car.f1") %>% select(rpert_out))
  
  # Generate annual demographic rates (subject to temporal variability)
  # Bound between 0 and 1, hence the transformation through the losjrgit link
  
  sj <- plogis(rnorm(nyears, qlogis(sj.sim), Sjuv.t))
  sjc <- plogis(rnorm(nyears, qlogis(sjc.sim), Sjuv.t))
  sjr <- plogis(rnorm(nyears, qlogis(sjr.sim), Sjuv.t))
  si <- plogis(rnorm(nyears, qlogis(si.sim), Simm.t))
  sir <- plogis(rnorm(nyears, qlogis(sir.sim), Simm.t))
  sa <- plogis(rnorm(nyears, qlogis(sa.sim), Sad.t))
  ha <- plogis(rnorm(nyears, h.sim, h.t))
  fla <- plogis(rnorm(nyears, fl.sim, fl.t))
  hac <- plogis(rnorm(nyears, hc.sim, h.t))
  flac <- plogis(rnorm(nyears, flc.sim, fl.t))
  
  
  Ka <- rnorm(nyears, K.sim, 0) # no temporal var in K
  Ma <- rnorm(nyears, M.sim, 0) # no temporal var in M
  ca <- rnorm(nyears, c.sim, 0) # no temporal var in clutch size
  Car <- rnorm(nyears, car.sim, 0) # to temporal var in carrying capacity
  
  # Project population (include demographic stochasticity)
  
  for (t in 1:nyears){ # Loop over years
  
    harvest <- ifelse(t <= nharv,0.5,0) # harvests go to 0 after captive rearing finishes
    learn <- ifelse(t <= nlearn, 0.5,1) # fecundity goes up to 100% predicted values after learning is finished
    
    Fa <- (Ka[t]*ca[t]*infer*((1-harvest)+(harvest*relay))*((Ma[t]*X)+(1-Ma[t]))*(ha[t]*fla[t]))/2  # fecundity equation 
    Fac <- learn*(Ka[t]*ca[t]*infer*harvest*(hac[t]*flac[t]))/2
 
    sjc6 <- (sjc[t]^(1/12))^6 # 6 months survival in captivity
    sjr6 <- (sjr[t]^(1/12))^6 # 6 months survival in the wild after captivity
    
    N[1,t+1,s] <- rpois(1, (min(Car[t], N[5,t,s]) * sa[t] * Fa))  # wild fledglings
    N[2,t+1,s] <- rpois(1, (min(Car[t], N[5,t,s]) * sa[t] * Fac)) # captive fledglings
    N[3,t+1,s] <- rbinom(1, (N[1,t,s]), sj[t])                 # wild juveniles survival
    N[4,t+1,s] <- rbinom(1, (N[2,t,s]), (sjc6 * sjr6))         # captive + released juveniles' survival
    N[5,t+1,s] <- rbinom(1, (N[3,t,s]), si[t]) + rbinom(1, (N[4,t,s]), sir[t]) + rbinom(1, (N[5,t,s]), sa[t])
    if (sum(N[,t+1,s]) == 0) break
    else {
      r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s])) # calculate r for each time step in each sim
      alive[t,s] <- t
    } # else
    
    persist.a3[t]<- sum(N[5,t,]>2)  / s # persistance probability at each timestep
    Fac.val.a3[t] = Fac
    Car.val.a3[s] = car.sim
    
  } # t
  
  sj.val.a3[s] = sj.sim
  si.val.a3[s] = si.sim
  sa.val.a3[s] = sa.sim
  Fa.val.a3[s] = Fa
  
  
  mean.r.a3[s] <- mean(r[min(alive[,s], na.rm = TRUE):max(alive[,s], na.rm = TRUE),s])
  lambda.a3[s] <- exp(mean.r.a3[s])
  
} # s


# mean intrinsic growth rate

mean(mean.r.a3) 
sd(mean.r.a3)

mean(lambda.a3)
sd(lambda.a3)

# # store lambda 
# 
 write.csv(lambda.a3, "lambda_A3_car.csv")

#  pop summary year 50
summary(N[5,51,])

# Extinction probability (after nyears)

sum(N[5,nyears,]==0) / nsim

# calculation of quasi-extinction 

sum(apply(N[5,,1:nsim],2,min)<3)/nsim # adults only

#Probability that population is equal to or smaller than the initial adult population size in year 1, 
# and probability that population is surviving, but equal to or less than the 2017 size.

sum(N[5,nyears,]<=N[5,1,1]) / nsim

sum((N[5,nyears,]<=14) & ( N[5,nyears,]>2)) / nsim

# #store persistence
# 
write.csv(persist.a3, "persist_A3_car.csv")

# population summary
summary(N[5,51,])

########################################################
# PLOTS - FIELD 1 + CAPTIVE 1 - ALTERNATIVE 3 PUNISHED

# Compare persistence curves of SQ and f1
x<- 1:nyears

plot(x, persist, type = "l", cex = 1.5, lwd = 2, ylim = c(0.75,1), xlab = "Years", ylab = "p (persistence)")+ 
  lines(x, persist.f1, type="l", col="green", lwd = 2) + 
  lines(x, persist.f2, type="l", col = "dodgerblue", lwd = 2) +
lines(x, persist.a3, type = "l", col = "red", lwd=2)

plot(x, persist.SQf1, type = "l", cex = 1.5, lwd = 2, ylim = c(0.5,1), xlab = "Years", ylab = "p (persistence)") +
  lines(x, persist.a3, type = "l", col = "red", lwd=2)

# prep data for plots

N.mean <- apply(N,c(1,2),mean)
N.sd <- apply(N, c(1,2), sd)
N.lci <- apply(N,c(1,2),quantile, probs=c(.025))
N.uci <- apply(N,c(1,2),quantile, probs=c(.975))
N.upr <- N.mean + 1.96*N.sd                       # upper 95% confidence level (1.96 * SD)
N.lwr <- N.mean - 1.96*N.sd                       # lower 95% confidence level
N.25 <- apply(N,c(1,2),quantile, probs=c(.25))    # lower quartile
N.50 <- apply(N,c(1,2),quantile, probs=c(.5))     # median quartile
N.75 <- apply(N,c(1,2),quantile, probs=c(.75))    # upper quartile

N50_a3<-N[5,50,1:nsim] # store the adult population sizes at the end of simulation
N25_a3<-N[5,25,1:nsim] # store the adult population sizes halfway through simulation


plot_results_a3 = data.frame(alt="three",years = 1:51, N.mean=N.mean[5,], N.uci = N.uci[5,], N.lci=N.lci[5,], N.upr = N.upr[5,], N.lwr = N.lwr[5,],
                          N.25 = N.25[5,],N.50 = N.50[5,], N.75 = N.75[5,])

#store projections
write.csv(plot_results_a3, "projA3_car.csv")

#store final population size
write.csv(N50_a3, "N50A3_car.csv")


# plot population projection with mean adult population size and uncertainty over 50 years
gg2<-ggplot(plot_results_a3, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.upr, ymin=N.lwr),fill = "lightgrey", alpha=0.35) +   
  geom_ribbon(aes(x = years, ymax=N.75, ymin=N.25), fill = "dimgrey", alpha=0.35) + 
  geom_line(aes(x = years, y=N.50), color="white", size = 2, alpha=0.6) + geom_line(color="black", size=1, alpha=0.6) + 
  theme_classic(base_size = 16)+  coord_cartesian(ylim=c(0,50))  + 
  labs(title = "Predicted mean tara iti population size", x = "Year", y = "Population size (adult females)")

gg2

gg.a3<-ggplot(plot_results_a3, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.upr, ymin=N.lwr),fill = "lightgrey", alpha=0.35) + geom_line(color="black", size=1, alpha=0.6) + theme_classic(base_size = 16)+  coord_cartesian(ylim=c(0,100))  + labs(title = "Predicted mean tara iti population size", x = "Year", y = "Population size (adult females)")

gg.a3

# Plot graph with population growth rates and population size

#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(1,1))
par(las = 1, cex = 1.1)
plot(r[,1], type = "l", lwd = 0.5, ylab = "Annual population growth rate", xlab = "Time", ylim = range(r[which(!is.na(alive))]), col = "lightgrey")
for (s in 2:nsim){
  lines(r[!is.na(alive[,s]),s], lwd = 0.5, col = "lightgrey")
}
lines(apply(r, 1, mean, na.rm = TRUE), lwd = 1.5)
legend("bottomright", lwd = c(1,1), col = c("lightgrey", "black"), legend = c("Individual", "Mean"), bty = "n")
text(x = 0, y = 1.3, "a")

a <- hist(mean.r, nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate")
text(x = a$mids[1], y = max(a$counts)*0.95, "b")
box()
b<- plot(density(mean.r))

a <- hist(mean.r[not.extinct], nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate")
text(x = a$mids[1], y = max(a$counts)*0.95, "c")
box()

a <- hist(lambda.a3, nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate", xlim = (c(0.8, 1.2)))
text(x = a$mids[1], y = max(a$counts)*0.95, "b")
box()
b<- plot(density(lambda.f1))


par(mfrow=c(1,1))

plot(N[1,,s], type = "l", lwd = 0.5, ylab = "Population size", xlab = "Time", ylim = c(0,200), col = "white")
for (s in 2:nsim){
  lines(N[3,,s], lwd = 0.5, col = "dodgerblue") # plot just the adults
}


hist(N[5,50,1:nsim], xlim=c(0,300), ylim=c(0,1000), col="dodgerblue",
     breaks=seq(0,10000, 5), xlab="Final population size at t=50", main='')

lambda.density<-plot(density(lambda.sq), xlim=c(0.8, 1.1), ylim = c(0, 60), main = "", 
                     xlab = expression(lambda), lwd = 2, cex = 2) + abline(v=1, col = "grey") + 
  lines(density(lambda.f1), col = "green", lwd = 2) + lines(density(lambda.f2), col = "dodgerblue", lwd = 2) +
  lines(density(lambda.a3), col = "purple", lwd = 2) 

#################################################################################################################################

## EXPLORE THE DIFFERENT MODEL OUTPUTS FOR VITAL RATES

#################################################################################################################################




