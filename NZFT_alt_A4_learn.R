 
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

biol_elicit_raw<-read.csv("Biol_elicit_1to9.csv", header=T)

biol_elicit_raw[,8]<-list(NULL)

## generate mean values for input into models

## filter for round 2

biol<-biol_elicit_raw %>%
  filter(round=="two") %>%
  group_by(parameter) %>%
  summarise(low.mean = mean(low, na.rm=T),best.mean = mean(best, na.rm=T), high.mean = mean(high, na.rm=T),
            low.sd = sd(low, na.rm=T), best.sd = sd(best, na.rm=T), high.sd = sd(high, na.rm=T))

biol<-as.data.frame(biol)

biol # data frame with all the values needed for the alternative strategies' models

# Function that draws the rpert value for each sim from the elicited parameter dataframe 

rperty<-function(param){
  
  df<-biol[biol$parameter == param,]
  rpert_out<-rpert(1,min = df[,2], mode = df[,3], max = df[,4]) # draw the rpert distribution, low, best, high
  df_out<-data.frame(param, rpert_out) # put the result of the rpert into a df
  return(df_out)
  
}


#####################################################################################

## FIELD 1 + CAPTIVE 2 - ALTERNATIVE 4 LEARNING ##

# SAME AS ALTERNATIVE 3 BUT INFERTILE MALES COME INTO CAPTIVITY

# Capture all infertile males when appropriate and bring them into captivity to live in purpose built flight aviary

# This new model also accounts for a learning period where captive productivity is reduced by 50% for the first 3 years

#####################################################################################

######## Model assumptions ###############

## Infertile males removed. I have made the assumption that this action will impact the rate of infertility of eggs laid, and that females will 
## pair up with males with no infertility issues (big assumption!).

## Additionally, if there are no known infertile pairs, "donate" will no longer be used as a form of management. I made the simple assumption
## that proportion of nests manage would then be current rate excluding 'donation'. 

## "Donate" is the most successful form of management (it has the smallest negative impact on productivity of all the management  types).
## If I were to do this more thoroughly, I would need to calculate the 'damage' caused by management excluding 'donate' alongside the 
## rate of management excluding "donate".


#### dealing with 'infertility' in the model

# Male IDs thought to be "infertile" in current population:
# "C41570","C69024","C69027"
# Age in 2017: 15, 10 and 9 years old
# mean age 11.333

# These males are in the current population and so to calculate the change in fertility rate of eggs, I only looked at the last decade (2007-2017)
# (given their mean age is 11.3)
# Code for calculation is in "breed_explore_22_02.R"

## Mean rate of fertility once breeding attempts with these males are removed from the calculation:
### 0.90 (+/- 0.08 SD) 

#### dealing with "donate" management type in the model

# remove 'donate' management type and re-calculate rate of management for the last decade (1997-2017)
## mean annual management rate (calculated as total proportion of nests managed per season)
### 0.44 (+/- 0.21 SD, 0.07 SE)

## ISSUE TO FIX 18.11 - RESOLVED?

# this doesn't account for the lowering of management rate due to field actions
# so should I calculate the elicted proportion reduction in management rate and then also lower the management
# by ten percent as per the calculation above?

## Potential solution - just use M.f1 for all alternatives
# There is not a lot of difference between M.f1 and M.a5. 

# It standardises everything across the alternatives to say that the basic new interventions of field 1 will
# reduce the need to manage clutches by ~20%

########################################################################



#########################################################################
# Define the number of years with predictions and the Monte Carlo setting

nyears <- 50  # Number of years (projection time frame)
nsim <- 10000 # Number of replicate populations simulated
nharv <- 11   # Number of years harvesting takes place
nlearn <- 3 # number of years of learning captive work

# Define empty arrays, vectors and initial stage-specific population sizes
# There will now be 5 life stages: captive managed juveniles and immatures, wild juvs and imms plus one adult stage

Ni <- c(0, 0, 3, 0, 14) # starting population vector
N <- array(0, dim = c(5, nyears + 1, nsim)) # edited for five life stages (two semi-captive stages)
N[,1,] <- Ni 

alive <- matrix(0, nrow = nyears, ncol = nsim) 
r <- matrix(0, nrow = nyears, ncol = nsim)
mean.r.a4 <- lambda.a4 <- persist.a4 <- numeric()   # empty vectors for r, lambda, p(persistence > 2 adult females)
sj.val.a4 <- si.val.a4 <- sa.val.a4 <- Fa.val.a4 <- Fac.val.a4 <- Car.val.a4 <- numeric()   # empty vectors for sensitivity



# Define constants or non-elicited, fixed values

X <- 0.47                        #'cost' of management-reduction in prob of fledging from egg, calc from GLMM outputs
infer <-0.9                    # mean proportion of fertile eggs 
infer.e <- 0.08                 # uncertainty expressed as SD of the mean
c <- 1.73                        # mean clutch size (first clutch attempt)
c.e <- 0.45                      # uncertainty expressed as standard deviation of the mean
# 
# M <-0.44                       # mean proportion of eggs managed (2007-2017) excluding donation
# M.e <- 0.21                      # uncertainty expressed as SD of the mean

harvest <- 0.5                   # harvest rule, set to 50% of breeding females
relay <- 0.52
relay <- as.numeric(relay)
learn <- numeric()

# Define temporal varition 

Sjuv.t <- 0.001                 # temporal variability of juvenile survival (would be expressed as SD on logit scale)
Simm.t <- 0.001                 # temporal variability of immature survival (would be expressed as SD on logit scale)
Sad.t <- 0.001                  # temporal variability of immature survival (would be expressed as SD on logit scale)
h.t <- 0.2046                  # temporal variability of hatch prob (from model m2) expressed as untransformed SD
fl.t <- 0.2389                 # temporal variability of fledge prob (from model f3) expressed as untransformed SD

# carrying capacity

#Car <- 17     # most likely value from field 1 elicitation
#Car<-100
#Car <- 37      # most likely by optimistic expert 2

## model

for (s in 1:nsim){                                          # Loop over replicate populations
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  
  p_out<-lapply(biol$parameter,rperty)
  p_out<- do.call(rbind, p_out)
  
  # Generate mean demographic values from elicited data
  
  # elicited survival parameters
  
  sj.sim <- as.numeric(filter(p_out, param == "sjuv.f1") %>% select(rpert_out)) # dplyr
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
  # M.sim <-  rbeta(1, beta.params(M, M.e)$alpha, beta.params(M, M.e)$beta)
  M.sim <- as.numeric(filter(p_out, param == "M.f1") %>% select(rpert_out))
  K.sim <- as.numeric(filter(p_out, param == "K.f1") %>% select(rpert_out))
  c.sim <- rtruncnorm(1, a=0, b=2, mean=c, sd=c.e)                              # not elicited for field 1
  
  # ELICITED carrying capacity
  
  car.sim <- as.numeric(filter(p_out, param=="Car.f1") %>% select(rpert_out))
  
  # Generate annual demographic rates (subject to temporal variability)
  # Bound between 0 and 1, hence the transformation through the logit link
  
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
 
  Ma <- rnorm(nyears, M.sim, 0) # no temporal var in M
  Ka <- rnorm(nyears, K.sim, 0) # no temporal var in K
  ca <- rnorm(nyears, c.sim, 0) # no temporal var in clutch size
  Car <- rnorm(nyears, car.sim, 0) # no temporal var in carrying capacity
  
  # Project population (include demographic stochasticity)
  
  for (t in 1:nyears){ # Loop over years
    
    harvest <- ifelse(t <= nharv,0.5,0)
    learn <- ifelse(t <= nlearn, 0.5,1)
    
    Fa <- (Ka[t]*ca[t]*infer*((1-harvest)+(harvest*relay))*((Ma[t]*X)+(1-Ma[t]))*(ha[t]*fla[t]))/2  # fecundity equation 
    Fac <- learn*(Ka[t]*ca[t]*infer*harvest*(hac[t]*flac[t]))/2
    

    sjc6 <- (sjc[t]^(1/12))^6 # 6 months survival in captivity
    sjr6 <- (sjr[t]^(1/12))^6 # 6 months survival in the wild after captivity
    
    N[1,t+1,s] <- rpois(1, min(Car[t], N[5,t,s]) * sa[t] * Fa)  # wild fledglings
    N[2,t+1,s] <- rpois(1, min(Car[t], N[5,t,s]) * sa[t] * Fac) # # captive fledglings
    N[3,t+1,s] <- rbinom(1, (N[1,t,s]), sj[t])       # wild juveniles survival
    N[4,t+1,s] <- rbinom(1, (N[2,t,s]), (sjc6 * sjr6)) # captive + released juveniles' survival
    N[5,t+1,s] <- rbinom(1, (N[3,t,s]), si[t]) + rbinom(1, (N[4,t,s]), sir[t]) + rbinom(1, (N[5,t,s]), sa[t])
    if (sum(N[,t+1,s]) == 0) break
    else {
      r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s])) # calculate r for each time step in each sim
      alive[t,s] <- t
    } # else
    
    persist.a4[t]<- sum(N[5,t,]>2)  / s # persistance probability at each timestep
    Fac.val.a4[t] = Fac
    
  } # t
  
  sj.val.a4[s] = sj.sim
  si.val.a4[s] = si.sim
  sa.val.a4[s] = sa.sim
  Fa.val.a4[s] = Fa
  Car.val.a4[s] = car.sim
  
  mean.r.a4[s] <- mean(r[min(alive[,s], na.rm = TRUE):max(alive[,s], na.rm = TRUE),s])
  lambda.a4[s] <- exp(mean.r.a4[s])
  
} # s


# mean intrinsic growth rate

mean(mean.r.a4) 
sd(mean.r.a4)

mean(lambda.a4)
sd(lambda.a4)

## look at productivity rates to check the code is working
Fa.val.a4[,1]
N[2,,8]
N[4,,10]

# carrying capacity
plot(hist(Car.val.a4))

# Extinction probability (after nyears)

sum(N[5,nyears,]==0) / nsim

# calculation of quasi-extinction 

sum(apply(N[5,,1:nsim],2,min)<3)/nsim # adults only

# work out how to remove top two lines of matrix?
#sum(apply(colSums(N[-1,,1:nsim]),2,min)<3)/nsim # immatures and adults only

#Probability that population is equal to or smaller than the initial adult population size in year 1, 
#and probability that population is surviving, but equal to or less than the 2017 size.

sum(N[5,nyears,]<=N[5,1,1]) / nsim

sum((N[5,nyears,]<=14) & ( N[5,nyears,]>2)) / nsim

## store persist and lambda

write.csv(persist.a4, "persist_A4_car.csv")
write.csv(lambda.a4, "lambda_A4_car.csv")

summary(N[5,51,])

########################################################
# PLOTS - FIELD 1 + CAPTIVE 2 - ALTERNATIVE 4 PLUS LEARNING

# Compare persistence curves of SQ and f1
x<- 1:nyears
plot(x, persist, type="l", ylim=c(0.5, 1), xlab = "Years", ylab = "Probability of persistence")
lines(x, persist.f1, type="l", col = "green")
lines(x, persist.f2, type = "l", col = "dodgerblue")
lines(x, persist.a3, type = "l", col = "red")
lines(x, persist.a4, type = "l", col = "purple")

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

N50_a4<-N[5,50,1:nsim] # store the adult population sizes at the end of simulation
N25_a4<-N[5,25,1:nsim] # store the adult population sizes halfway through simulation


plot_results_a4_learn = data.frame(alt = "four", years = 1:51, N.mean=N.mean[5,], N.uci = N.uci[5,], N.lci=N.lci[5,], N.upr = N.upr[5,], N.lwr = N.lwr[5,],
                          N.25 = N.25[5,],N.50 = N.50[5,], N.75 = N.75[5,])

#store projections
write.csv(plot_results_a4_learn, "projA4_car.csv")

#store final population size
write.csv(N50_a4, "N50A4_car.csv")

# plot population projection with mean adult population size and uncertainty over 50 years
gg2<-ggplot(plot_results, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.upr, ymin=N.lwr),fill = "lightgrey", alpha=0.35) +   
  geom_ribbon(aes(x = years, ymax=N.75, ymin=N.25), fill = "dimgrey", alpha=0.35) + 
  geom_line(aes(x = years, y=N.50), color="white", size = 2, alpha=0.6) + geom_line(color="black", size=1, alpha=0.6) + 
  theme_classic(base_size = 16)+  coord_cartesian(ylim=c(0,100))  + 
  labs(title = "Predicted mean tara iti population size", x = "Year", y = "Population size (adult females)")

gg2

gg.a4<-ggplot(plot_results_a4, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.upr, ymin=N.lwr),fill = "lightgrey", alpha=0.35) + geom_line(color="black", size=1, alpha=0.6) + theme_classic(base_size = 16)+  coord_cartesian(ylim=c(0,250))  + labs(title = "Predicted mean tara iti population size", x = "Year", y = "Population size (adult females)")

gg.a4

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

a <- hist(mean.r.a4, nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate")
text(x = a$mids[1], y = max(a$counts)*0.95, "b")
box()
b<- plot(density(mean.r))

a <- hist(mean.r.a4[not.extinct], nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate")
text(x = a$mids[1], y = max(a$counts)*0.95, "c")
box()

a <- hist(lambda.a4, nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate", xlim = (c(0.8, 1.2)))
text(x = a$mids[1], y = max(a$counts)*0.95, "b")
box()
b<- plot(density(lambda.a5))


par(mfrow=c(1,1))

plot(N[5,,s], type = "l", lwd = 0.5, ylab = "Population size", xlab = "Time", ylim = c(0,200), col = "white")
for (s in 2:nsim){
  lines(N[5,,s], lwd = 0.5, col = "dodgerblue") # plot just the adults
}


hist(N[5,50,1:nsim], xlim=c(0,300), ylim=c(0,1000), col="dodgerblue",
     breaks=seq(0,10000, 5), xlab="Final population size at t=50", main='')

lambda.density<-plot(density(lambda.sq), xlim=c(0.9, 1.1), ylim = c(0, 60), main = "", 
                     xlab = expression(lambda), lwd = 2, cex = 2) + abline(v=1, col = "grey") + 
  lines(density(lambda.f1), col = "green", lwd = 2) + lines(density(lambda.f2), col = "dodgerblue", lwd = 2) +
  lines(density(lambda.a3), col = "red", lwd = 2) + lines(density(lambda.a4), col = "purple", lwd = 2) +  lines(density(lambda.a5), col = "orange", lwd = 2)

#################################################################################################################################

## EXPLORE THE DIFFERENT MODEL OUTPUTS FOR VITAL RATES

#################################################################################################################################

# 
# corr_df<-data.frame(sj<-sj.val, si<-si.val, sa<- sa.val, Fa<-Fa.val, r<-mean.r,
#                     sj1<-sj.val.f1, si1<-si.val.f1, sa1<- sa.val.f1, Fa1<-Fa.val.f1, r1<-mean.r.f1,
#                     sj2<-sj.val.f2, si2<-si.val.f2, sa2<- sa.val.f2, Fa2<-Fa.val.f2, r2<-mean.r.f2,
#                     sj3<-sj.val.a3, si3<-si.val.a3, sa3<- sa.val.a3, Fa3<-Fa.val.a3, Fac3<-Fac.val.a3, r3<-mean.r.a3)

cor5<-ggplot(corr_df, aes(Fa3, r3)) 
+ geom_point(shape = 21, size = 1.5) 
+ labs(title = "Correlation of juvenile survival + mean population growth rate")
+ theme_classic(base_size = 12) 
+ labs(title = "Correlation of fecundity + mean population growth rate") 
+ theme_classic(base_size = 12) 

cor5+ geom_point(Fa1, r1, color = "#207394", shape = 21, size = 1.5)

corr_df$Fac3 # issue - for some reason, the Fac values are not being stored! Always showing 0.




