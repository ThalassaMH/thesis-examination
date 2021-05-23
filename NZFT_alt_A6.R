
library(popbio)
library(mc2d) # Load betaPERT distribution to describe uncertain parameters
library(sjmisc)
library(tidyr)
library(plyr)
library(reshape2) 
library(dplyr)
library(ggplot2)
library(lme4)
library(Hmisc)
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

# Function that draws the rpert value for each sim from the elicited parameter dataframe (Fi Spooner)

rperty<-function(param){
  
  df<-biol[biol$parameter == param,]
  rpert_out<-rpert(1,min = df[,2], mode = df[,3], max = df[,4]) # draw the rpert distribution, low, best, high
  df_out<-data.frame(param, rpert_out) # put the result of the rpert into a df
  return(df_out)
  
}


#####################################################################################
# Define the number of years with predictions and the Monte Carlo setting

nyears <- 50 # Number of years (projection time frame)
nsim <- 100 # Number of replicate populations simulated
nsupp <- 11 #Number of years nests are supplemented with OZFT


# Define empty arrays, vectors and initial stage-specific population sizes

Ni <- c(0,0,3,14) # starting population vector
N <- array(0, dim = c(4, nyears + 1, nsim)) # edited for four life stages
N[,1,] <- Ni 

alive <- matrix(0, nrow = nyears, ncol = nsim) 
r <- matrix(0, nrow = nyears, ncol = nsim)
mean.r.a6 <- numeric()
sj.val.a6 <- si.val.a6 <- sa.val.a6 <- Fa.val.a6 <- persist.a6<- lambda.a6<-Faz.val.a6<-sjz.val<-Car.val.a6<- numeric() 


## s is parametric uncertainty

# Define constants or non-elicited, fixed values

X <- 0.47                        #'cost' of management-reduction in prob of fledging from egg, calc from GLMM outputs
infer <-0.679                    # mean proportion of fertile eggs 
infer.e <- 0.156                 # uncertainty expressed as SD of the mean


# Define temporal varition 

Sjuv.t <- 0.01                 # temporal variability of juvenile survival (would be expressed as SD on logit scale)
Simm.t <- 0.01                 # temporal variability of immature survival (would be expressed as SD on logit scale)
Sad.t <- 0.01                  # temporal variability of immature survival (would be expressed as SD on logit scale)
h.t <- 0.2046                  # temporal variability of hatch prob (from model m2) expressed as untransformed SD
fl.t <- 0.2389                 # temporal variability of fledge prob (from model f3) expressed as untransformed SD

# Define carrying capacity, if using a fixed value

#Car <- 18        # most likely from elicitation

# Run the projection:
  

# Project population

for (s in 1:nsim){                                          # Loop over replicate populations
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  
  p_out<-lapply(biol$parameter,rperty)
  p_out<- do.call(rbind, p_out)
  
  # Generate mean demographic values from elicited data
  
  # elicited survival parameters
  
  sj.sim <- as.numeric(filter(p_out, param == "sjuv.f2") %>% select(rpert_out)) 
  si.sim <- as.numeric(filter(p_out, param == "simm.f2") %>% select(rpert_out))
  sa.sim <- as.numeric(filter(p_out, param == "sad.f2") %>% select(rpert_out))
  sjz.sim<- as.numeric(filter(p_out, param== "sjuv.oz") %>% select(rpert_out)) # ozft juvenile surv
  
  # elicited breeding parameters
  
  h.sim <- as.numeric(filter(p_out, param == "h.f2") %>% select(rpert_out))
  fl.sim <- as.numeric(filter(p_out, param == "fl.f2") %>% select(rpert_out))
  hz.sim <- as.numeric(filter(p_out, param == "h.oz2") %>% select(rpert_out)) # oz egg hatch prob
  flz.sim <- as.numeric(filter(p_out, param == "fl.oz2") %>% select(rpert_out)) # oz egg fledge prob
  
  M.sim <- as.numeric(filter(p_out, param == "M.f1") %>% select(rpert_out)) 
  K.sim <- as.numeric(filter(p_out, param == "K.f2") %>% select(rpert_out))
  c.sim <- as.numeric(filter(p_out, param == "c.f2") %>% select(rpert_out))
  
  car.sim <- as.numeric(filter(p_out, param=="Car.f2") %>% select(rpert_out))
  
  # Generate annual demographic rates (subject to temporal variability)
  
  sj <- plogis(rnorm(nyears, qlogis(sj.sim), Sjuv.t))
  si <- plogis(rnorm(nyears, qlogis(si.sim), Simm.t))
  sa <- plogis(rnorm(nyears, qlogis(sa.sim), Sad.t))
  sjz <- plogis(rnorm(nyears, qlogis(sjz.sim), Sjuv.t))
  ha <- plogis(rnorm(nyears, h.sim, h.t))
  fla <- plogis(rnorm(nyears, fl.sim, fl.t))
  haz <- plogis(rnorm(nyears, hz.sim, h.t))
  flaz <- plogis(rnorm(nyears, flz.sim, fl.t))
  
  Ma <- rnorm(nyears, M.sim, 0) # no temporal var in M
  Ka <- rnorm(nyears, K.sim, 0) # no temporal var in prop of fem breeding
  ca <- rnorm(nyears, c.sim, 0) # no temporal var in clutch size
  Car <- rnorm(nyears, car.sim, 0) # to temporal var in carrying capacity
  
  # Project population (include demographic stochasticity)
  
  for (t in 1:nyears){ # Loop over years
    supp <- ifelse(t <= nsupp,1,0) # when time gets above the ten year period, supplementing turns to zero
    
    Fa <- (Ka[t]*ca[t]*infer*((Ma[t]*X)+(1-Ma[t]))*(ha[t]*fla[t]))/2  # fecundity equation 
    Faz <- supp*(Ka[t]*(2-(ca[t]*infer))*((Ma[t]*X)+(1-Ma[t]))*(haz[t]*flaz[t]))/2 # fecundity of OZFT 
    
    N[1,t+1,s] <- rpois(1, (min(Car[t], N[4,t,s]) * sa[t] * Fa))
    N[2,t+1,s] <- rpois(1, (min(Car[t], N[4,t,s]) * sa[t] * Faz))
    N[3,t+1,s] <- rbinom(1, (N[1,t,s]), sj[t]) + rbinom(1, (N[2,t,s]), sjz[t]) # surv of NZ and OZ juvs
    N[4,t+1,s] <- rbinom(1, (N[3,t,s]), si[t]) + rbinom(1, (N[4,t,s]), sa[t]) # surv of NZ and OZ the same from immature onwards
    
    if (sum(N[,t+1,s]) == 0) break
    else {
      r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s])) # calculate r for each time step in each sim
      alive[t,s] <- t
    } # else
    
    persist.a6[t]<- sum(N[4,t,]>2)  / s # persistance probability at each timestep
    
  } # t
  
  sj.val.a6[s] = sj.sim
  sjz.val[s] = sjz.sim
  si.val.a6[s] = si.sim
  sa.val.a6[s] = sa.sim
  Fa.val.a6[s] = Fa
  Faz.val.a6[s] = Faz
  Car.val.a6[s] = car.sim
  
  mean.r.a6[s] <- mean(r[min(alive[,s], na.rm = TRUE):max(alive[,s], na.rm = TRUE),s])
  lambda.a6[s] <- exp(mean.r.a6[s])
  
} # s


# mean intrinsic growth rate

mean(mean.r.a6) 
sd(mean.r.a6)

mean(lambda.a6)
sd(lambda.a6)

# store lambda

write.csv(lambda.a6, "lambda_A6_car.csv")

# Extinction probability (after nyears)

sum(N[4,nyears,]==0) / nsim

# calculation of quasi-extinction 

sum(apply(N[4,,1:nsim],2,min)<3)/nsim # adults only

# store persistence 

write.csv(persist.a6, "persist_A6_car.csv")

#Probability that population is equal to or smaller than the initial adult population size in year 1, 

sum(N[4,nyears,]<=N[4,1,1]) / nsim

summary(N[4,51,])

########################################################
# PLOTS - FIELD 2 + OZFT  -ALTERNATIVE 6

# Compare persistence curves of SQ and f1
x<- 1:nyears
plot(x, persist, type="l", ylim=c(0.5, 1), xlab = "Years", ylab = "Probability of persistence", main="Alternative persistance probabilities")
lines(x, persist.f1, type="l", col = "green")
lines(x, persist.f2, type = "l", col = "dodgerblue")
lines(x, persist.a3, type = "l", col = "red")
lines(x, persist.a4, type = "l", col = "purple")
lines(x, persist.a5, type = "l", col = "orange")
lines(x, persist.a6, type = "l", col = "brown")
legend(1, 0.8, legend=c("Status Quo", "A1", "A2", "A3", "A4", "A5", "A6"), 
       col=c("black", "green", "dodgerblue", "red","purple", "orange","brown"), lty=1, cex=0.8)


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

N50_a6<-N[4,50,1:nsim] # store the adult population sizes at the end of simulation
N25_a6<-N[4,25,1:nsim] # store the adult population sizes halfway through simulation

plot_results_a6 = data.frame(alt = "six", years = 1:51, N.mean=N.mean[4,], N.uci = N.uci[4,], N.lci=N.lci[4,], N.upr = N.upr[4,], N.lwr = N.lwr[4,],
                              N.25 = N.25[4,],N.50 = N.50[4,], N.75 = N.75[4,])

# store population size and projections

write.csv(N50_a6, "N50A6_car.csv")
write.csv(plot_results_a6, "projA6_car.csv")
 
# plot population projection with mean adult population size and uncertainty over 50 years
gg2<-ggplot(plot_results_a6, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.uci, ymin=N.lci),fill = "lightgrey", alpha=0.35) +   
  geom_ribbon(aes(x = years, ymax=N.75, ymin=N.25), fill = "dimgrey", alpha=0.35) + 
  geom_line(aes(x = years, y=N.50), color="white", size = 2, alpha=0.6) + geom_line(color="black", size=1, alpha=0.6) + 
  theme_classic(base_size = 16)+  coord_cartesian(ylim=c(0,100))  + 
  labs(title = "Predicted mean tara iti population size", x = "Year", y = "Population size (adult females)")

gg2

gg.a6<-ggplot(plot_results_a6.1, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.upr, ymin=N.lwr),fill = "lightgrey", alpha=0.35) + 
  geom_line(color="black", size=1, alpha=0.6) + theme_classic(base_size = 16)+  
  coord_cartesian(ylim=c(0,100))  + 
  labs(title = "Predicted mean tara iti population size", x = "Year", y = "Population size (adult females)")

gg.a6

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

a <- hist(mean.r.a6, nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate")
text(x = a$mids[1], y = max(a$counts)*0.95, "b")
box()
b<- plot(density(mean.r.a6))

a <- hist(mean.r.a6[not.extinct], nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate")
text(x = a$mids[1], y = max(a$counts)*0.95, "c")
box()

a <- hist(lambda.a6, nclass = 25, col = "dodgerblue", main = "", xlab = "Population growth rate", xlim = (c(0.8, 1.2)))
text(x = a$mids[1], y = max(a$counts)*0.95, "b")
box()
b<- plot(density(lambda.a6))


par(mfrow=c(1,1))

plot(N[4,,s], type = "l", lwd = 0.5, ylab = "Population size", xlab = "Time", ylim = c(0,200), col = "white")
for (s in 2:nsim){
  lines(N[4,,s], lwd = 0.5, col = "dodgerblue") # plot just the adults
}


hist(N[4,50,1:nsim], xlim=c(0,300), ylim=c(0,1000), col="dodgerblue",
     breaks=seq(0,10000, 5), xlab="Final population size at t=50", main='')

lambda.density<-plot(density(lambda.sq), xlim=c(0.8, 1.2), ylim = c(0, 50), main = "", 
                     xlab = expression(lambda), lwd = 2, cex = 4) + abline(v=1, col = "grey") + 
  lines(density(lambda.f1), col = "green", lwd = 2) + lines(density(lambda.f2), col = "dodgerblue", lwd = 2) +
  lines(density(lambda.a3), col = "red", lwd = 2) + lines(density(lambda.a4), col = "purple", lwd = 2) +
 lines(density(lambda.a5), col = "orange", lwd = 2) +  lines(density(lambda.a6), col = "brown", lwd = 2)
legend(0.8, 50, legend = c("SQ", "Field 1","Field 2", "Field 1 + Captive 1", "Field 1 + Captive 2", "Field 2 + Captive 3", "Field 2 + OZFT"), 
       col = c("black", "green", "dodgerblue", "red", "purple", "orange", "brown"), lty = 1, lwd = 2, cex = 1)

#####################################################################################################

#################################################################################################################################

## EXPLORE THE DIFFERENT MODEL OUTPUTS FOR VITAL RATES

#################################################################################################################################


corr_df<-data.frame(sj<-sj.val, si<-si.val, sa<- sa.val, Fa<-Fa.val, r<-mean.r,
                    sj6<-sj.val.a6, si6<-si.val.a6, sa6<- sa.val.a6, Fa6<-Fa.val.a6, 
                   r6<-mean.r.a6)

cor5<-ggplot(corr_df, aes(Fa6, r6)) + geom_point(shape = 21, size = 1.5) + labs(title = "Correlation of juvenile survival + mean population growth rate") +   theme_classic(base_size = 12) + labs(title = "Correlation of fecundity + mean population growth rate") +
  theme_classic(base_size = 12) 

cor5+ geom_point(Fa6, r6, color = "#207394", shape = 21, size = 1.5)

corr_df$Fac3 # issue - for some reason, the Fac values are not being stored! Always showing 0.


#########################################################################################################################

## plot all the population size projections into one graph #

plot_all_results<-rbind(plot_results, plot_results_a1, plot_results_a2, plot_results_a3, plot_results_a4,
                        plot_results_a5, plot_results_a6, plot_results_a8)
head(plot_all_results)
write.csv(plot_all_results, "N_projections_all_alts_elicitedK.csv")
plot_all_results<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\Fairy Terns & SDM\\Decision Making\\TARA ITI WORKSHOP\\Phase 4\\Population modelling results\\N_projections_all_alts_elicitedK.csv", header=T)

# project_all<-ggplot(plot_all_results, aes(x=years, y=N.mean)) + 
#    geom_line(aes(color=alt), size=1, alpha=0.6) + theme_classic(base_size = 16)+  
#   coord_cartesian(ylim=c(0,150))  + theme(legend.position = "right") + geom_ribbon(aes(ymax=N.upr, ymin=N.lwr, fill = alt), alpha=0.2)+
#   labs(title = "Predicted mean tara iti population size (max. territories 100)", x = "Year", y = "Population size (adult females)") 
# project_all

project_all<-ggplot(plot_all_results, aes(x=years, y=N.mean)) + 
  geom_line(aes(color=alt), size=1.5) + theme_classic(base_size = 16)+  
  coord_cartesian(ylim=c(0,75))  + theme(legend.position = "right") + geom_ribbon(aes(ymax=N.upr, ymin=N.lwr,fill = alt),  alpha=0.1) +
  labs(title = "Predicted mean tara iti population size (elicited territories)", x = "Year", y = "Population size (adult females)") 
project_all

plot(density(N[3,50,1:nsim]), xlab="Final population size at t=50", main='', lwd = 2)
hist(N[3,50,1:nsim], xlim=c(0,300), ylim=c(0,1000), col="dodgerblue",
     breaks=seq(0,10000, 5), xlab="Final population size at t=50", main='')



