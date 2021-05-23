####################################################################################

### STATUS QUO PRODUCTIVITY WITH FIELD 1 SURVIVAL ESTIMATES ALL AGE CLASSES


# replace status quo parameters with field 1 elicited estimates:

          #   ad surv, juv surv, imm surv


#####################################################################################
## only run this bit if needed

# Function to translate mean and sd of survival into parameters of a beta distribution

beta.params <- function(x_bar, sd.x){
  u <- x_bar*(1-x_bar)/sd.x^2-1
  alpha <- x_bar*u
  beta <- (1-x_bar)*u
  return(list(alpha = alpha, beta = beta))
}

## Load expert elicitation data

biol_elicit_raw<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\Fairy Terns and SDM\\NZFT\\DATA\\FAIRY-TERN-R\\df\\Biol_elicit_1to9.csv", header=T)
biol_elicit_raw[,8]<-list(NULL) # removing an odd column

## generate mean values for input into models

## filter for round 2

biol<-biol_elicit_raw %>%
  filter(round=="two") %>%
  group_by(parameter) %>%
  summarise(low.mean = mean(low, na.rm=T),best.mean = mean(best, na.rm=T), high.mean = mean(high, na.rm=T),
            low.sd = sd(low, na.rm=T), best.sd = sd(best, na.rm=T), high.sd = sd(high, na.rm=T))

biol<-as.data.frame(biol) # data frame with all the values needed for the alternative strategies' models

# Function that draws the rpert value for each sim from the elicited parameter dataframe

rperty<-function(param){
  
  df<-biol[biol$parameter == param,]
  rpert_out<-rpert(1,min = df[,2], mode = df[,3], max = df[,4]) # draw the rpert distribution, low, best, high
  df_out<-data.frame(param, rpert_out) # put the result of the rpert into a df
  return(df_out)
}


######################################################################################################

# Define the number of years with predictions and the Monte Carlo setting

nyears <- 50 # Number of years (projection time frame)
nsim <- 10000 # Number of replicate populations simulated


# Define population matrix and initial stage-specific population sizes

Ni <- c(0, 3, 14) # starting population vector
N <- array(0, dim = c(3, nyears + 1, nsim)) # edited for three life stages
N[,1,] <- Ni 

alive <- matrix(0, nrow = nyears, ncol = nsim) 
r <- matrix(0, nrow = nyears, ncol = nsim)
mean.r.SQf1 <- lambda.SQf1 <- persist.SQf1 <- numeric()   # empty vectors for r, lambda, p(persistence > 2)
sj.val.SQf1 <- si.val.SQf1 <- sa.val.SQf1 <- Fa.val.SQf1 <-car.val.SQf1<- numeric()   # empty vectors for sensitivity


# Define mean, measurement error and temporal variability of the demographic parameters
# NB temporal variability for survival has been generically assigned as 0.01 - this is *not* based on real data


Sjuv <- 0.8053
Sjuv.e <- 0.097                # Uncertainty of juv. survival expressed as std error from RMark output
Sjuv.t <- 0.001                # temporal variability of juvenile survival (would be expressed as SD on logit scale)
Simm <- 0.927
Simm.e <- 0.061                # Uncertainty of imm survival expressed as std error from RMark output
Simm.t <- 0.001                # temporal variability of immature survival (would be expressed as SD on logit scale)
Sad <- 0.919
Sad.e <- 0.023                 # Uncertainty of imm survival expressed as std error from RMark output
Sad.t <- 0.001                 # temporal variability of immature survival (would be expressed as SD on logit scale)

h <- qlogis(0.809)              # untransformed prob of eggs hatching (backtransformed from ggpredict output)
h.e <- 0.346                    # uncertainty expressed as untransformed std error from ggpredict 
h.t <- 0.2046                   # temporal variability of hatch prob (from model m2) expressed as untransformed SD in the random effects for season start var
fl <-qlogis(0.704)              # untransformed prob of hatchlings fledging (backtransformed from ggpredict output)
fl.e <-0.410                    # uncertainty expressed as untransformed std error from ggpredict 
fl.t <- 0.2389                  # temporal variability of fledge prob (from model f3) expressed as untransformed SD
K <-0.665                       # mean prop of females breeding (1997-2017)
K.e <-0.18                      # uncertainty expressed as standard deviation of the mean
c <- 1.73                      # mean clutch size 
c.e <- 0.45                    # uncertainty expressed as standard deviation of the mean

## constants 

X <- 0.47                        #'cost' of management-reduction in prob of fledging from egg
infer <-0.679                    # mean proportion of fertile eggs 
infer.e <- 0.156                 # uncertainty expressed as SD of the mean
M <-0.542                        # mean proportion of eggs managed (2007-2017)
M.e <- 0.21                      # uncertainty expressed as SD of the mean

#Car <- 17
# set a carrying capacity for number of territories

for (s in 1:nsim){                                          # Loop over replicate populations
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  
  p_out<-lapply(biol$parameter,rperty)
  p_out<- do.call(rbind, p_out)
  
  # Generate mean demographic values from data
  
  # ELICITED survival parameters
  
  sj.sim <- as.numeric(filter(p_out, param == "sjuv.f1") %>% select(rpert_out)) # FIELD 1 SURV
  si.sim <- as.numeric(filter(p_out, param == "simm.f1") %>% select(rpert_out)) # FIELD 1 SURV
  sa.sim <- as.numeric(filter(p_out, param == "sad.f1") %>% select(rpert_out)) # FIELD 1 SURV
  
  # ELICITED carrying capacity
  car.sim <- as.numeric(filter(p_out, param=="Car.f1") %>% select(rpert_out))
  
  # status quo, data driven breeding parameters
  
  h.sim <-  rnorm(1, h, h.e)             # transformed in the next param simulation 
  fl.sim <- rnorm(1, fl, fl.e)           # transformed in the next param simulation 
  K.sim <- rbeta(1, beta.params(K, K.e)$alpha, beta.params(K, K.e)$beta)
  c.sim <- rtruncnorm(1, a=0, b=2, mean=c, sd=c.e) 
  M.sim <-  rbeta(1, beta.params(M, M.e)$alpha, beta.params(M, M.e)$beta) 
  
 
  
  # Generate annual demographic rates (subject to temporal variability)
  
  sj <- plogis(rnorm(nyears, qlogis(sj.sim), Sjuv.t))
  si <- plogis(rnorm(nyears, qlogis(si.sim), Simm.t))
  sa <- plogis(rnorm(nyears, qlogis(sa.sim), Sad.t))
  ha <- plogis(rnorm(nyears, h.sim, h.t))
  fla <- plogis(rnorm(nyears, fl.sim, fl.t))
  Ka <- rnorm(nyears, K.sim, 0) # no temporal var in K
  Ma <- rnorm(nyears, M.sim, 0) # no temporal var in M
  ca <- rnorm(nyears, c.sim, 0) # no temporal var in clutch size
  Car <- rnorm(nyears, car.sim, 0) # to temporal var in carrying capacity
  
  # Project population (include demographic stochasticity)
  
  
  for (t in 1:nyears){ # Loop over years
    
    Fa <- (Ka[t]*ca[t]*infer*((Ma[t]*X)+(1-Ma[t]))*(ha[t]*fla[t]))/2 # fecundity equation
    
    N[1,t+1,s] <- rpois(1, min(Car[t], N[3,t,s]) * sa[t] * Fa)
    N[2,t+1,s] <- rbinom(1, (N[1,t,s]), sj[t])
    N[3,t+1,s] <- rbinom(1, (N[2,t,s]), si[t]) + rbinom(1, (N[3,t,s]), sa[t])
    if (sum(N[,t+1,s]) == 0) break
    else {
      r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s])) # calculate r for each time step in each sim
      alive[t,s] <- t
    } # else
    
    persist.SQf1[t]<- sum(N[3,t,]>2)  / s # persistance probability at each timestep
    
  } # t
  
  sj.val.SQf1[s] = sj.sim
  si.val.SQf1[s] = si.sim
  sa.val.SQf1[s] = sa.sim
  Fa.val.SQf1[s] = Fa
  car.val.SQf1[s] = car.sim
  
  mean.r.SQf1[s] <- mean(r[min(alive[,s], na.rm = TRUE):max(alive[,s], na.rm = TRUE),s])
  lambda.SQf1[s] <- exp(mean.r.SQf1[s])
  
} # s



######################################################################################################  
# calculation of quasi-extinction (Stefanos code)

sum(apply(N[3,,1:nsim],2,min)<3)/nsim # adults only

sum(apply(colSums(N[-1,,1:nsim]),2,min)<3)/nsim # immatures and adults only

sum(N[3,nyears,]==0) / nsim

sum(N[3,nyears,]<=N[3,1,1]) / nsim


mu.lambda.SQf1<-mean(lambda.SQf1)
mean(lambda.SQf1)
sd(lambda.SQf1)
# 
# #store lambda
write.csv(lambda.SQf1, "lambdaSQF1_car.csv")

summary(N[3,nyears,])
SQ.f1.summary<-summary(N[3,nyears,])

# # store persistence 
write.csv(persist.SQf1, "persistSQF1_car.csv")

######################################################################################################

# all sims are 1000, p(ext) adults only - 
# Original status quo model = 9.8
# Original field 1 model = 29.4

# 1. Status quo model with clutch rate allowed to vary *e.g. as per field 1 = 20.8

# 2. replace the clutch size >> status quo code = 17.8

# 3. keep clutch size change + changing the adult survival from elicited >> status quo vals
#  = 7.3

# 4. change up juvenile survival for status quo = 13.1

# 5. elicited juv surv estimates, change clutch size back to variable +status quo adult surv
#   = 14.4


######################################################################################################
# look at parameter distributions produced by the rbeta code (status quo)

K.sim <- rbeta(10000, beta.params(0.72, 0.12)$alpha, beta.params(0.72, 0.12)$beta)
Krnorm <-rnorm(10000, K, K.e)
plot(density(K.sim))+lines(density(Krnorm), col="blue")

summary(Krnorm)
quantile(Krnorm, probs = c(0.025, 0.975))
quantile(K.sim, probs = c(0.025, 0.975))

y.sim <- rbeta(10000, 10,5)
plot(density(K.sim))
plot(density(y.sim))
sa.sim <- rbeta(10000, beta.params(Sad, Sad.e)$alpha, beta.params(Sad, Sad.e)$beta)

# Using the stored data from the models, it is clear that SQ adult survival is better than th
# elicited data. Since adult survival is the dominant eigenvalue this will affect the models
# disproportionately

plot(density(sa.val), col="black", main = 'Discrepancy between estimated adult survival')+ lines(density(sa.val.f1), col="blue")+lines(density(sa.val.f2), col = "green") + 
  legend("topright", legend=c("Field 2", "Field 1", "SQ"),col=c("green", "blue", "black"), lty=1, cex=0.8) 
plot(density(si.val.f2), col="green")+lines(density(si.val.f1), col="blue")+lines(density(si.val))
plot(density(sj.val.f2), col="green")+lines(density(sj.val.f1), col="blue")+lines(density(sj.val))
plot(density(Fa.val.f2), col="green", xlim=c(0,0.5))+lines(density(Fa.val.f1), col="blue")+lines(density(Fa.val))+
  lines(density(Fa.val.a3), col="purple", lwd=2) + lines(density(Fac.val.a3), col="yellow", lwd=2)

# how can captive rearing fecundity values be negative???
plot(density(Fac.val.a3))


######################################

# What happens when we use Field 1 elicited values for all survival rates, with the status quo code?

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


plot_results.sqf1 = data.frame(alt = "SQf1.car", years = 1:51, N.mean=N.mean[3,], N.uci = N.uci[3,], N.lci=N.lci[3,], N.upr = N.upr[3,], N.lwr = N.lwr[3,],
                          N.25 = N.25[3,],N.50 = N.50[3,], N.75 = N.75[3,])

gg.sq.with.f1.surv<-ggplot(plot_results.sqf1, aes(x=years, y=N.mean)) + 
  geom_ribbon(aes(x = years, ymax=N.upr, ymin=N.lwr),fill = "dimgrey", alpha=0.35) + 
  geom_line(color="mediumblue", size=1.5, alpha=0.6) + theme_classic(base_size = 16)+  
  coord_cartesian(ylim=c(0,max(N.upr+10)))  + 
  scale_x_continuous(expand=c(0,0)) +  
  scale_y_continuous(expand=c(0,0)) +
  labs(title = "Predicted tara iti population size, SQ productivity + f1 survival", x = "Year", y = "Population size (no. of adult females)")

gg.sq.with.f1.surv

lambda.density.editedSQ<-plot(density(lambda.SQf1), lwd = 3, xlim=c(0.8, 1.2), main = "",  xlab = expression(lambda)) + abline(v=1, col = "grey")
plot(density(N[3,50,1:nsim]))
# 
write.csv(plot_results.sqf1, "projsqf1_car.csv")

N50_sqf1<-N[3,51,1:nsim] # store the adult population sizes at the end of simulation
N25_sqf1<-N[3,25,1:nsim] # store the adult population sizes halfway through simulation
# 
write.csv(N50_sqf1, "N50SQF1_car.csv")

#N50_sqf1<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\Fairy Terns and SDM\\Decision Making\\TARA ITI WORKSHOP\\Phase 4\\Population modelling results\\Final pop sizes 5000 sims elicited K\\Final_pop_sizes_t=50_SQF1.csv", header=T)
#N50_sqf1[,3]<-"SQ + F1 Surv"
