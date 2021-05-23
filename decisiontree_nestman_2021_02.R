
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
library(wesanderson)
library(ggeffects)
library(emmeans)
library(jagsUI)
library(MultinomialCI)
library(mc2d)
library(binom)
library(RColorBrewer)
library(Cairo)
install.packages('stringr')
library(stringr)

#################################################################################################################

##### DECISION TREE ANALYSIS FOR EGG MANAGEMENT #####

#invisible(lapply(c("jagsUI","jagsUI","RColorBrewer","Cairo"), library, character.only = TRUE))

# Create list with data for all nodes

all.nodes <- list(n.A=c(114,12,23), # original tree
     n.B=c(83,17,14), # original tree
     n.C=c(27,5,11),
     n.D=c(14,10,3),
     n.E=c(64,35), # binomial (successes, fails - check this)
     n.F=c(18,44), # natal return or donate
     n.G=c(36,2,6),
     n.H=c(18,5,13),
     n.I=c(12,2,4),
     n.J=c(5,1,6),
     n.K=c(8,0,6),
     n.L=c(2,0,6)
     )


# JAGS model
cat("
model{
  # likelihood
  
    for(i in 1:n.eggs){
      y[i]~ dcat(p[]) # categorical distribution
  }
  
  # prior
  p[1:n.categories]~ ddirch(alpha[]) # multivariate distrib Dirichlet, 
  # alpha is a vector of prior weights (here, evenly split over the three outcomes, eg no prior)
}

", file="mn.bug")

# Run and extract results
fit.mn <- post.p <- list()
for(i in 1:length(all.nodes)){
  # Data
  datalst <- list(y=rep(1:length(all.nodes[[i]]),all.nodes[[i]]),n.categories=length(unique(all.nodes[[i]])),
                  n.eggs=sum(all.nodes[[i]]),alpha=rep(1/3,length(unique(all.nodes[[i]]))))
  # alpha is a vector of prior weights (here, evenly split over the three outcomes, eg no prior)
  # Inits
  inits<-function(){list(p=rep(1/max(datalst$y),max(datalst$y)))}
  # Monitor
  params=c("p")
  # Settings
  ni=10000;nb=5000;nt=1;nc=3
  # Run and store results as a list element
  fit.mn[[i]] <- jagsUI(data=datalst, inits=inits, parameters.to.save=params,model.file="mn.bug",
                  n.chains=nc, n.iter=ni,n.burnin=nb,n.thin=nt,parallel=TRUE)
  # Extract posteriors of probabilities for this node
  post.p[[i]] <- fit.mn[[i]]$sims.list$p
  # Track progress of loop
  print(paste("Node",names(all.nodes)[i],"complete"))
  
}
# Name elements of list
names(post.p) <- paste0("p",substr(names(all.nodes),2,3))
# Attach list (so you don't have to specify "post.p") - don't worry about
attach(post.p)

# Decision tree (same n of iterations as jags model)
tree.out <- array(NA,dim=c(ni,5),dimnames=list(as.character(1:ni),c("A","B","C.a","C.b","D")))
solution <- rep(NA,ni)
for(i in 1:ni){
  tree.out[i,1] <- p.A[i,1] * p.B[i,1]              # No management
  tree.out[i,2] <- p.C[i,1] * p.D[i,1]              # Donate
  tree.out[i,3] <- p.E[i,1] * p.G[i,1] * p.H[i,1]   # AI plus Donate
  tree.out[i,4] <- p.E[i,1] * p.I[i,1] * p.J[i,1]   # AI plus Natal return
  tree.out[i,5] <- p.K[i,1] * p.L[i,1]              # Shift
  solution[i] <- colnames(tree.out)[which.max(tree.out[i,])]
}
table(solution)

# Cumulative distribution functions for all curves
c.dist <- array(NA,dim=c(ni,5))
for(i in 1:ncol(tree.out)){
  c.dist[,i] <- ecdf(unname(tree.out[,i]))(v=seq(0,1,length.out=ni))
}


## Density plot

fig1.df <- data.frame(Action=rep(c("None", "Donate","A.Inc plus Donate", "A.Inc plus Natal Return", "Shift"),each=ni),
                      Outcome=as.numeric(tree.out),
                      CDF.y=as.numeric(c.dist),
                      CDF.x=rep(seq(0,1,length.out=ni),ncol(tree.out))
                      )
# Colours for plot
my.cols <- brewer.pal(dim(tree.out)[2], "Set1") 

# Plot 1: distribution of outcomes
p1 <- ggplot(fig1.df, aes(x=Outcome, group=Action, fill=Action))+
  geom_density(alpha=0.5,colour="transparent")+
  scale_colour_manual(name="",values=my.cols)+
  scale_fill_manual(name="",values=my.cols)+
  theme_bw()+
  theme(legend.position=c(0.9,0.9),
        legend.text=element_text(size=9),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(size=16, color="black"),
        axis.title.x=element_text(vjust=-1.2, size=16, color="black"),
        axis.ticks.length=unit(0.5,"mm"),
        axis.text = element_text(size=14),
        strip.text.x = element_text(size=14),
        plot.margin=unit(c(5,5,5,7), "mm")
        #panel.margin=unit(c(4,4,4,4),"mm")
  )+
  scale_x_continuous(name="Probability of chick hatching and fledging",breaks=seq(0,0.8,0.2),limits=c(0,0.8),expand=c(0,0))+
  scale_y_continuous(name="Density",breaks=seq(0,18,by=3),limits=c(0,14),expand=c(0,0))
p1 # Check
    # Export to pdf

# cairo_pdf(file="Decision_tree_fig1_update03.pdf",width=110/25.4, height=110/25.4, pointsize=5,antialias = "subpixel")
# p1
# dev.off()

CairoPNG(file = "decision_prob_dist3.png",units="in", width=10, height=8, pointsize=6, dpi = 300, antialias = "subpixel")
p1
dev.off()

# Plot 2: cdf of outcomes
p2 <- ggplot(fig1.df, aes(x=CDF.x, y=CDF.y, group=Action, colour=Action))+
  geom_line()+
  scale_colour_manual(name="",values=my.cols)+
  theme_bw()+
  theme(legend.position=c(0.9,0.9),
        legend.text=element_text(size=9),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(size=16, color="black"),
        axis.title.x=element_text(vjust=-1.2, size=16, color="black"),
        axis.ticks.length=unit(0.5,"mm"),
        axis.text = element_text(size=14),
        strip.text.x = element_text(size=14),
        plot.margin=unit(c(5,5,5,7), "mm")
  )+
  scale_x_continuous(name="Probability of chick hatching and fledging",breaks=seq(0,0.8,0.2),limits=c(0,1),expand=c(0,0))+
  scale_y_continuous(name="Density",breaks=seq(0,1,by=0.2),limits=c(0,1),expand=c(0,0))
p2 # Check
# Export to pdf
# cairo_pdf(file="Decision_tree_fig2_update03.pdf",width=110/25.4, height=110/25.4, pointsize=5,antialias = "subpixel")
# p2
# dev.off()

CairoPNG(file = "decision_CDF2.png",units="in", width=10, height=8, pointsize=6, dpi = 300, antialias = "subpixel")
p2
dev.off()


detach(post.p)

######

## reporting mean outcomes with error ##

# no management

mu.NM<-mean(tree.out[,1])
sd.NM<-sd(tree.out[,1])
NM.uci <- mu.NM + 1.96*sd.NM                       # upper 95% confidence level (1.96 * SD)
NM.lci <- mu.NM - 1.96*sd.NM                       # lower 95% confidence level

NM.lwr <- quantile(tree.out[,1], 0.025)
NM.upr <- quantile(tree.out[,1], 0.975)

# donate

mu.don<-mean(tree.out[,2])
sd.don<-sd(tree.out[,2])
don.uci <- mu.don + 1.96*sd.don                       
don.lci <- mu.don - 1.96*sd.don      

don.lwr <- quantile(tree.out[,2], 0.025)
don.upr <- quantile(tree.out[,2], 0.975)


# ai plus donate

mu.aid<-mean(tree.out[,3])
sd.aid<-sd(tree.out[,3])
aid.uci <- mu.aid + 1.96*sd.aid                      
aid.lci <- mu.aid - 1.96*sd.aid    

aid.lwr <- quantile(tree.out[,3], 0.025)
aid.upr <- quantile(tree.out[,3], 0.975)

# ai plus natal return

mu.ainr<-mean(tree.out[,4])
sd.ainr<-sd(tree.out[,4])
ainr.uci <- mu.ainr+ 1.96*sd.ainr                      
ainr.lci <- mu.ainr - 1.96*sd.ainr

ainr.lwr <- quantile(tree.out[,4], 0.025)
ainr.upr <- quantile(tree.out[,4], 0.975)


# shift

mu.shift<-mean(tree.out[,5])
sd.shift<-sd(tree.out[,5])
shift.uci <- mu.shift + 1.96*sd.shift                     
shift.lci <- mu.shift - 1.96*sd.shift

shift.lwr <- quantile(tree.out[,5], 0.025)
shift.upr <- quantile(tree.out[,5], 0.975)


# zoo death

mu.zd<-mean(p.E[,1])
sd.zd<-sd(p.E[,1])
zd.uci<-mean(p.E[,1]) + 1.96*sd(p.E[,1])
zd.lci<-mean(p.E[,1]) - 1.96*sd(p.E[,1])

zd.lwr <- quantile(p.E[,1], 0.025)
zd.upr <- quantile(p.E[,1], 0.975)


#################################################################################################################

## Counterfactuals ##

## counterfactual 1): manager perception of threat is 100% accurate, and if eggs that were managed had not been, 
## there would be a 100% failure rate

## total eggs are 83 + 0 / 305
## 83/305 # 0.27

## counterfactual 2): managers perceived a threat but it was not really there, in which case unmanaged eggs would actually 
## experience a background rate of survival to fledging, same as 'no nest management' branch

## chicks fledged  n = 156 * mu.NM (lci - uci, 156*NM.lci, 156*NM.uci), + 83 fledged from 'none' branch
## 86.3 (73.8-98.7) + 83/305 = 169.3 (156.8-181.7)/305 = 0.56 (0.51-0.6)

##########################################################

## Scenarios of improvement ##

## scenario 1): managers improve threat detection and no eggs are lost to predators or tides because they are managed before
## they die

## n=12 eggs were (known to be) lost to tide, storm or predator (the actual number will be higher) - these eggs are 
## experience probability of survival as a mean of the managed branches (not shift)

((mu.don+mu.aid+mu.ainr)/3) * 12 # 2.99 extra chicks fledge
((don.lwr+aid.lwr+ainr.lwr)/3) * 12 # 1.69 lower
((don.upr+aid.upr+ainr.upr)/3) * 12 # 4.52 upper

## scenario 2): managers improve 'artificial incubation and return/donate' success rate
## in this scenario, eggs that were managed through Auckland Zoo (n=99) have an improved success rate 

# A) in line with 'donate' management (50% improvement)
((mu.aid+mu.ainr)/2)*1.5 # 0.32 - same as donate mean

((mu.don*99)+14+2+83)/305
((don.lwr*99)+14+2+83)/305
((don.upr*99)+14+2+83)/305

((mu.don*99)+14+2+83)
((don.lwr*99)+14+2+83)
((don.upr*99)+14+2+83)

# B) 75% improvement in efficacy

((mu.aid+mu.ainr)/2)*1.75 # 0.38 

((0.38*99)+14+2+83)/305 # = 0.45 proportional improvement in fledging success

## scenario 3): managers improve habitat management, reducing threat presence

## Calculating the threshold when you should start managing
# 29 points between -13 and +16. We would start managing when the increase
# in fledglings goes to 0
16/29

# 55% certain before managing

#####################################################################################################

## decision tree scenarios graph

theme_set(theme_bw(base_size=20, base_family='Times New Roman')+ 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# load the dataframe

scenario_df<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\02_chapter-two\\02-decision-tree\\dataframes\\dectree-scenarios.csv", header = T)


# make a bar plot with error bars

# #reorder
scenario_df$scenario = with(scenario_df, reorder(scenario, fledglings))



# another way to reorder
scenario_df$scenario<- factor(scenario_df$scenario, levels = c("Improve AI efficacy", 
                                                               "Improve true positive threat detection",
                                                               "Counterfactual 2","Counterfactual 1", "Observed outcomes"))
scenario_df$scenarionew<-str_wrap(scenario_df$scenario, width = 16)

f3<-ggplot(scenario_df, aes(scenario, fledglings))+
  geom_bar(fill = "orangered", stat = "identity", width = 0.35, alpha = 0.5) +
  geom_linerange(aes(scenario,ymin =f.lci, ymax = f.uci))+  
  scale_y_continuous(name="Number of fledglings",limits = c(0,195),expand=c(0,0))+
  scale_x_discrete(name = "Scenario")+
  theme_bw()+
  theme(legend.position=c(0.8,0.6),
        legend.text=element_text(size=9),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(size=20, color="black"),
        axis.title.x=element_text(vjust=-1.2, size=20, color="black"),
        axis.ticks.length=unit(0.5,"mm"),
        axis.text = element_text(size=16),
        strip.text.x = element_text(size=18),
        plot.margin=unit(c(5,5,5,7), "mm"))+
          coord_flip() 

cairo_pdf(file="scenario_test.pdf",width=110/25.4, height=110/25.4, pointsize=5,antialias = "subpixel")
f3
dev.off()

CairoPNG(file = "scenario_test6.png",units="in", width=12, height=10, pointsize=6, dpi = 300, antialias = "subpixel")
f3
dev.off()

## code with text wrap 
# f3<-ggplot(scenario_df, aes(stringr::str_wrap(scenario,16), fledglings))+
#   geom_bar(fill = "orangered", stat = "identity", width = 0.35, alpha = 0.5) +
#   geom_linerange(aes(x=stringr::str_wrap(scenario,16),ymin =f.lci, ymax = f.uci))+  
#   scale_y_continuous(name="Number of fledglings",limits = c(0,195),expand=c(0,0))+

###################################################################################################################

## PART ONE - GET THE DATA READY FOR EACH UNCERTAINTY NODE ##

### skip to end of part one to upload the dataframe from hard drive

### 1. load the ready-made df's

NZFT_df<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\02_chapter-two\\02-decision-tree\\dataframes\\NZFT.R.csv",header=T, strip.white=TRUE)
NZFT_fertile_df<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\02_chapter-two\\02-decision-tree\\dataframes\\NZFT.fertile.R.csv",header=T, strip.white=TRUE)
NZFT_df_INC_MULTI<-read.csv("C:\\Users\\thmcm\\Documents\\PhD - seabird conservation management\\02_chapter-two\\02-decision-tree\\dataframes\\dec_tree_aimulti.csv",header=T, strip.white=TRUE)

### remove the one nest that has "add" as management from all datasets ###

NZFT_fertile_df<-NZFT_fertile_df[!NZFT_fertile_df$egg.management=="add",]
NZFT_fertile_df<-NZFT_fertile_df[!NZFT_fertile_df$nest_ID=="288",]
NZFT_fertile_df<-NZFT_fertile_df[!NZFT_fertile_df$nest_ID=="287",]

NZFT_df_INC_MULTI<-NZFT_df_INC_MULTI[!NZFT_df_INC_MULTI$nest_ID=="288",]
NZFT_df_INC_MULTI<-NZFT_df_INC_MULTI[!NZFT_df_INC_MULTI$nest_ID=="287",]

## nest 247 was marked as NA for no_eggs_fail - corrected to '0'
## R17 corrected to '1'
## AH35 corrected to '0'
## AB57 corrected to '2'
## R26 corrected to '0'

NZFT_df_INC_MULTI %>%
  filter(!is.na(no_eggs_fail)) %>% # there are NAs in the fertile database
  summarise(n = sum(fertile), nhatch = sum(no_eggs_hatched), nfail = sum(no_eggs_fail))

#################################################################################################################

## PART TWO - RETRIEVE THE VALUES FOR EACH UNCERTAINTY NODE ##

# All eggs are candled and presumed to be fertile

## prob A. Management 'NONE' egg survival fate ##

unman_egg_outcomes<-NZFT_fertile_df %>%
  filter(eggman !="1") %>%
  filter(!is.na(no_eggs_fail)) %>% # there are NAs in the fertile database
  summarise(n = sum(total_no_eggs), nhatch = sum(no_eggs_hatched), nfail = sum(no_eggs_fail))

## total number of unmanaged eggs = 149
## total number of unmanaged eggs hatch = 113
## total fail = 33

unman_egg_outcomes

unman_egg_threats_byegg<-
  NZFT_fertile_df %>%
  filter(!is.na(eggfailcause)) %>% # remove NAs
  filter(eggman !="1") %>%
  group_by(eggfailcause) %>%
  summarise(n = sum(no_eggs_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

unman_egg_threats_byegg
# First node proportions

# hatch = 0.774
113/149

# fail (predator, tide or storm) = 0.082
# pred = 2+2=4, tide/storm =8, total = 12
12/149

# fail, other or unknown  = 0.144
21/149

0.774+0.082+0.144

#####################
## prob B. Management 'NONE' chick survival ##

## hatched eggs  = 113

unman_chick_outcomes<-NZFT_fertile_df %>%
  filter(eggman !="1") %>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(no_eggs_hatched), nfledge = sum(no_chicks_fledge), nfail = sum(no_chicks_fail))
unman_chick_outcomes

# Second node proportions

# fledge = 0.72
81/112

unman_egg_threats_bychick<-
  NZFT_fertile_df %>%
  filter(!is.na(chickfailcause)) %>% # remove NAs
  filter(eggman !="1") %>%
  group_by(chickfailcause) %>%
  summarise(n = sum(no_chicks_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

unman_egg_threats_bychick

# fail (predator, tide, storm) = 0.152
(8 + 9)/112

# fail other, unknown = (n = 14) 0.128
11+1+1+1
1-(0.152+0.72)

####################################################################

######## DONATE branch
## prob C. Management 'DONATE' egg survival fate ##

# overall outcomes, eggs

donate_egg_outcomes<-NZFT_fertile_df %>%
  filter(egg.management =="donate") %>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(total_no_eggs), nhatch = sum(no_eggs_hatched), nfail = sum(no_eggs_fail))

## total number of donated eggs = 43
## total number of donated eggs hatch = 27
## total fail = 16
donate_egg_outcomes

donate_egg_threats_byegg<-
  NZFT_fertile_df %>%
  filter(!is.na(eggfailcause)) %>% # remove NAs
  filter(egg.management =="donate") %>%
  group_by(eggfailcause) %>%
  summarise(n = sum(no_eggs_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

donate_egg_threats_byegg

# fail, tide, storm, pred = 5

# fail, unknown, other = 11

27/43 # 0.628 hatch
5/43 # 0.116
11/43 # 0.256

########################
## prob D. Management 'DONATE' chick survival fate ##

# overall outcomes, chicks

donate_chick_outcomes<-NZFT_fertile_df %>%
  filter(egg.management =="donate") %>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(no_eggs_hatched), nfledge = sum(no_chicks_fledge), nfail = sum(no_chicks_fail))

# number of chicks 27
# number fledge 14
# number fail 13
donate_chick_outcomes

donate_chick_threats_byegg<-
  NZFT_fertile_df %>%
  filter(!is.na(chickfailcause)) %>% # remove NAs
  filter(egg.management =="donate") %>%
  group_by(chickfailcause) %>%
  summarise(n = sum(no_chicks_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

donate_chick_threats_byegg

# fledged

14/27 # 0.519

# fail cause tide/storm, pred = 6+3+1 = 10

10/27 # 0.370 

# fail cause unknown / other = 3

3/27 # 0.111

## multinomial confidence interval D

# prob.d<-c(14,10,3)
# multinomialCI(prob.d, 0.05, verbose = FALSE)


## Management AI and MULTI branch #######################
# 61 lines of data

inc_multi_egg_outcomes<-NZFT_fertile_df %>%
  filter(egg.management %in% c("ai", "multi")) %>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(total_no_eggs), nhatch = sum(no_eggs_hatched), nfail = sum(no_eggs_fail))
inc_multi_egg_outcomes

NZFT_df_INC_MULTI%>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(total_no_eggs), nhatch = sum(no_eggs_hatched), nfail = sum(no_eggs_fail))


## prob E. Management 'AI' egg survival fate  - zoo death

## total number of eggs that were managed either  A.INC or MULTI = 99


#############################
## ZOO DEATH uncertainty node E


## For eggs marked zoo death, this will have to include some clutches that died hatching, or lived briefly in the zoo as a chick as well,
## even if they are down in the database as chick? ... think... they will get double counted as dying at egg and dying at chick

sum(NZFT_df_INC_MULTI$no_eggs_zoodeath) # 35

## therefore 64 eggs make it to natal return or donate

n.incmulti<-99
n.zoodeath<-35

#############################

## prob F. Management 'AI' egg decision node - NATAL RETURN or DONATE

sum(NZFT_df_INC_MULTI$no_eggs_natalreturn) # 18
sum(NZFT_df_INC_MULTI$no_eggs_fostered) # 43 
sum(NZFT_df_INC_MULTI$no_eggs_donated) # 48 'up for donation'
sum(NZFT_df_INC_MULTI$total_no_eggs) # 99 eggs

#############################
## prob. G. Management AI plus DONATE egg survival outcomes


NZFT_df_INC_MULTI %>%
  filter(egg.management=="multi") %>%
  filter(no_eggs_donated >0) %>%
  summarise(n = sum(no_eggs_fostered), nhatch = sum(no_hatchlings), nfledge = sum(no_chicks_fledge))

NZFT_df_INC_MULTI %>%
  filter(egg.management=="multi") %>%
  filter(no_eggs_donated >0) %>%
  summarise(n = sum(no_eggs_donated), nhatch = sum(no_hatchlings), nfledge = sum(no_chicks_fledge))

### n = 44 fostered
### n = 48 donated (4 eggs must be zoo death)

## IF binomial node:
## 36 hatched
36/43

# 7 failed

## IF Multinomial node:
# p.J<-(36,2,5)

NZFT_df_INC_MULTI %>%
  filter(egg.management=="multi") %>%
  filter(no_eggs_fostered >0) %>% # or 'donated'
  filter(!is.na(eggfailcause)) %>% # remove NAs
  filter(no_eggs_zoodeath=="0") %>%
  group_by(eggfailcause) %>%
  summarise(n = sum(no_eggs_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

# tide/storm = 2
# other/unknown = 5

##########################
## prob. H. Management AI plus DONATE chick survival outcomes

NZFT_df_INC_MULTI %>%
  filter(egg.management=="multi") %>%
  filter(no_eggs_donated >0) %>%
  filter(!is.na(chickfailcause)) %>% # remove NAs
  group_by(chickfailcause) %>%
  summarise(n = sum(no_chicks_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

## IF binomial node:

## 18 fledged
18/36

## IF Multinomial node:

# 18 fledged
# 5 died tide/storm/pred
# 13 died unknown other

# p.K<-(18,8,13)

############################
## prob. I Management AI plus NATAL RETURN egg survival outcomes

NZFT_df_INC_MULTI %>%
  filter(no_eggs_natalreturn >0) %>%
  summarise(n = sum(no_eggs_natalreturn), nhatch = sum(no_hatchlings), nfledge = sum(no_chicks_fledge))

## 12 eggs hatched after natal return = 0.667
12/18

NZFT_df_INC_MULTI %>%
  filter(no_eggs_natalreturn >0) %>%
  filter(!is.na(eggfailcause)) %>% # remove NAs
  group_by(eggfailcause) %>%
  summarise(n = sum(no_eggs_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

## of eggs returned to the wild after AI 
# 2 embryo death, 1 abandon, 1 unknown
# 2 pred unknown

## prob. J Management A.INC plus NATAL RETURN chick survival outcomes

## 5 chicks fledged after natal return = 0.417

5/12

NZFT_df_INC_MULTI %>%
  filter(no_eggs_natalreturn >0) %>%
  filter(!is.na(chickfailcause)) %>% # remove NAs
  group_by(chickfailcause) %>%
  summarise(n = sum(no_chicks_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)


###############################################################################

## Management 'SHIFT' branch in decision tree

# prob. K Management SHIFT chick outcomes

NZFT_fertile_df %>% filter(egg.management=="shift") %>% summarise(n = sum(total_no_eggs))
NZFT_df_INC_MULTI %>% summarise(n = sum(no_eggs_shift))

# 14 eggs were shifted with no further action taken

shift_outcomes<-NZFT_fertile_df %>%
  filter(egg.management=="shift") %>%
  #filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(no_eggs_hatched), nfledge = sum(no_chicks_fledge), nfail = sum(no_chicks_fail))

shift_outcomes

# Of 14 eggs shifted with no further management:
# 8 hatched
# 2 fledged

## fail outcomes

NZFT_fertile_df %>% 
  filter(egg.management=="shift") %>% 
  filter(!is.na(eggfailcause)) %>% # remove NAs
  group_by(eggfailcause) %>%
  summarise(n = sum(no_eggs_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

##############

# prob. L Management SHIFT chick outcomes
# all 6 are unknown/other (embryo death)

NZFT_fertile_df %>% 
  filter(egg.management=="shift") %>% 
  filter(!is.na(chickfailcause)) %>% # remove NAs
  group_by(chickfailcause) %>%
  summarise(n = sum(no_chicks_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

## again, all 6 are unknown / other (health)

# 22 were shifted then brought into captivity - these have gone into MULTI

#################################################################################################################

## managed egg summaries

man_egg_outcomes<-NZFT_fertile_df %>%
  filter(eggman =="1") %>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(total_no_eggs), nhatch = sum(no_eggs_hatched), nfail = sum(no_eggs_fail))
man_egg_outcomes

man_chick_outcomes<-NZFT_fertile_df %>%
  filter(eggman =="1") %>%
  filter(!is.na(no_eggs_fail)) %>%
  summarise(n = sum(no_eggs_hatched), nfledge = sum(no_chicks_fledge), nfail = sum(no_chicks_fail))
man_chick_outcomes


man_egg_threats_byegg<-
  NZFT_fertile_df %>%
  filter(!is.na(eggfailcause)) %>% # remove NAs
  filter(eggman =="1") %>%
  group_by(eggfailcause) %>%
  summarise(n = sum(no_eggs_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

man_egg_threats_byegg

## if 159 managed eggs all fail, they would be divided as 26% due to tide/pred and 74% other
## which equals 41 (tide pred) and 118 (other)

man_egg_threats_bychick<-
  NZFT_fertile_df %>%
  filter(!is.na(chickfailcause)) %>% # remove NAs
  filter(eggman =="1") %>%
  group_by(chickfailcause) %>%
  summarise(n = sum(no_chicks_fail)) %>%
  mutate(freq = n / sum(n)) %>% # create a proportion
  arrange(desc(freq)) # arrange by most to least common)

man_egg_threats_bychick


####################################

sum(NZFT_df$total_no_eggs)
sum(NZFT_fertile_df$total_no_eggs, na.rm=T)
