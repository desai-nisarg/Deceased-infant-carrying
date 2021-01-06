####----INTRODUCTION----#####
# Manuscript title: "Deceased-Infant Carrying in Nonhuman Anthropoids: Insights From Systematic Analysis and Case Studies of Bonnet Macaques (Macaca radiata) and Lion-Tailed Macaques (Macaca silenus)"
# Manuscript authors: Sayantan Das, Joseph J. Erinjery, Nisarg Desai, Kamaraj Mohan, Honnavalli N. Kumara, Mewa Singh
# Code author: Nisarg Desai, desai054[at]umn[dot]edu
# Last update: 16 March 2018

rm(list=ls())
library(ape)
library(MCMCglmm)
library(phytools)
library(VGAM)
setwd("...to the directory containing nexus tree file and data...")

# Importing raw data
data<-read.csv("DIC_data.csv",header=TRUE) # Read raw data
data$animal<-gsub(" ", "_", data$species) # Make species names readable by R
data$species<-gsub(" ", "_", data$species) # Make another identical column for additional random effect in MCMCglmm
data$lDIC<-log(data$DIC) # log transform the dependent variable
View(data)

data1<-data[complete.cases(data),] # Remove rows with missing values
View(data1) # Pan_paniscus is removed

# Importing the tree
phylo<-read.nexus("DIC_consensusTree_1000_Version3.nex")
plot.phylo(phylo, edge.color="red", tip.color="blue", edge.width=2, font=3) # Plot with polytomy
phylo<- drop.tip(phylo, "Pan_paniscus") # Get rid of Pan_paniscus in the tree
plot.phylo(phylo, edge.color="red", tip.color="blue", edge.width=2, font=3)
is.ultrametric(phylo) # Need ultrametric tree

Resolved2<- multi2di(phylo, random=FALSE) # Break polytomy
plot(Resolved2)
finaltree <- compute.brlen(Resolved2, method="Grafen") # Make tree ultrametric
plot.phylo(finaltree, edge.color="red", tip.color="blue", edge.width=2, font=3)
is.ultrametric(finaltree) # Yaay!



# Prep for MCMCglmm

# Number of interations
nitt<-240000

# Length of burnin
burnin<-40000

# Amount of thinning
thin<-100

# Priors
prior<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),
             R=list(V=1,nu=0.02))

# Ready to fit the models

formula_1 <- lDIC~Age.Cat+Parity # Model 1 for maternal characteristics
formula_2 <- lDIC~Sex.of.infant+Age.at.death # Model 2 for infant characteristics
formula_3 <- lDIC~Degree.of.arboreality+Condition+Context+T.2 # Model 3 for ecological conditions

# Model 1

fit1a <- MCMCglmm(fixed = formula_1, random=~animal+species, data=data1, pedigree=finaltree,nitt = nitt,
                  burnin = burnin,thin = thin, prior = prior)
fit1b <- MCMCglmm(fixed = formula_1, random=~animal+species, data=data1, pedigree=finaltree,nitt = nitt,
                  burnin = burnin,thin = thin, prior = prior) # Re-run the model to check for convergence
summary(fit1a)


gelman.diag(mcmc.list(fit1a$Sol, fit1b$Sol)) # Check convergence for fixed effects
gelman.diag(mcmc.list(fit1a$VCV, fit1b$VCV)) # Check convergence for random effects

H2_fit1<- var(fit1a$VCV[,'animal'])/
  (var(fit1a$VCV[,'animal'])+var(fit1a$VCV[,'species'])+
     var(fit1a$VCV[,'units'])) #H square is similar to Pagel's lambda

# Model 2

fit2a <- MCMCglmm(fixed = formula_2, random=~animal+species, data=data1, pedigree=finaltree,nitt = nitt,
                  burnin = burnin,thin = thin, prior = prior)
fit2b <- MCMCglmm(fixed = formula_2, random=~animal+species, data=data1, pedigree=finaltree,nitt = nitt,
                  burnin = burnin,thin = thin, prior = prior) # Re-run the model to check for convergence
summary(fit2a)


gelman.diag(mcmc.list(fit2a$Sol, fit2b$Sol)) # Check convergence for fixed effects
gelman.diag(mcmc.list(fit2a$VCV, fit2b$VCV)) # Check convergence for random effects

H2_fit2<- var(fit2a$VCV[,'animal'])/
  (var(fit2a$VCV[,'animal'])+var(fit2a$VCV[,'species'])+
     var(fit2a$VCV[,'units'])) #H square is similar to Pagel's lambda

# Model 3

fit3a <- MCMCglmm(fixed = formula_3, random=~animal+species, data=data1, pedigree=finaltree,nitt = nitt,
                  burnin = burnin,thin = thin, prior = prior)
fit3b <- MCMCglmm(fixed = formula_3, random=~animal+species, data=data1, pedigree=finaltree,nitt = nitt,
                  burnin = burnin,thin = thin, prior = prior) # Re-run the model to check for convergence
summary(fit3a)


gelman.diag(mcmc.list(fit3a$Sol, fit3b$Sol)) # Check convergence for fixed effects
gelman.diag(mcmc.list(fit3a$VCV, fit3b$VCV)) # Check convergence for random effects

H2_fit3<- var(fit3a$VCV[,'animal'])/
  (var(fit3a$VCV[,'animal'])+var(fit3a$VCV[,'species'])+
     var(fit3a$VCV[,'units'])) #H square is similar to Pagel's lambda

