#this is my data analysis for my dissertation
#Sarita Whitehead - 01/10/2022


install.packages("readxl")

library(nimble)
library(tidyverse)
library(MCMCvis)
library(readxl)
library(dplyr)

#loading my data
Fulmar <- read_xlsx("Fulmar.xlsx")
View(Fulmar)

#set NA to 0
is.na(Fulmar)
#which values in my object are true
is.na (Fulmar)
Fulmar[is.na(Fulmar)] <- 0
View(Fulmar)

#change all values of 2 in data to 1
View(Fulmar)
#change values of 2 to 1 in one column 
Fulmar["1973"][Fulmar["1973"] == 2] <- 1
View(Fulmar)
#change values of 2 to 1 in columns 4 to 50 (1967-2018)
Fulmar [, 4:55][Fulmar[, 4:55] == 2] <- 1
View(Fulmar)
Fulmar [, 4:55]
Fulmar [ (29),50:55]

#move first three columns to end
copyofdata <- Fulmar
Fulmar<- Fulmar[,-c(1:3)]
View(Fulmar)
View(copyofdata)
Fulmar[,52 +1]<-copyofdata[,2]
View(Fulmar)

#Rename the years
colnames(Fulmar) <- c("year_1967", "year_1968", "year_1969", 
                    "year_1970", "year_1971", "year_1972",
                    "year_1973", "year_1974", "year_1975", "year_1976", 
                    "year_1977", "year_1978", "year_1979", "year_1980",
                    "year_1981", "year_1982", "year_1983", "year_1984",
                    "year_1985", "year_1986", "year_1987", "year_1988",
                    "year_1989", "year_1990", "year_1991", "year_1992",
                    "year_1993", "year_1994", "year_1995","year_1996",
                    "year_1997","year_1998", "year_1999", "year_2000",
                    "year_2001", "year_2002", "year_2003", "year_2004",
                    "year_2005", "year_2006", "year_2007", "year_2008",
                    "year_2009", "year_2010", "year_2011", "year_2012",
                    "year_2013", "year_2014", "year_2015", "year_2016",
                    "year_2017", "year_2018", "sex")
##### phi(.)p(.)--------------
hmm.phip <- nimbleCode({
  
  #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  phi ~ dunif(0, 1) # prior survival
  #Survival matrix
  gamma[1,1] <- phi # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
 
   #Recapture
  p ~ dunif(0, 1) # prior detection
  #Recapture matrix
  omega[1,1] <- 1 - p # Pr(alive t -> non-detected t)
  omega[1,2] <- p # Pr(alive t -> detected t)
  omega[2,1] <- 1 # Pr(dead t -> non-detected t)
  omega[2,2] <- 0 # Pr(dead t -> detected t)
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

# Format data for CJS
y <- Fulmar %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)


#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first)
my.constants

#Now the data in a list. Note that we add 1 to the data
#to have 1 for non-detections and 2 for detections.
#You may use the coding you prefer of course, you will just need to
#adjust the $\Omega$ and $\Gamma$ matrices in the model above.
my.data <- list(y = y + 1)
View(my.data)

#Specify initial values. For the latent states, we go for the easy way,
#and say that all individuals are alive through the study period.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
initial.values()

#Some information that we now pass as initial value
#info (observations of alive) are actually known states,
#and could also be passed as data in which case the initial values have to be 0.
#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p")
parameters.to.save

# MCMC details
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#Let's run nimble.
mcmc.phip <- nimbleMCMC(code = hmm.phip,
                        constants = my.constants,
                        data = my.data,
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin,
                        nchains = n.chains)
#Examine the results.
MCMCsummary(mcmc.phip, round = 2)
MCMCtrace(mcmc.phip, pdf = F)


#model 1 - survival is constant and probability of recapture is dependent on sex

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(.)p(s)--------------
hmm.phips <- nimbleCode({
 
   #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  
  #Survival
  phi ~ dunif(0,1) # prior survival
  
  #Survival matrix
  gamma[1,1] <- phi # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
  #Recapture depends on sex
  for(i in 1:N){
    logit(p[i]) <- beta[sex[i]]
   
     #Observation matrix
    omega[1,1,i] <- 1 - p[i] # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i] # Pr(alive t -> detected t)
    omega[2,1,i] <- 1 # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0 # Pr(dead t -> detected t)
  }
  
  # Priors for beta (recapture changes with sex, so we need two betas;
  #beta[sex[i] -> beta[1] and beta[2]])
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  
  # inverse logit for transforming p estimate
  p_male <- ilogit(beta[1])
  p_female <- ilogit(beta[2])
 
   #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)

#Data in a list. We add 1 to the data to have 1 for non-detections and 
#2 for detections.
my.data <- list(y = y + 1)
View(my.data)

#Initial values
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  phi = runif(1,0,1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("beta", "phi", "p_male", "p_female")
parameters.to.save

#MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#At last, let's run nimble.
mcmc.phips <- nimbleMCMC(code = hmm.phips,
                         constants = my.constants,
                         data = my.data,
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin,
                         nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.phips, round = 2)
MCMCtrace(mcmc.phips, pdf=F)

#model 2 survival is constant and recapture is dependent on time

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(.)p(t)--------------
hmm.phipt <- nimbleCode({
  #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  #Survival
  phi ~ dunif(0, 1) # Prior for survival
  #Survival matrix
  gamma[1,1] <- phi # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
  
  #Recapture
  for(t in 1:(T-1)){
    p[t] ~ dunif(0,1) # Prior for p.
    #Recapture matrix
    omega[1,1,t] <- 1 - p[t] # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t] # Pr(alive t -> detected t)
    omega[2,1,t] <- 1 # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0 # Pr(dead t -> detected t)
  }
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first)
my.constants

#Now the data in a list.
my.data <- list(y = y + 1)

#Initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p")
parameters.to.save

#MCMC details.
n.iter <- 8000
n.burnin <- 1000
n.chains <- 2

#At last, let's run nimble.
mcmc.phipt <- nimbleMCMC(code = hmm.phipt,
                         constants = my.constants,
                         data = my.data,
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin,
                         nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.phipt, round = 2)
MCMCtrace(mcmc.phipt, pdf=F)

#model 3 survival is constant and recapture = s+t

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(.)p(s+t)--------------
hmm.phips_t <- nimbleCode({
  #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
#Survival
  phi ~ dunif(0,1) # prior survival

#Survival matrix
  gamma[1,1] <- phi # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
  
# Recapture
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta[sex[i]]+ lambda[t] #additive time + sex

#Recapture matrix
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
  
#Priors for beta
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)

  # Time fixed effect.
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(0, sd = 1.5)
  }
  
  # ilogit for p.
  for (t in 1:(T-1)){
    p_male[t] <- ilogit(beta[1]+ lambda[t])
    p_female[t] <- ilogit(beta[2] + lambda[t])
  }
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#Now the data in a list.
my.data <- list(y = y + 1)
View(my.data)

#Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  phi = runif(1,0,1),
                                  lambda = rnorm(my.constants$T-1, 0, 1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p_male", "p_female")
parameters.to.save

#MCMC details.
n.iter <- 15000
n.burnin <- 5000
n.chains <- 2

#At last, let's run nimble.
mcmc.phips_t <- nimbleMCMC(code = hmm.phips_t,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = parameters.to.save,
                           niter = n.iter,
                           nburnin = n.burnin,
                           nchains = n.chains)


#Examine the results.
MCMCsummary(mcmc.phips_t, round = 2)
MCMCtrace(mcmc.phips_t,pdf=F)

#Model 4: Survival is constant and recapture = s*t
#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(.)p(s*t)--------------
hmm.phipst <- nimbleCode({
  #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  phi ~ dunif(0,1) # prior survival
  
  #Survival matrix
  gamma[1,1] <- phi # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
  
  # Recapture
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta[sex[i]]+ lambda[t] + kappa[sex[i],t] #interaction sex * time
      
      #Recapture matrix
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
  
  #Priors for beta
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  
  #Time fixed effect.
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(0, sd = 1.5)
  }
  
  #Time as random effect in the interaction
  lambda.sigma ~ dunif(0, 10)
  for(i in 1:2){
    for (t in 1:(T-1)){
      kappa[i,t] ~ dnorm(0, sd = lambda.sigma)
    }
  }
  
  # ilogit for p.
  for (t in 1:(T-1)){
    p_male[t] <- ilogit(beta[1]+ lambda[t] + kappa[1,t])
    p_female[t] <- ilogit(beta[2] + lambda[t] + kappa[2,t])
  }
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants
str(my.constants)
#Now the data in a list.
my.data <- list(y = y + 1)
str(my.data)
#Specify initial values
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  phi = runif(1,0,1),
                                  lambda = rnorm(my.constants$T-1, 0, 1),
                                  lambda.sigma = runif(1,0,1),
                                  kappa = matrix(rnorm(102, 0, 1), 2, 51),
                                  z = zinits)
initial.values() 
str(initial.values)
#note kappa is a 2-row matrix: [1,] -> male; [2,] -> female
# with number columns = time occasions.

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p_male", "p_female")
parameters.to.save

#MCMC details.
n.iter <- 15000
n.burnin <- 5000
n.chains <- 2

#At last, let's run nimble.
mcmc.phipst <- nimbleMCMC(code = hmm.phipst,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

#Examine the results.
MCMCsummary(mcmc.phipst, round = 2)
MCMCtrace(mcmc.phipst,pdf=F)

#Survival (phi) is deponent on time, probability of recapture (p) is constant 

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(t)p(.)--------------
hmm.phitp <- nimbleCode({
  #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  # Survival
  for(t in 1:(T-1)){
    phi[t] ~ dunif(0,1) # prior survival
    #Survival matrix
    gamma[1,1,t] <- phi[t] # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t] # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1 # Pr(dead t -> dead t+1)
  }

  #Recapture matrix
  p ~ dunif(0, 1) # prior detection
  omega[1,1] <- 1 - p # Pr(alive t -> non-detected t)
  omega[1,2] <- p # Pr(alive t -> detected t)
  omega[2,1] <- 1 # Pr(dead t -> non-detected t)
  omega[2,2] <- 0 # Pr(dead t -> detected t)
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first)
my.constants

#Now the data in a list.
my.data <- list(y = y + 1)

#'Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p")
parameters.to.save

#MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#At last, let's run nimble.
mcmc.phitp <- nimbleMCMC(code = hmm.phitp,
                         constants = my.constants,
                         data = my.data,
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin,
                         nchains = n.chains)
#Examine the results.
MCMCsummary(mcmc.phitp, round = 2)
MCMCtrace(mcmc.phitp, params = "p", pdf=F)

#model 6: survival(phi) is dependent on sex and recapture (P) is constant

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(s)p(.)--------------
hmm.phisp <- nimbleCode({
 
   # Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  for(i in 1:N){
    logit(phi[i])<- beta[sex[i]]
   
     #Survival matrix
    gamma[1,1,i] <- phi[i] # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i] # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1 # Pr(dead t -> dead t+1)
  }
 
   # Priors for beta
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
 
   # ilogit for phi
  phi_male <- ilogit(beta[1])
  phi_female <- ilogit(beta[2])
  
  #Recapture
  p ~ dunif(0,1) # prior for p
  
  # Recapture matrix
  omega[1,1] <- 1 - p # Pr(alive t -> non-detected t)
  omega[1,2] <- p # Pr(alive t -> detected t)
  omega[2,1] <- 1 # Pr(dead t -> non-detected t)
  omega[2,2] <- 0 # Pr(dead t -> detected t)
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("beta", "phi_male", "phi_female", "p")
parameters.to.save

#' MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

mcmc.phisp <- nimbleMCMC(code = hmm.phisp,
                         constants = my.constants,
                         data = my.data,
                         inits = initial.values,
                         monitors = parameters.to.save,
                         niter = n.iter,
                         nburnin = n.burnin,
                         nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.phisp, round = 2)
MCMCtrace(mcmc.phisp, params = "all", pdf=F)

#model 7: survival (phi) is dependent on t and recapture (p) is dependent on s

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(t)p(s)--------------
hmm.phitps <- nimbleCode({
 
   #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  for(t in 1:(T-1)){
    phi[t] ~ dunif(0,1) # prior for phi
  
      #Survival matrix
    gamma[1,1,t] <- phi[t] # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t] # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1 # Pr(dead t -> dead t+1)
  }
  
  #Recapture
  for(i in 1:N){
    logit(p[i]) <- beta[sex[i]]
   
     #Observation matrix
    omega[1,1,i] <- 1 - p[i] # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i] # Pr(alive t -> detected t)
    omega[2,1,i] <- 1 # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0 # Pr(dead t -> detected t)
  }
 
   #Priors for beta
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
 
   # ilogit for p
  p_male <- ilogit(beta[1])
  p_female <- ilogit(beta[2])
  
  #likelihood 
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#Now the data in a list.
my.data <- list(y = y + 1)

#Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  beta = rnorm(2,0,1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p_male", "p_female", "beta")
parameters.to.save

#MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#At last, let's run nimble.
mcmc.phitps <- nimbleMCMC(code = hmm.phitps,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

#Examine the results.
MCMCsummary(mcmc.phitps, round = 2)
MCMCtrace(mcmc.phitps, params = "all",pdf=F)

#model 8: survival (phi) is dependent on t and recapture (p) is dependent on t

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(t)p(t)--------------
hmm.phitpt <- nimbleCode({
  
  #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  for(t in 1:(T-1)){
    phi[t] ~ dunif(0,1) # prior for phi
  
      #Survival matrix
    gamma[1,1,t] <- phi[t] # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t] # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1 # Pr(dead t -> dead t+1)
  }
  
  #Recapture
  for(t in 1:(T-1)){
    p[t] ~ dunif(0,1) # prior for p

        # Recapture matrix
    omega[1,1,t] <- 1 - p[t] # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t] # Pr(alive t -> detected t)
    omega[2,1,t] <- 1 # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0 # Pr(dead t -> detected t)
  }
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

#get the occasion of first capture for all individuals
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first)
my.constants

#Now the data in a list. N
my.data <- list(y = y + 1)

#Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(my.constants$T-1,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p")
parameters.to.save

#MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#At last, let's run nimble.
mcmc.phitpt <- nimbleMCMC(code = hmm.phitpt,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

#Examine the results.
MCMCsummary(mcmc.phitpt, round = 2)
MCMCtrace(mcmc.phipt, params = "all",pdf=F)

#Model 9: Phi is dependent on time, P is sex +time

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(t)p(s+t)--------------
hmm.phitps_t <- nimbleCode({
  
  #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0

    #Survival
  for(t in 1:(T-1)){
    phi[t] ~ dunif(0,1) # prior for phi
  
      #Survival matrix
    gamma[1,1,t] <- phi[t] # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t] # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1 # Pr(dead t -> dead t+1)
  }
 
   #Recapture
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta[sex[i]] + lambda[t]
    
        #Recapture matrix
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
 
   #Priors for beta
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  
  # Time fixed effect.
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(0, sd = 1.5)
  }
  
  # ilogit for p.
  for (t in 1:(T-1)){
    p_male[t] <- ilogit(beta[1]+ lambda[t])
    p_female[t] <- ilogit(beta[2] + lambda[t])
  }
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})

#Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#Now the data in a list.
my.data <- list(y = y + 1)

#Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  phi = runif(my.constants$T-1,0,1),
                                  lambda = rnorm(my.constants$T-1, 0, 1),
                                  z = zinits)
initial.values()

#Specify the parameters we wish to monitor.
parameters.to.save <- c("phi", "p_male", "p_female")
parameters.to.save

#MCMC details.
n.iter <- 8000
n.burnin <- 1000
n.chains <- 2

#At last, let's run nimble.
mcmc.phitps_t <- nimbleMCMC(code = hmm.phitps_t,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin,
                            nchains = n.chains)

#Examine the results.
MCMCsummary(mcmc.phitps_t, round = 2)
MCMCtrace(mcmc.phitps_t, params = "p_female", pdf=F)


#model 10: phi dependent on sex, p dependent on sex 
#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(s)p(s)--------------
hmm.phisps <- nimbleCode({
  
  # Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
 
   #Survival
  for(i in 1:N){
    logit(phi[i])<- beta[sex[i]]
    
    #Survival matrix
    gamma[1,1,i] <- phi[i] # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i] # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1 # Pr(dead t -> dead t+1)
  }
 
   # Prior for b1
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  
  #ilogit for phi
  phi_male <- ilogit(beta[1])
  phi_female <- ilogit(beta[2])
 
   #Recapture
  for(i in 1:N){
    logit(p[i]) <- beta2[sex[i]]
   
     #Recapture matrix
    omega[1,1,i] <- 1 - p[i] # Pr(alive t -> non-detected t)
    omega[1,2,i] <- p[i] # Pr(alive t -> detected t)
    omega[2,1,i] <- 1 # Pr(dead t -> non-detected t)
    omega[2,2,i] <- 0 # Pr(dead t -> detected t)
  }
  
  # Priors for b2
  beta2[1] ~ dnorm(mean = 0, sd = 1.5)
  beta2[2] ~ dnorm(mean = 0, sd = 1.5)
 
   #ilogit for p
  p_male <- ilogit(beta2[1])
  p_female <- ilogit(beta2[2])
  
#likelihood 
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  beta2 = rnorm(2,0,1),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("beta", "phi_male", "phi_female", "p_male", "p_female")
parameters.to.save

#' MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#' At last, let's run nimble.
mcmc.phisps <- nimbleMCMC(code = hmm.phisps,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.phisps, round = 2)
MCMCtrace(mcmc.phisps, params = "all", pdf=F)

#model 11: phi is dependent on sex, p is dependent on time

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(s)p(t)--------------
hmm.phispt <- nimbleCode({
  
  #Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  for(i in 1:N){
    logit(phi[i])<- beta[sex[i]]
  
      # Survival matrix
    gamma[1,1,i] <- phi[i] # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i] # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1 # Pr(dead t -> dead t+1)
  }
  
  # Prior for b1
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  
  #ilogit for phi
  phi_male <- ilogit(beta[1])
  phi_female <- ilogit(beta[2])
  
  #Recapture
  for(t in 1:(T-1)){
    p[t] ~ dunif(0,1)
    
    #Recapture matrix
    omega[1,1,t] <- 1 - p[t] # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t] # Pr(alive t -> detected t)
    omega[2,1,t] <- 1 # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0 # Pr(dead t -> detected t)
  }
 
   # Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(my.constants$T-1,0,1),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("beta", "phi_male", "phi_female", "p")
parameters.to.save

#' MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#' At last, let's run nimble.
mcmc.phispt <- nimbleMCMC(code = hmm.phispt,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.phispt, round = 2)
MCMCtrace(mcmc.phispt, pdf=F)

#model 12: phi dependent on sex, p dependent on sex+time

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

##### phi(s)p(s+t)--------------
hmm.phisps_t <- nimbleCode({
  
  # Transition matrix
  for(i in 1:N){
    logit(phi[i])<- beta[sex[i]]
  
      #Transition matrix
    gamma[1,1,i] <- phi[i] # Pr(alive t -> alive t+1)
    gamma[1,2,i] <- 1 - phi[i] # Pr(alive t -> dead t+1)
    gamma[2,1,i] <- 0 # Pr(dead t -> alive t+1)
    gamma[2,2,i] <- 1 # Pr(dead t -> dead t+1)
  }
  
  # Priors for b1
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  
  # ilogit for phi
  phi_male <- ilogit(beta[1])
  phi_female <- ilogit(beta[2])
  
  # Initial state prob.
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  ## Observation matrix
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta2[sex[i]] + lambda[t]
      
      #Observation matrix
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
  
  # Priors for b2
  beta2[1] ~ dnorm(mean = 0, sd = 1.5)
  beta2[2] ~ dnorm(mean = 0, sd = 1.5)
  
  # Time fixed effect
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(0, 1.5)
  }
  
  # Recapture probability.
  for(t in 1:(T-1)){
    p_male[t] <- ilogit(beta2[1] + lambda[t])
    p_female[t] <- ilogit(beta2[2] + lambda[t])
  }
  
  #Likelihood 
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  beta2 = rnorm(2,0,1),
                                  lambda = rnorm(51, 0, 1),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("beta", "lambda", "phi_male", "phi_female", "p_male", "p_female")
parameters.to.save

#' MCMC details.
n.iter <- 2500
n.burnin <- 1000
n.chains <- 2

#' At last, let's run nimble.
mcmc.phisps_t <- nimbleMCMC(code = hmm.phisps_t,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin,
                            nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.phisps_t, round = 2)
MCMCtrace(mcmc.phisps_t, params = "p", pdf=F)

#results for graph
results <- MCMCsummary(mcmc.phisps_t, round = 2)
View(results)

#formating data for a line graph
results_2 <- results[ ,-c(3:7)]
View(results_2)
results_p_phi <- results_2 [-c(1:53),]
view(results_p_phi)
results_p <- results_p_phi [-c(103:104),]
view(results_p)

#setting a sex column in results
results_p$sex <- c("female", "male") 
view(results_p)
results_p$sex[1:52]="female"
view(results_p)
results_p$sex[52:102]="male"



#this bit of code here doesn't work and i don't know why!
results_p$year <- c( "1968", "1969", "1970", "1971", "1972", "1973", 
                    "1974", "1975", "1976", "1977", "1978", "1979", "1980",
                    "1981", "1982", "1983", "1984","1985", "1986", "1987", 
                    "1988","1989", "1990", "1991", "1992","1993", "1994", 
                    "1995","1996","1997","1998", "1999", "2000","2001", 
                    "2002", "2003", "2004", "2005", "2006", "2007", "2008",
                    "2009", "2010", "2011", "2012","2013", "2014", "2015", 
                    "2016", "2017", "2018", "1968", "1969", "1970", "1971", "1972", "1973", 
                    "1974", "1975", "1976", "1977", "1978", "1979", "1980",
                    "1981", "1982", "1983", "1984","1985", "1986", "1987", 
                    "1988","1989", "1990", "1991", "1992","1993", "1994", 
                    "1995","1996","1997","1998", "1999", "2000","2001", 
                    "2002", "2003", "2004", "2005", "2006", "2007", "2008",
                    "2009", "2010", "2011", "2012","2013", "2014", "2015", 
                    "2016", "2017", "2018") 
View(results_p)

#creating a line graph for recapture probability over time by sex
ggplot(results_p, aes(x=year, y=mean, group=factor(sex)))+
  #stat_summary(geom = "line", fun.y = "mean", size = 2, aes(group = factor(sex)))+
  geom_point()+
  geom_line(aes(color=sex))+
  scale_color_manual(values = c("#00AFBB", "#FC4E07"))+
  theme(axis.text.x = element_text(size=9, angle=45))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))

#edited by ana
cols <- c("#1170AA",  "#EF6F6A")

myplot=ggplot(results_p, aes(x=year, y=mean, color=factor(sex)))+
  geom_point(position = position_dodge(width = 0.5))+
  scale_color_manual(values = cols)+
  theme(axis.text.x = element_text(size=9, angle=45))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                position = position_dodge(width = 0.5))
  
#theme 
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

#formating data for survival plot
results_phi <- results_p_phi[-c(1:102),]
view(results_phi)
results_phi$sex <- c("female", "male")
View(results_phi)  
results_phi$year <- c()

#colour
cols <- c("#1170AA",  "#EF6F6A")

#creating a graph for survival by sex
myplot_2=ggplot(results_phi, aes( x= sex, y= mean, colour=factor(sex))) +
  geom_line()+
  geom_point()+
  scale_color_manual(values = cols)+
  scale_y_continuous(name="mean", limits=c(0, 1))+
  geom_errorbar(aes(ymin= mean-sd , ymax= mean+sd), width=0.3)

myplot_2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))               
rlang::last_error()
#Survival with interaction of sex and time ~ ϕ(sex·time)
##### phi(s*t)p(.)--------------
#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

#start of code

hmm.phistps <- nimbleCode({
 
   #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #Survival
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- beta[sex[i]] + lambda[t] + kappa[sex[i],t]
     
       #Survival matrix
      gamma[1,1,i,t] <- phi[i,t] # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t] # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0 # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1 # Pr(dead t -> dead t+1)
    }
  }
  
  #Recapture
  p ~ dunif(0, 1) # prior recapture
  
  #Recapture matrix
  omega[1,1] <- 1 - p # Pr(alive t -> non-detected t)
  omega[1,2] <- p # Pr(alive t -> detected t)
  omega[2,1] <- 1 # Pr(dead t -> non-detected t)
  omega[2,2] <- 0 # Pr(dead t -> detected t)

    ## Priors for b1
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
 
   #Time fixed effect
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(mean = 0, sd = 1.5)
  }
 
   # Time as random for the interaction
  t.sigma ~ dunif(0, 10)
  for(i in 1:2){
    for(t in 1:(T-1)){
      kappa[i,t] ~ dnorm(mean = 0, sd = t.sigma)
    }
  }
  
  # ilogit for phi
  for (t in 1:(T-1)){
    phi_male[t] <- ilogit(beta[1]+ lambda[t] + kappa[1,t])
    phi_female[t] <- ilogit(beta[2] + lambda[t] + kappa[2,t])
  }
  
  # Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  p = runif(1,0,1),
                                  lambda = rnorm(51, 0, 1),
                                  t.sigma = runif(1,0,1),
                                  kappa = matrix(rnorm(102, 0, 1), 2, 51),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("phi_male", "phi_female", "p")
parameters.to.save

#' MCMC details.
n.iter <- 10000
n.burnin <- 1000
n.chains <- 2

#' At last, let's run nimble.
mcmc.hmm.phistps <- nimbleMCMC(code = hmm.phistps,
                               constants = my.constants,
                               data = my.data,
                               inits = initial.values,
                               monitors = parameters.to.save,
                               niter = n.iter,
                               nburnin = n.burnin,
                               nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.hmm.phistps, round = 2)
MCMCtrace(mcmc.hmm.phistps, params = "all", pdf=F)

#ϕ(sex·time) p(sex)
###### phi(s*t)p(s)--------------

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

#start model code
hmm.phistps <- nimbleCode({
  #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  
  #survival
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- beta[sex[i]] + lambda[t] + kappa[sex[i],t]
      
      #Survival matrix
      gamma[1,1,i,t] <- phi[i,t] # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t] # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0 # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1 # Pr(dead t -> dead t+1)
    }
  }
  
   #Recapture
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i]) <- beta2[sex[i]]
      
      #Recapture matrix
      omega[1,1,i] <- 1 - p[i] # Pr(alive t -> non-detected t)
      omega[1,2,i] <- p[i] # Pr(alive t -> detected t)
      omega[2,1,i] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i] <- 0 # Pr(dead t -> detected t)
    }
  }
  ## Priors for b1 and b2
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  beta2[1] ~ dnorm(mean = 0, sd = 1.5)
  beta2[2] ~ dnorm(mean = 0, sd = 1.5)
  
  #Time fixed effect
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(mean = 0, sd = 1.5)
  }
 
   # Time as random for the interaction
  t.sigma ~ dunif(0, 10)
  for(i in 1:2){
    for(t in 1:(T-1)){
      kappa[i,t] ~ dnorm(mean = 0, sd = t.sigma)
    }
  }
 
   # ilogit for phi
  for (t in 1:(T-1)){
    phi_male[t] <- ilogit(beta[1]+ lambda[t] + kappa[1,t])
    phi_female[t] <- ilogit(beta[2] + lambda[t] + kappa[2,t])
  }
  
  #ilogit for p
  p_male <- ilogit(beta2[1])
  p_female <- ilogit(beta2[2])
  
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  beta2 = rnorm(2,0,1),
                                  lambda = rnorm(51, 0, 1),
                                  t.sigma = runif(1,0,1),
                                  kappa = matrix(rnorm(102, 0, 1), 2, 51),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("phi_male", "phi_female", "p")
parameters.to.save

#' MCMC details.
n.iter <- 10000
n.burnin <- 1000
n.chains <- 2

#' At last, let's run nimble.
mcmc.hmm.phistps <- nimbleMCMC(code = hmm.phistps,
                               constants = my.constants,
                               data = my.data,
                               inits = initial.values,
                               monitors = parameters.to.save,
                               niter = n.iter,
                               nburnin = n.burnin,
                               nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.hmm.phistps, round = 2)
MCMCtrace(mcmc.hmm.phistps, params = "all", pdf=F)

#ϕ(sex·time) p(sex+time)
##### phi(s*t)p(s+t)--------------

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

#start model code
hmm.phistpst <- nimbleCode({
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- beta[sex[i]] + lambda[t] + kappa[sex[i],t]
      gamma[1,1,i,t] <- phi[i,t] # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t] # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0 # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1 # Pr(dead t -> dead t+1)
    }
  }
  
  #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta2[sex[i]] + lambda2[t]
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
  
  ## Priors for b1 b2 b3 and b4
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  beta2[1] ~ dnorm(mean = 0, sd = 1.5)
  beta2[2] ~ dnorm(mean = 0, sd = 1.5)
 
   #Time fixed effect
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(mean = 0, sd = 1.5)
    lambda2[t] ~ dnorm(mean = 0, sd = 1.5)
  }
  
  # Time as random for the interaction
  t.sigma1 ~ dunif(0, 10)
  for(i in 1:2){
    for(t in 1:(T-1)){
      kappa[i,t] ~ dnorm(mean = 0, sd = t.sigma1)
    }
  }
  
  # ilogit for phi and p
  for (t in 1:(T-1)){
    phi_male[t] <- ilogit(beta[1]+ lambda[t] + kappa[1,t])
    phi_female[t] <- ilogit(beta[2] + lambda[t] + kappa[2,t])
    p_male[t] <- ilogit(beta2[1] + lambda2[t])
    p_female[t] <- ilogit(beta2[2] + lambda2[t])
  }
  
  # Pr(dead t -> detected t)
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})

#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first

#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants

#' Now the data in a list.
my.data <- list(y = y + 1)

#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  beta2 = rnorm(2,0,1),
                                  lambda = rnorm(51, 0, 1),
                                  lambda2 = rnorm(51, 0, 1),
                                  t.sigma1 = runif(1,0,1),
                                  kappa = matrix(rnorm(102, 0, 1), 2, 51),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("phi_male", "phi_female", "p_male", "p_female")
parameters.to.save

#' MCMC details.
n.iter <- 10000
n.burnin <- 1000
n.chains <- 2

#' At last, let's run nimble.
mcmc.phistpst <- nimbleMCMC(code = hmm.phistpst,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin,
                            nchains = n.chains)

#' Examine the results.
MCMCsummary(mcmc.hmm.phistps, round = 2)
MCMCtrace(mcmc.hmm.phistps, params = "all", pdf=F)

#ϕ(sex·time) p(sex·time)
##### phi(s*t)p(s*t)--------------

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

#start model code

hmm.phistpst <- nimbleCode({
  #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
 
   #Survival
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- beta[sex[i]] + lambda[t] + kappa[sex[i],t]
      
      #Survival matrix
      gamma[1,1,i,t] <- phi[i,t] # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t] # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0 # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1 # Pr(dead t -> dead t+1)
    }
  }
  #Recapture
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta2[sex[i]] + lambda2[t] + kappa2[sex[i],t]
      #Recapture matrix
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
  ## Priors for betas
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  beta2[1] ~ dnorm(mean = 0, sd = 1.5)
  beta2[2] ~ dnorm(mean = 0, sd = 1.5)
  #Time fixed effect
  for(t in 1:(T-1)){
    lambda[t] ~ dnorm(mean = 0, sd = 1.5)
    lambda2[t] ~ dnorm(mean = 0, sd = 1.5)
  }
  # Time as random for the interaction
  t.sigma1 ~ dunif(0, 10)
  t.sigma2 ~ dunif(0, 10)
  for(i in 1:2){
    for(t in 1:(T-1)){
      kappa[i,t] ~ dnorm(mean = 0, sd = t.sigma1)
      kappa2[i,t] ~ dnorm(mean = 0, sd = t.sigma2)
    }
  }
  # ilogit for phi and p
  for (t in 1:(T-1)){
    phi_male[t] <- ilogit(beta[1]+ lambda[t] + kappa[1,t])
    phi_female[t] <- ilogit(beta[2] + lambda[t] + kappa[2,t])
    p_male[t] <- ilogit(beta2[1] + lambda2[t] + kappa2[2,t])
    p_female[t] <- ilogit(beta2[2] + lambda2[t] + kappa2[2,t])
  }
  # Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})
#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first
#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants
#' Now the data in a list.
my.data <- list(y = y + 1)
#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  beta2 = rnorm(2,0,1),
                                  lambda = rnorm(51, 0, 1),
                                  lambda2 = rnorm(51, 0, 1),
                                  t.sigma1 = runif(1,0,1),
                                  t.sigma2 = runif(1,0,1),
                                  kappa = matrix(rnorm(102, 0, 1), 2, 51),
                                  kappa2 = matrix(rnorm(102, 0, 1), 2, 51),
                                  z = zinits)
initial.values()

#' Specify the parameters we wish to monitor.
parameters.to.save <- c("phi_male", "phi_female", "p_male", "p_female")
parameters.to.save

#' MCMC details.
n.iter <- 10000
n.burnin <- 1000
n.chains <- 2
#' At last, let's run nimble.
mcmc.phistpst <- nimbleMCMC(code = hmm.phistpst,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin,
                            nchains = n.chains)
#' Examine the results.
MCMCsummary(mcmc.phistpst, round = 2)
MCMCtrace(mcmc.phistpst, params = "all", pdf=F)

##ϕ(sex+time) p(sex+time)
##### phi(s+t)p(s+t)--------------

#removing incorrect and missing sex data
which(Fulmar == 3, arr.ind=TRUE)
which(Fulmar == -99, arr.ind=TRUE)
Fulmar_2 <- na.omit(Fulmar)
Fulmar_2 <- Fulmar %>% slice(-c(121, 129, 137, 142, 145, 146, 147,
                                153, 157, 159, 171, 174, 355, 371, 362, 364))

# Format data for CJS
y <- Fulmar_2 %>%
  select(year_1967:year_2018) %>%
  as.matrix()
head(y)

#defining sex
sex <- Fulmar_2 [, 53]
sex_fulmar<-sex
str(as.vector(sex_fulmar))
sex_fulmar<-unlist(as.vector(sex_fulmar))
names(sex_fulmar)<-NULL
sex<-sex_fulmar

#start model code
hmm.phistpst <- nimbleCode({
  #Initial state prob
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  #Survival
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- beta[sex[i]] + lambda[t]
      #Survival matrix
      gamma[1,1,i,t] <- phi[i,t] # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t] # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0 # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1 # Pr(dead t -> dead t+1)
    }
  }
  #Recapture
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- beta2[sex[i]] + lambda2[t]
      #Recapture matrix
      omega[1,1,i,t] <- 1 - p[i,t] # Pr(alive t -> non-detected t)
      omega[1,2,i,t] <- p[i,t] # Pr(alive t -> detected t)
      omega[2,1,i,t] <- 1 # Pr(dead t -> non-detected t)
      omega[2,2,i,t] <- 0 # Pr(dead t -> detected t)
    }
  }
  ## Priors for b1 b2
  beta[1] ~ dnorm(mean = 0, sd = 1.5)
  beta[2] ~ dnorm(mean = 0, sd = 1.5)
  beta2[1] ~ dnorm(mean = 0, sd = 1.5)
  beta2[2] ~ dnorm(mean = 0, sd = 1.5)
  #Time fixed effect
  for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mean = 0, sd = 1.5)
    lambda2[t] ~ dnorm(mean = 0, sd = 1.5)
  }
  #ilogit for phi and p
  for(t in 1:(T-1)){
    phi_male[t] <- ilogit(beta[1]+ lambda[t])
    phi_female[t] <- ilogit(beta[2] + lambda[t])
    p_male[t] <- ilogit(beta2[1] + lambda2[t])
    p_female[t] <- ilogit(beta2[2] + lambda2[t])
  }
  #Likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, i, j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i, j-1])
    }
  }
})
#' Get the occasion of first capture for all individuals.
first <- apply(y, 1, function(x) min(which(x !=0)))
first
#' A list with constants.
my.constants <- list(N = nrow(y),
                     T = ncol(y),
                     first = first,
                     sex = sex)
my.constants
#' Now the data in a list.
my.data <- list(y = y + 1)
#' Specify initial values.
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(beta = rnorm(2,0,1),
                                  beta2 = rnorm(2,0,1),
                                  lambda = rnorm(51,0,1),
                                  lambda2 = rnorm(51,0,1),
                                  z = zinits)
initial.values()
#' Specify the parameters we wish to monitor.
parameters.to.save <- c("phi_male", "phi_female", "p_male", "p_female")
parameters.to.save
#' MCMC details.
n.iter <- 10000
n.burnin <- 1000
n.chains <- 2
#' At last, let's run nimble.
mcmc.phistpst <- nimbleMCMC(code = hmm.phistpst,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = parameters.to.save,
                            niter = n.iter,
                            nburnin = n.burnin,
                            nchains = n.chains)
#' Examine the results.
MCMCsummary(mcmc.phistpst, round = 2)
MCMCtrace(mcmc.phistpst, params = "all", pdf=F)
