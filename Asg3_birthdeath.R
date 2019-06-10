######### The log-likelihood of the birthdeath model #########
## This function calculates log-likelihood of a tree under a birthdeath model given the input parameters 
## The function takes a birth rate (birthRate), a death rate (deathRate) and a phylogenetic tree (an object class phylo) 
## and outputs log-likelihood 
## It is extracted from the likelihood equation (dev function) in ape::birthdeath
## and is a log form of EQUATION 21, Nee et al. (1994). "The reconstructed evolutionary
## process." Philosophical Transactions of the Royal Society B, Biological Sciences. 344, 
## 305-311. doi: 10.1098/rstb.1994.0068

# R package
library(ape)

# read in data
phy <- read.nexus("https://www.r-phylo.org/w/images/0/02/Geospiza.nex")
# plot phylogenetic tree
plot(phy)
title("phylogeny of Darwin's finches (from APE package)") # add title
axisPhylo()
mtext("millions of years ago (Ma)", side=1, line=2.5) # add axis

# Create a function called bd_loglik()
bd_loglik <- function(birthRate, deathRate, phy)
  {
  # r is diversification rate (r = birth rate - death rate)
  r <- birthRate - deathRate
  # a is relative extriction rate (a = death rate/ birth rate)
  a <- deathRate / birthRate
  # input phylogenetic tree
  phy <- phy
  #	x is the length of time interval between the present and the birth of nth lineage (age of internal nodes)
  # because we start from 2th lineage so we have NA for the first 
  x <- c(NA, branching.times(phy))
  # N is the number of lineages in a tree
  N <- length(phy$tip.label)
  # calculate the log likelihood of birthdeath model 
  # log(gamma(N)) is n-1! possible topologies for any set of n-1 waiting times
  # waiting time is the time interval between events (successive speciation)
  log_L = (lfactorial(N - 1) 
           # log(diversification rate) = expected waiting time
           # N-2 expected waiting times for 2 internal branches 
           + (N - 2) * log(r) 
           # obtain expected number of diversification on 2 internal branches 
           # by multiplying diversification rate with sum of x (age of internal nodes)(excluding root)
           + r * sum(x[3:N]) + 
             # N is the number of tips (lineages that do not give births) 
             # N * log(1 - a) is waiting time of no events
             N * log(1 - a) 
           # Divide by two lineages times amount of waiting time from each node (including root) 
           # exp(r * x[2:N]) = expected number of diversification from each node
           # and sum(log(exp(r * x[2:N]) - a)) = total amount of waiting time from each node
           - 2 * sum(log(exp(r * x[2:N]) - a)))
  return(log_L)
  }
# Note: The likelihood of this model (for a tree with exant tips) is the likelihood that each lineage speciate 
# times the likelihood that lineages do not give birth 
# times the likelihood that the first 2 lineages are not extinct at the present (Nee et al., 1994)
# putting together, we can derive the likelihood from multiplying together the waiting times of speciation events in a tree with parameters (a and r)

# testing bd_loglik function
bd_loglik(birthRate=0.5, deathRate=0.3, phy=phy)
bd_loglik(birthRate=0.5, deathRate=0.2, phy=phy)
bd_loglik(birthRate=0.5, deathRate=0.1, phy=phy)
bd_loglik(birthRate=0.5, deathRate=0.0, phy=phy)


# plot log-likelihood curve with a fixed value of death rate
death_rates <- rep(1.7 , times=100)
birth_rates <- seq(from=1, to=10, length.out=100)
Log.likelihood = rep(NA, 100)
for (i in 1:100)
{
  Log.likelihood[i] <- bd_loglik(birthRate=birth_rates[i], deathRate=death_rates[i], phy=phy)
}
plot(birth_rates, Log.likelihood, main = "Log-likelihood vs. birth rate", xlab = "birthrates (death rate of 1.7)", ylab = "Log-likelihood")


# plot likelihood curve with a fixed value of death rate
Likelihood = rep(NA, 100)
for (i in 1:100)
{
  Likelihood[i] <- exp(bd_loglik(birthRate=birth_rates[i], deathRate=death_rates[i], phy=phy))
}
plot(birth_rates, Likelihood, main = "Likelihood vs. birth rate", xlab = "birthrates (death rate of 1.7)", ylab = "Likelihood")
abline(h= 253249654, col = "blue") # draw a line for 95% cutoff

# Test if the output from bd_loglik (given the parameters) matches the output of the function birthdeath() from the R package ape 
# The birthdeath() outputs maximum log-likelihood and parameter estimates of relative extinction rate (d / b) and diversification rate (b - d)
# we can obtain a birth rate and death rate from these parameter estimates 
ape_bd <- birthdeath(phy)
names(ape_bd) #get names of birthdeath() output
ape_bd$para["d/b"]

# get diversification and relative extinction rates from the birthdeath function
diversification_rate <- ape_bd$para["b-d"]
extinction_rate <- ape_bd$para["d/b"] 
# converting diversification and relative extrinction rate to birth rate and death rate 
b_rate_ape <- unname(diversification_rate / (1-extinction_rate))
d_rate_ape <- unname(extinction_rate * b_rate_ape)
# use these parameters to calculate log-likelihood  
My_Loglik <- bd_loglik(birthRate=b_rate_ape, deathRate=d_rate_ape, phy=phy)
# see if both functions give the same log-likelihood 
# The maximum likelihood of the birthdeath function is dev/-2
ape_bd$dev/-2
My_Loglik
ape_bd$dev/-2 == My_Loglik # same log-likelihood!

### model comparison ###
## Yule model ##
# Log-likelihood of yule model 
# The function takes a phylogenetic tree and a birth rate and outputs the log-likelihood under the Yule model
# It is extracted from the likelihood equation in ape::yule
yule_loglik <- function(birthRate, phy)
{
  # (-birthRate * X) = the probability of no extra speciation event occurred along the branches 
  # The probability of the nth speciation event is proportional to n times birthRate
  # X is sum of the branch lengths
  X <- sum(phy$edge.length)
  # for n lineages, it is lfactorial(phy$Nnode) + (phy$Nnode - 1) * log(birthRate) 
  Log_L <- -birthRate * X + lfactorial(phy$Nnode) + (phy$Nnode - 1) * log(birthRate) 
  return(Log_L)
}

# comparing yule_loglik function with yule() from the ape library
yule_output <- yule(phy)
yule_output
# lambda = 2.61642
# log-likelihood = 22.09385
yule_loglik(birthRate = 2.61642, phy)
# matched log-likelihood!

# Maximizing the likelihood and get parameter estimates for birth death and yule model
yule_optim <- optimise(yule_loglik, phy = phy, interval = c(0,100), maximum = TRUE)
bd_optim <- optim(c(2.5,2), function(x) bd_loglik(x[1],x[2],phy), control = list(fnscale=-1))
