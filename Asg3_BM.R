######### The log likelihood of the Brownian motion (BM) model #########
# This is a modified function from the function BMlk from http://www.phytools.org/***SanJuan2016/ex/5/Fitting-BM.html
# the log-likelihood equation is based on EQ2 on the same page
# This function takes an object class phylo, two parameters (root.state and variance) and trait values 
# and outputs log-likelihood under a Brownian motion model 
# Note: under the Brownian motion model, the trait evolution over time can be modelled as 
# E[z(t)] ~ N(z(0),sigma^2*t1) where z(0) is the root state, t is branch length and ??^2 is rate parameter
# (eq. 3.17) from Harmon's Phylogenetic Comparative Methods: https://lukejharmon.github.io/pcm/chapter3_bmintro/, chapter 3.4

# packages requires 
library(phytools)
library(ape)
library(geiger)

# data 
# read in the Geospiza phylogeny and tip data 
geotree <- read.nexus("https://www.r-phylo.org/w/images/0/02/Geospiza.nex")
geodata <- read.table("https://www.r-phylo.org/w/images/5/5c/Geospiza.txt")
# drop "olivacea" from the analysis
geotree <- drop.tip(geotree, "olivacea")
geodata <- geodata[geotree$tip.label, ]
# use wingL as trait value (data)
wingL <- geodata$wingL
names(wingL) <- row.names(geodata)

# BM loglikelihood function 
BM_loglik <- function(tree, variance, root.state, data) {
  
  # C, the phylogenetic variance-covariance matrix 
  # The calculation of this matrix is based on the branch length and tree topology
  # we use vcvPhylo() to obtain C of tree 
  C <- vcvPhylo(tree, anc.nodes = F)
  inv.C <- solve(C) # inverted C 
  N <- length(data); # the number of tips
  # creates a vector of the expected trait value  
  # this is equivalent to the root state in BM model 
  Expected_trait <- rep(root.state, N) 
  # Rate is rate matrix representing rate of evolution 
  # multipling C by variance
  Rate <- C * variance; 
  # also convert inv.C to the inverse of the rate
  inv.Rate <- inv.C * variance ^-1; 
  # calculate likelihood 
  # The BM likelihood function uses a standard formula for the likelihood of a multivariate normal distribution
  # lnlNum represents the distance of data from the mean, which is the root state in BM model  
  lnlNum<- -0.5*(data - Expected_trait) %*% inv.Rate %*% (data - Expected_trait) # %*% is matrix multiplication
  # det(Rate) calculates the determinant of Rate
  lnlDen<- log(sqrt((2*pi)^N*det(Rate)))
  Loglik<-lnlNum-lnlDen
  return(Loglik);
}

# testing BM_loglik function
BM_loglik(geotree, variance=0.07, root.state = 4.2, wingL)
BM_loglik(geotree, variance=0.06, root.state = 4.2, wingL)
BM_loglik(geotree, variance=0.05, root.state = 4.2, wingL)
BM_loglik(geotree, variance=0.04, root.state = 4.2, wingL)
BM_loglik(geotree, variance=0.03, root.state = 4.2, wingL)

# plot log likelihood curve with fixed value of root state
root.states <- rep(4.2 , times=100)
variances <- seq(from=0.01, to=0.25, length.out=100)

Log.likelihood.BM = rep(NA, 100)
for (i in 1:100)
{
  Log.likelihood.BM[i] <- BM_loglik(geotree, variance=variances[i], root.state=root.states[i], wingL)
}

plot(variances, Log.likelihood.BM, main = "Log-likelihood vs. variance (BM model)", xlab = "variances (root state of 4.2)", ylab = "Log-likelihood")

# plot likelihood of death rates with fixed values of birthrate
# convert log-likelihood to likelihood using exp()  
root.states <- rep(4.2 , times=100)
variances <- seq(from=0.01, to=0.25, length.out=100)

Likelihood.BM <- rep(NA, 100)
for (i in 1:100)
{
  Likelihood.BM[i] <- exp(BM_loglik(geotree, variance=variances[i], root.state=root.states[i], wingL))
}

plot(variances, Likelihood.BM, main = "likelihood vs. variance (BM model)", xlab = "variances (root state of 4.2)", ylab = "Likelihood")
abline(h = 189.6986, col = "blue") # draw a line on 95% cutoff

# testing BM_loglik function
# we will use wingL as data (trait values) and geotree 
# compare with fitContinuous function from the R package geiger  
# the function takes an object class phylo and trait values
# and output maximum likelihood and parameter estimates based on selected model
# find ML and estimate parameters under BM model 
bm.ml <- fitContinuous(phy= geotree, wingL, model="BM")
# bm.ml, sigsq = 0.070546, z0 = 4.205953, log-likelihood = 8.243269
# assign these parameters to BM_loglik function
my_BM_loglik <- as.numeric(BM_loglik(geotree, variance = 0.070546, root.state = 4.205953, wingL))
# see if both functions give the same log-likelihood
my_BM_loglik 
bm.ml$opt$lnL # same log-likelihood


###### model comparison ######
# This function is a modified version of BM_loglik function 
# It allows a transformation of the branch lengths of the tree under three Pagel tree transformations (lambda, kappa and delta)
# this function uses rescale() from the the R package geiger to transform branch lengths 
# This function takes tree (an object class phylo), variance, root.state, trait value (data), Pagel model and transformation parameter
# and outputs log-likelihood under BM model given the parameters  

Pagel_loglik <- function(tree, variance, root.state, data, model = NULL, transform.para) {
  # transforms branch length according to Pagel model 
  # and the branch transformation parameter (transform.para)
  switch(model, lambda = {
    tree <- rescale(tree, "lambda", transform.para)
  }, kappa = {
    tree <- rescale(tree, "kappa", transform.para)
  }, delta = {
    tree <- rescale(tree, "delta", transform.para)
  })
  # create the phylogenetic variance-covariance matrix
  C <- vcvPhylo(tree, anc.nodes = F)
  inv.C <- solve(C) # inverted C 
  N <- length(data); # the number of tips
  # creates a vector of the expected trait value  
  # this is equivalent to the root state in BM model 
  Expected_trait <- rep(root.state, N) 
  # Rate is rate matrix representing rate of evolution 
  Rate <- C * variance; 
  # also convert inverted C to the inverse of the rate
  inv.Rate <- inv.C * variance ^-1; 
  # calculate likelihood 
  # %*% is matrix multiplication
  lnlNum<- -0.5*(data - Expected_trait) %*% inv.Rate %*% (data - Expected_trait)
  lnlDen<- log(sqrt((2*pi)^N*det(Rate)))
  Loglik<-lnlNum-lnlDen
  return(Loglik);
}

# comparing Pagel_loglik with available R function from the R package geiger
lambda.ml <- fitContinuous(phy= geotree, wingL, model="lambda")
# lambda.ml, Lambda = 0.000000, variance = 0.02405703, root.state = 4.235734, ML = 9.805223
# put all parameters in Pagel_loklik
my_lambda <- Pagel_loglik(geotree, variance = 0.022206, root.state = 4.235734, wingL, model = "lambda", 0.000000)
my_lambda # give the same log-likelihood! # ML = 9.805223

# kappa model
kappa.ml <- fitContinuous(phy= geotree, wingL, model="kappa")
# kappa.ml, kappa = 0.610305, variance = 0.027173, root.state = 4.185320, ML = 8.888181
my_kappa <- Pagel_loglik(geotree, variance = 0.027173, root.state = 4.185320, wingL, model = "kappa", 0.610305)
my_kappa # give the same log-likelihood! # ML = 8.888181

# delta model 
delta.ml <- fitContinuous(phy= geotree, wingL, model="delta")
# delta.ml, delta = 2.999999, variance = 0.035363, root.state = 4.227527, ML = 9.362256
my_delta <- Pagel_loglik(geotree, variance = 0.035363, root.state = 4.227527, wingL, model = "delta", 2.999999)
my_delta # give the same log-likelihood! # ML = 9.362256


## Getting ML and parameter estimates for each model using optim function 
bm_optim <- optim(c(0.01,4), control = list(fnscale=-1),
                  method=c("L-BFGS-B"), lower = c(0.01,3.5), upper = c(1,5), 
                  fn = function(x) BM_loglik(geotree, x[1],x[2],data = wingL))

lambda_optim <- optim(c(0.01,4,1e-08), control = list(fnscale=-1),
                      method=c("L-BFGS-B"), lower = c(0.01,3.5,1e-09), upper = c(1,5,1), 
                      fn = function(x) Pagel_loglik(geotree, x[1],x[2],data = wingL, model = "lambda", x[3]))

kappa_optim <- optim(c(0.01,4,0.1), control = list(fnscale=-1),
                     method=c("L-BFGS-B"), lower = c(0.01,3.5,1e-09), upper = c(1,5,1), 
                     fn = function(x) Pagel_loglik(geotree, x[1],x[2],data = wingL, model = "kappa", x[3]))

delta_optim <- optim(c(0.01,4,1), control = list(fnscale=-1),
                     method=c("L-BFGS-B"), lower = c(0.01,3.5,1e-09), upper = c(1,5,20), 
                     fn = function(x) Pagel_loglik(geotree, x[1],x[2],data = wingL, model = "delta", x[3]))




