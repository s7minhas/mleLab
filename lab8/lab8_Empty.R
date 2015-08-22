# Lab 7
# SM

###################################################
# Set up workspace
rm(list=ls()) 
setwd("~/Dropbox/Duke/Spring 2015/PS 733/lab7")
set.seed(6886)

# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
	if(! lib %in% installed.packages()[,1])
	  { install.packages(lib, repos='http://cran.rstudio.com/') }
	suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c('foreign', 'MASS', 'ggplot2', 'reshape2')
loadPkg(packs)
###################################################

###################################################
# Load Data
data = read.dta("http://www.indiana.edu/~jslsoc/stata/spex_data/binlfp2.dta")
	# DV: female labor force participation (lfp) 
	# IVs: age (age), log wage (lwg), household income (inc)
table(data$lfp)
data$lfp = as.numeric(data$lfp)-1
###################################################


###################################################
### BUILDING YOUR OWN MLE ESTIMATOR, STEP BY STEP

# 1. what do we need?
# input - data: y, X 
# input - model components: stochastic/family, systematic/link
# output - vector of betas

binreg = function(y, X){
	X = cbind(1,X)
	neg.loglik = function(b, X, y){
		 pi = 1/(1 + exp(-X %*% b)) # the link function
		 -sum(y * log(pi) + (1-y)*log(1-pi)) # the neg log likelihood
	}
	results = optim(rep(0,ncol(X)), neg.loglik, hessian=T, method="BFGS", X=X, y=y)
	list(coefficients=results$par,
		varcovariance=solve(results$hessian), 
		converged=results$convergence==0)
}

# Run function on data
dv = data$lfp
iv = cbind(data$age, data$lwg, data$inc)
m1 = binreg(dv, iv)	

m1$coefficients
sqrt(diag(m1$varcovariance))
m1$converged
###################################################		
	
###################################################
### GLM function
m2 = glm(lfp ~ age + lwg + inc, data=data, family=binomial(link="logit"))
summary(m2)$'coefficients'
###################################################
	
###################################################
### Interpretation
## Predicted probabilities

# let's say we are interested in the effect of lwg
# we hold age and inc at their mean
lwgRange = seq(min(data$lwg), max(data$lwg), length.out=100)
scen = with(data=data, 
	cbind(1, mean(age), lwgRange, mean(inc))
	)

# Calculate predicted log odd values
predVals = scen %*% coef(m2)

# Apply link function to get predicted probabilities
predProbs = 1/(1 + exp(-predVals))

# Plot result
ggData = data.frame(lwgRange, predProbs)
names(ggData) = c('Wage', 'LFP')
ggplot(data=ggData, aes(x=Wage, y=LFP)) + geom_line() + ylim(0,1)

# Add confidence intervals
predScen = data.frame(scen)
names(predScen) = names(coef(m2))
preds = predict(m2, 
	type='link', se.fit=TRUE, newdata=predScen
	)
critVal = qnorm(0.975)
predInts = data.frame(
	LFP = preds$fit, 
	up95 = preds$fit + critVal*preds$se.fit,
	lo95 = preds$fit - critVal*preds$se.fit
	)

# Convert to predicted probabilities
predInts = apply(predInts, 2, function(x) 1/(1+exp(-x)) )

# Now plot with conf int
ggData = data.frame(Wage=lwgRange, predInts)
ggplot(data=ggData, aes(x=Wage, y=LFP, ymin=lo95, ymax=up95)) + 
	geom_line() + geom_ribbon(alpha=.5) + ylim(0,1)

## Simulation based approach
sims = 1000
draws = mvrnorm(sims, coef(m2), vcov(m2))
predUncert = draws %*% t(scen)

# Convert to predicted probabilities
predUncert = apply(predUncert, 2, function(x) 1/(1+exp(-x)) )

# Spaghetti plot
ggData = melt(predUncert)
names(ggData) = c('Scenario', 'Wage', 'LFP')
ggData$Wage = rep(lwgRange, each=sims)
ggplot(ggData, aes(x=Wage, y=LFP, group=Scenario)) + geom_line(alpha=.1)

# Alternatively just pull out the mean and 95% interval of the simulations



###################################################

###################################################
# AIC and BIC
# To calculate either of these we first need to determine
	# the log likelihood. We derived this earlier and found:
	# log-likelihood = sum(y*log(pi) + (1-y)*log(1-pi)), where 
	# pi = 1/(1+exp(-X %*% b))
# So first lets calculate pi:
beta = coef(m2)
X = data.matrix(cbind(1, data[,c('age', 'lwg', 'inc')]))
predProbs = 1/(1+ exp(-X %*% beta))

# Lets plug and chug now
logLike = sum(data$lfp * log(predProbs) + (1-data$lfp)*log(1-predProbs))

# Using this we can calculate the AIC and BIC
# AIC: -2 * log-likelihood + 2 * npar
-2*logLike + 2 * 4
AIC(m2)

# BIC: -2 * log-likelihood + log(nrow(data)) * npar
-2*logLike + log(nrow(data)) * 4
AIC(m2, k=log(nrow(data)))

# The smaller the AIC and BIC the better the fit

## ROC Plot
# Choose a threshold
threshold = 0.5

# Classify predictions as zero or one based on threshold
predRoc = as.numeric( predProbs > threshold )

# Calculate False positive rate and True positive rate
# FPR: probability to be predicted positive, given that someone is negative
# TPR: probability to be predicted positive, given that someone is positive
FPR = sum( (predRoc==1)*(data$lfp==0) ) / sum(data$lfp == 0)
TPR = sum( (predRoc==1)*(data$lfp==1) ) / sum(data$lfp == 1)
FPR; TPR

# To generate a ROC plot we do this process for multiple threshold values
roc = function(threshold){
	# Set up output matrix
	pr=matrix(NA, ncol=3, nrow=length(threshold), 
		dimnames=list(NULL, c('Threshold', 'FPR', 'TPR') ) )

	# Loop through thresholds
	for(ii in 1:length(threshold)){
		predRoc = as.numeric( predProbs > threshold[ii] )
		FPR = sum( (predRoc==1)*(data$lfp==0) ) / sum(data$lfp==0)
		TPR = sum( (predRoc==1)*(data$lfp==1) ) / sum(data$lfp==1)
		pr[ii,1] = threshold[ii]
		pr[ii,2] = FPR
		pr[ii,3] = TPR
	}

	# Return output
	return( data.frame( pr ) )
}

# Plot roc plot
rocCurve = roc( seq(0, 1, by=.001) )
ggplot(rocCurve, aes(x=FPR, y=TPR)) + geom_line() + geom_abline(intercept=0, slope=1)

# Approximating area under the curve
i = 2:nrow(rocCurve)
auc = (rocCurve$FPR[i] - rocCurve$FPR[i - 1]) %*% (rocCurve$TPR[i] + rocCurve$TPR[i - 1])/2
auc

## Alternatively, .... separation plots
library(separationplot)
separationplot(pred=as.vector(predProbs), actual=as.vector(data$lfp), 
	newplot=FALSE)
###################################################




