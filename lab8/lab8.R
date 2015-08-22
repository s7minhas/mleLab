# Lab 8
# SM

###################################################
# Set up workspace
rm(list=ls()) 
setwd("~/Dropbox/Duke/Spring 2015/PS 733/lab8")
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
		pi = as.vector(1/(1+exp(-X %*% b))) # the link function
		- sum(y*log(pi) + (1-y)*log(1-pi)) # the neg log likelihood
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
summary(m2)
###################################################
	
###################################################
### Interpretation
## Predicted probabilities

# let's say we are interested in the effect of lwg
# we hold age and inc at their mean
lwgRange = seq(min(data$lwg), max(data$lwg), length.out=100)
scen = with(data=data, cbind(1, mean(age), lwgRange, mean(inc)) )

# Calculate predicted log odd values
predVals = scen %*% coef(m2)

# Apply link function to get predicted probabilities
predProbs = 1/(1 + exp(-predVals))

# Plot result
ggData = data.frame(lwgRange, predProbs)
names(ggData) = c('Wage', 'LFP')
ggplot(ggData, aes(x=Wage, y=LFP)) + geom_line() + ylim(0, 1)

# Add confidence intervals
predScen = data.frame(scen)
names(predScen) = names(coef(m2))
preds = predict(m2, 
	type='link', se.fit=TRUE, 
	newdata=predScen )
critVal = qnorm(0.975)
predInts = data.frame(
	LFP=preds$fit,
	up95=preds$fit + (critVal * preds$se.fit),
	lo95=preds$fit - (critVal * preds$se.fit)
	)

# Convert to predicted probabilities
predInts = apply(predInts, 2, function(x){ 1/(1 + exp(-x)) })

# Now plot with conf int
ggData = data.frame(Wage=lwgRange, predInts)
ggplot(ggData, aes(x=Wage, y=LFP, ymin=lo95, ymax=up95)) + geom_line() + geom_ribbon(alpha=.5)

## Simulation based approach
sims = 1000
draws = mvrnorm(sims, coef(m2), vcov(m2))
predUncert = draws %*% t(scen)

# Convert to predicted probabilities
predUncert = apply(predUncert, 2, function(x){ 1/(1+exp(-x)) })

# Spaghetti plot
ggData = melt(predUncert)
names(ggData) = c('Scenario', 'Wage', 'LFP')
ggData$Wage = rep(lwgRange, each=sims)
ggplot(ggData, aes(x=Wage, y=LFP, group=Scenario)) + geom_line(alpha=.1)

# Alternatively just pull out the mean and 95% interval of the simulations
ggData = t(apply(predUncert, 2, function(x){ 
	mean = mean(x)
	qlo95 = quantile(x, 0.025)
	qhi95 = quantile(x, 0.975)
	rbind(mean, qlo95, qhi95) } ))
ggData = data.frame(ggData)
colnames(ggData) = c('LFP', 'Lo95', 'Hi95')
ggData$Wage = lwgRange

ggplot(ggData, aes(x=Wage, y=LFP, ymin=Lo95, ymax=Hi95)) + geom_line() + geom_ribbon(alpha=.5)
###################################################

###################################################
### Performance

# AIC and BIC
# To calculate either of these we first need to determine
	# the log likelihood. We derived this earlier and found:
	# log-likelihood = sum(y*log(pi) + (1-y)*log(1-pi)), where 
	# pi = 1/(1+exp(-X %*% b))
# So first lets calculate pi:
beta = coef(m2)
X = data.matrix( cbind(1, data[, names(beta)[2:length(beta)] ] ) )
predProbs = 1/( 1 + exp(-X %*% beta) )
# Lets plug and chug now
logLike = sum(data$lfp * log(predProbs) + (1 - data$lfp)*log(1-predProbs))

# Using this we can calculate the AIC and BIC
# AIC: -2 * log-likelihood + 2 * npar
-2 * logLike + 2 * length(coef(m2)) # AIC(m2)
# BIC: -2 * log-likelihood + log(nrow(data)) * npar
-2 * logLike + log(nrow(data)) * length(coef(m2)) # AIC(m2, k=log(nrow(data)))
# The smaller the AIC and BIC the better the fit

## ROC Plot
# Choose a threshold
threshold = 0.5
# Classify predictions as zero or one based on threshold
predRoc = as.numeric(predProbs > threshold)

# Calculate False positive rate and True positive rate
# FPR: probability to be predicted positive, given that someone is negative
# TPR: probability to be predicted positive, given that someone is positive
FPR = sum( (predRoc==1)*(data$lfp==0) ) / sum(data$lfp==0)
TPR = sum( (predRoc==1)*(data$lfp==1) ) / sum(data$lfp==1)
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

roc.curve=function(threshold){
Ps=(predProbs>threshold)*1
FP=sum((Ps==1)*(data$lfp==0))/sum(data$lfp==0)
TP=sum((Ps==1)*(data$lfp==1))/sum(data$lfp==1)
vect=c(FP,TP)
names(vect)=c("FPR","TPR")
return(vect)
}

roc = function(prediction, actual){
	library(ROCR)
	pred = prediction(prediction, actual)
	perf = performance(pred,"tpr","fpr")
	rocData = data.frame(attributes(perf)$x.values[[1]], attributes(perf)$y.values[[1]])
	names(rocData) = c('FPR', 'TPR')
	return(rocData)
}

getAUC = function(prediction, actual){
	library(ROCR)
	pred = prediction(prediction, actual)	
	attributes(performance(pred,"auc"))$y.values[[1]]
}

# Plot roc plot
# rocCurve = roc( seq(0, 1, .0001) )
rocCurve = roc(predProbs, data$lfp)
ggplot(rocCurve, aes(x=FPR, y=TPR)) + geom_line() + geom_abline(intercept=0, slope=1)

auc = getAUC(predProbs, data$lfp)
auc

library(ROCR)

pred <- prediction(predProbs, data$lfp)
perf <- performance(pred,"tpr","fpr")
# plotting the ROC curve
par(mfrow=c(1,1))
plot(perf,col="black",lty=3, lwd=3)
plot(attributes(perf)$x.values[[1]], attributes(perf)$y.values[[1]], type='l')

# Approximating area under the curve
i = 2:nrow(rocCurve)
auc = (rocCurve$FPR[i] - rocCurve$FPR[i - 1]) %*% (rocCurve$TPR[i] + rocCurve$TPR[i - 1])/2
auc

## Alternatively, .... separation plots
sepData = data.frame(prob = predProbs, dv = data$lfp)
sepData = sepData[order(sepData$prob), ]
col = c(rgb(red = 254, green = 232, blue = 200, max = 255), 
	rgb(red = 227, green = 74, blue = 51, max = 255)) 

sepPlot = ggplot(data=sepData, aes(ymin=0, ymax=1, xmin=0, xmax=1:length(dv)))
# Create a rectangular plot to fill
sepPlot = sepPlot + geom_rect(fill = "#FEE8C8")
# Add predicted probability line
sepPlot = sepPlot + geom_line(aes(y=prob, x=1:length(dv)), lwd=.8)
# Add lines to denote events
sepPlot = sepPlot + geom_linerange(aes(color=factor(dv), x=1:length(dv)), alpha=.5)
# Color event lines
sepPlot = sepPlot + scale_color_manual(values = col)
sepPlot = sepPlot + scale_y_continuous("", breaks = c(0, 0.25, 0.5, 0.75, 1.0))
sepPlot = sepPlot + scale_x_continuous("", breaks = NULL)
sepPlot = sepPlot + theme(
	axis.ticks = element_blank(),
	legend.position = "none", 
	panel.background = element_blank(), 
	panel.grid = element_blank() )
sepPlot
###################################################