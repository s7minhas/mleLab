# Ordinal, Nominal, Conditional Logit Lab
# SM

###################################################
# Set up workspace
rm(list=ls()) 
setwd("~/Dropbox/Duke/Spring 2015/PS 733/lab10")
set.seed(6886)

# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
	if(! lib %in% installed.packages()[,1])
	  { install.packages(lib, repos='http://cran.rstudio.com/') }
	suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c('foreign', 'MASS', 'nnet', 'separationplot',
	'ggplot2', 'reshape2', 'rms')
loadPkg(packs)
theme_set(theme_bw())

# Helpful functions
char = function(x){ as.character(x) }
###################################################

###################################################
#### Ordered Logistic Regression
# Load data
krain<-read.dta("isq05.dta",convert.dates = TRUE,
	convert.factors = TRUE, missing.type = FALSE,
	convert.underscore=TRUE, warn.missing.labels=TRUE)

# Condense categories in magnitud vector
krain$mag2 = krain$magnitud
krain$mag2[krain$magnitud < 2] = 1
krain$mag2[krain$magnitud >= 2 & krain$magnitud < 4] = 2
krain$mag2[krain$magnitud >= 4] = 3
table(krain$mag2, krain$magnitud)

# Convert to factor
krain$mag2 = factor(krain$mag2, levels = sort( unique(krain$mag2) ) )

# Model parameters
dv = 'mag2'
ivs = c('intrvlag', 'icntglag', 'maglag', 'genyr', 
	'stfl', 'regtype', 'ethkrain', 'marg', 'coldwar')
form = formula(paste0(
	dv, ' ~ ', paste0(ivs, collapse=' + ')
	) )

# Run model
ordMod = polr(form, data=krain, method="logistic")

# Model output
summary(ordMod)

# Extract parameter estimates
beta = coef(ordMod)
tau = ordMod$zeta

# Set up scenario, vary the level of ethnic fractionalization
scenRange = with(krain, seq(min(ethkrain), max(ethkrain), by=.01) )
X = with(krain, 
	cbind( median(intrvlag), median(icntglag), mean(maglag), mean(genyr),
		median(stfl), mean(regtype), scenRange, mean(marg), 1
	 ) )
colnames(X) = names( beta )

# Calculate predicted probabilities
probLow = plogis(tau[1] - X %*% beta)
probMed = plogis(tau[2] - X %*% beta) - plogis(tau[1] - X %*% beta)
probHigh = 1 - plogis(tau[2] - X %*% beta)

# Set up dataframe for plotting
ggData = data.frame( scenRange, probLow, probMed, probHigh )
ggData = melt(ggData, id='scenRange')

# Plot
tmp=ggplot(ggData, aes(x=scenRange, y=value, color=variable)) + geom_line()
tmp=tmp + xlab('Ethnic Fractionalization') + ylab('Probability')
tmp=tmp + theme(legend.position='top', legend.title=element_blank())
tmp

# predicted values w uncertainty
sims=1000
draws = mvrnorm(sims, c(beta, tau), vcov(ordMod))
betaDraws = draws[, 1:length(ivs) ]
tauDraws = draws[, (length(ivs) + 1):ncol(draws) ]
preds = betaDraws %*% t(X)

# Calculate predicted probabilities
probLow = plogis(tauDraws[,1] - preds)
probMed = plogis(tauDraws[,2] - preds) - plogis(tauDraws[,1] - preds)
probHigh = 1 - plogis(tauDraws[,2] - preds)

# Pull out mean and 95% interval
info = function(x){ c( mean(x), quantile(x, probs=c(0.025, 0.975)) ) }
probLowSumm = t( apply(probLow, 2, info) )
probMedSumm = t( apply(probMed, 2, info) )
probHighSumm = t( apply(probHigh, 2, info) )

# Set up dataframe for plotting
ggData = data.frame( rbind( 
	cbind(scenRange, probLowSumm),
	cbind(scenRange, probMedSumm),
	cbind(scenRange, probHighSumm)
 ) )
colnames(ggData) = c('scenRange', 'mu', 'lo', 'hi')
ggData$scen = rep( c('probLow', 'probMed', 'probHigh'), each=nrow(X) )

# Plot
tmp=ggplot(ggData, aes(x=scenRange, y=mu, ymin=lo, ymax=hi, color=scen, fill=scen))
tmp=tmp + geom_ribbon(alpha=.2) + geom_line()
tmp

# Parallel lines assumption
	# Ordinal logistic regression assumes that the coefficients 
	# are equal across cuts of the dependent variable.
	# Because the relationship between all pairs of groups is the same, there is only 
	# one set of coefficients. If this was not the case, we would need different sets of 
	# coefficients in the model to describe the relationship between each pair of outcome groups. 

# Ward, Ahlquist pg. 107
	# If the included regressors are able to differentiate one category 
	# from another then we should expect to see a strong trend across the levels 
	# of the dependent variable. If the parallel regressions assumption holds then 
	# this trend should be linear and the conditional means should line up neatly 
	# along the dotted trend line.
par(mfrow=c(3,3))
plot.xmean.ordinaly(form, data=krain)
par(mfrow=c(1,1))

# Performance
sp.categorical(pred=ordMod$fitted.values, actual=krain$mag2)
###################################################

###################################################
# Multinomial Logistic
	# The multinomial logit are designed to model nominal outcomes in such a way
	# that the effects of the independent variable vary for each outcome.

# Can use this to test parallel lines assumption
multMod = multinom(form, data=krain)
summary(multMod)
# Ward, Ahlquist pg. 128
	# We can use a likelihood-based test by comparing the deviance in each model. 
	# The difference should be distributed χ2 with degrees of freedom equal to the 
	# difference in the degrees of freedom in each model
pchisq(deviance(ordMod) - deviance(multMod), multMod$edf - ordMod$edf, lower.tail=FALSE)

# Multinomial logit analysis
datahsb = read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")

# Convert to factor
datahsb$female = char(datahsb$female)
datahsb$female[datahsb$female=='female']=1
datahsb$female[datahsb$female=='male']=0
datahsb$female = factor(datahsb$female)

# Set parameters for model
dv='prog'
ivs=c('female', 'math')
form = formula(paste0(
	dv, ' ~ ', paste0(ivs, collapse=' + ')
	) )

# Run model
	# DV category left out is general
multMod = multinom(form, data=datahsb, Hess = TRUE)

# Get relevant parameter estimates
coefs = t( summary(multMod)$coefficients )
serrors = t( summary(multMod)$standard.errors )
zscores = coefs/serrors
pvals = (1 - pnorm(abs(zscores), 0, 1)) * 2

# Lets do some substantive effects analysis
beta = multMod$wts[multMod$wts!=0]

# Set up scenario
mathRange = with(datahsb, seq( min(math), max(math), by=3) ) 
xFemale = cbind(1, 1, mathRange)
xMale = cbind(1, 0, mathRange)
X = rbind(xFemale, xMale)

# Pull some draws
sims=1000 
set.seed(6886); draws = mvrnorm(sims, beta, solve(multMod$Hess) )

# Calculate parameters
acadPred = exp( draws[ ,1:( length(ivs) + 1 ) ] %*% t(X) )
vocPred = exp( draws[ ,( length(ivs) + 2 ):ncol(draws) ] %*% t(X) )

# Get denominator
simDenom = (1 + acadPred + vocPred ) 

# Get simulated probabilities for each category
genProb = 1/simDenom 
acadProb = acadPred / simDenom
vocProb = vocPred / simDenom

# Pull out mean and 95% interval
info = function(x){ c( mean(x), quantile(x, probs=c(0.025, 0.975)) )  }
genProbSumm = t( apply(genProb, 2, info) )
acadProbSumm = t( apply(acadProb, 2, info) )
vocProbSumm = t( apply(vocProb, 2, info) )

# Set up dataframe for plotting
ggData = data.frame( rbind( 
	cbind(mathRange, rep(c(1,0), each=15), genProbSumm),
	cbind(mathRange, rep(c(1,0), each=15), acadProbSumm),
	cbind(mathRange, rep(c(1,0), each=15), vocProbSumm)
 ) )
colnames(ggData) = c('mathRange', 'female', 'mu', 'lo', 'hi')
ggData$scen = rep( c('General', 'Academic', 'Vocational'), each=nrow(X) )

# Plot
tmp=ggplot(ggData, aes(x=mathRange, y=mu, ymin=lo, ymax=hi, color=scen, fill=scen))
tmp=tmp + geom_ribbon(alpha=.2) + geom_line()
tmp=tmp + facet_wrap(~female, scales='free') + theme(legend.position='top')
tmp

# Performance
sp.categorical(pred=multMod$fitted.values, actual=datahsb$prog)
###################################################

###################################################
# Conditional logit
# dta file
mycnd = read.dta('CES 1988-2011.dta')
# subset to 2000 & Pull out relevant vars

###################################################