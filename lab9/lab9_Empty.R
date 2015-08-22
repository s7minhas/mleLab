# Lab 8
# SM

###################################################
# Set up workspace
rm(list=ls()) 
setwd("~/Dropbox/Duke/Spring 2015/PS 733/lab9")
set.seed(6886)

# Function to load packages
loadPkg=function(toLoad){
  for(lib in toLoad){
	if(! lib %in% installed.packages()[,1])
	  { install.packages(lib, repos='http://cran.rstudio.com/') }
	suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c('AER', 'MASS', 'pscl', 'ggplot2', 'reshape2')
loadPkg(packs)

# Helpful functions
char = function(x){ as.character(x) }
###################################################

###################################################
# Load Pakistan Protest Data
load('pakProtestData.rda') # Loads object called pakData
###################################################

###################################################
# Examine DV
dv='protest.tALL'

plot( density( pakData[,dv] ) )
# Examine protests over time
plot( pakData$date, pakData[,dv], type='l' )
###################################################

###################################################
# Model count variables using a poisson process
ivs=c(
	'ProxElection', # Proximity to election
	'NY.GDP.PCAP.KD.l1', # GDP per capita in constant 2000 dollars, lagged by 1 month
	'FP.CPI.TOTL.ZG.l1', # Inflation lagged by one month
	'intratension.l1', # Number of conflictual actions and events internal to the government, lagged by 1 month
	'W.centdist.std.protest.tALL.l1' # Spatial temporal lag of dv
	)
form = as.formula( paste0(
	dv, ' ~ ',
	paste(ivs, collapse=' + ') 
	) )

# Poisson
poisMod = glm(form, data=pakData, family=poisson)
summary(poisMod)
exp( coef(poisMod) )

# Check for overdispersion
AER::dispersiontest( poisMod )

# Gelman & Hill test for overdispersion
# standard deviation of the Poisson is equal to the square root of the mean
predVals = data.matrix(cbind(1, pakData[,ivs])) %*% coef(poisMod)
predCnts = exp(predVals)
z=(pakData[,dv] - predCnts)/sqrt(predCnts)
df=nrow(pakData) - (length(ivs) + 1)
dispRat = sum(z^2)/df
print(paste0('Overdispersion ratio is ', round(dispRat,2) ) )
pval = pchisq(sum(z^2), df)
print(paste0('p-value ', round(pval, 4)) )
# p-value is 1, indicating that the probability is essentially zero 
# that a random variable from this chi sq distribution would be as large
# as what we found
###################################################

###################################################
# How to deal with overdispersion

# Manual quasipoisson
# One quick fix is to multiply all regression standard errors by
# square root of dispersion parameter
dispAdj = sqrt(dispRat)
coefTable = summary( poisMod )$'coefficient'
coefTable[,2] = coefTable[,2] * dispAdj # Adjust standard errors
coefTable[,3] = coefTable[,1]/coefTable[,2] # Recalculate z-statistic
coefTable[,4] = 2*pnorm( -abs(coefTable[,3]) ) # Recalculate p values
coefTable

# Quasilikelihood
qpoisMod = glm(form, data=pakData, family=quasipoisson)
print( round( summary(qpoisMod)$'coefficient', 3) )
print( round( coefTable, 3) )

# Negative binomial model
negbinMod = MASS::glm.nb(form, data=pakData)
print( round( summary(negbinMod)$'coefficient', 3) )
###################################################

###################################################
# One limitation of standard count models is that the zeros 
# and the nonzeros (positives) are assumed to come from the 
# same data-generating process.
# This assumption can be removed using hurdle and zero-inflated models. 

# Introducing some more zeros
sum(pakData[,dv]==0)
pakData2 = pakData
pakData2[ which(pakData2[,dv] < 4), dv] = 0
sum(pakData2[,dv]==0)

# Hurdle & Zero-inflated models
twoForm = as.formula( paste0(
	dv, ' ~ ',
	paste(ivs[1:3], collapse=' + '), ' | ', 
	paste(ivs[4:5], collapse=' + ')	
	) )

# Fitting these types of models
hpoisMod = hurdle(twoForm, data=pakData, dist='poisson')
summary(hpoisMod)$'coefficient'

zpoisMod = zeroinfl(twoForm, data=pakData, dist='poisson')
summary(zpoisMod)$'coefficient'
###################################################