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
packs=c('foreign', 'lmtest', 'sandwich', 'Amelia', 'sbgcop')
loadPkg(packs)
###################################################

###################################################
# Matrix algebra in R
A = matrix(c(2,4,6,8),ncol=2,byrow=T); A
B = matrix(c(1,3,5,7),ncol=2,byrow=T); B

# Multiplying a matrix with scalar
A; 5 * A

# Matrix multiplication
C = A %*% B; C

2*1 + 4*5 # First row, 	First column in C
2*3 + 4*7 # First row, 	Second column in C
6*1 + 8*5 # Second row, First column in C
6*3 + 8*7 # Second row, Second column in C

# Matrix multiplication is not always commutative
D = B %*% A; D 

# Invering a matrix
Ainv = solve(A); Ainv
I = Ainv %*% A; round(I)

# Matrix Multiplication with Identity Matrix
B %*% I
I %*% B
###################################################

###################################################
# Real world example
apsr = read.dta("apsrfinaldata.dta")

data = apsr[,c('gini_net_std', 'polity2', 'ELF_ethnic')]
n=nrow(data)
set.seed(6886)

# Run regression with LM
form = formula(pg ~ betweenstd + lngdpstd)
lm1 = lm(form, data=apsr)
summary(lm1)
###################################################

###################################################
# Now lets do this manually

###################
# Obtain datasets from inputs

# Retrieve vars from formula input
dv = all.vars(form)[1]
ivs = all.vars(form)[ 2:length(all.vars(form)) ]

# Create matrix with column for intercept and
## data from independent variables
y = apsr[,dv]
x = data.matrix(cbind(1, apsr[,ivs]))

# General parameters
n = nrow(x) # Number of observations
p = length(ivs) # Number of parameters
df = n-p-1 # degrees of freedom
###################

###################
# Calculate coefficient estimates
# Beta = (X'X)^-1 %*% X'Y
# calculating X'X
xTx = t(x) %*% x

# calculating X'y
xTy = t(x) %*% y

# calculating Beta
beta = solve(xTx) %*% xTy
beta
###################

###################
# Calculating standard errors
# se(Beta) = sqrt( sigma^2 * (X'X)^-1 )

# First lets calculate our yhat
yhat = x %*% beta

# First get out residuals
e = y - yhat

# calculating e'e (sum of squared residuals): assuming homoskedasticity
eTe = sum( e^2 )
eTe

# Adjust by degrees of freedom
sigma2 = eTe/(df)
sigma2

# Pull out the standard errors for the coefficient estimates
seBeta = sqrt( sigma2 * diag( solve(xTx) ) )
seBeta

# Calculate t values
tStat = beta/seBeta
tStat

# Calculate p-values
pval = 2 * pt(abs(tStat), df = df, lower.tail = FALSE)
pval

# Compare against lm
summary(lm1)$'coefficients'

# Also...the variance-covariance matrix of beta
varcovBeta = sigma2 * solve(xTx)
varcovBeta
vcov(lm1) # Compare with lm output
###################

###################
# Robust standard errors
# se(Beta) = sqrt( (X'X)^-1 ( X' Omega X ) (X'X)^-1 )

# Note that when we calculated eTe above we are essentially assuming
## that we have constant error variance, however, what if our model
## suffers from heteroskedasticity

# Lets check
plot( yhat, e ) # There does appear to be some fan-shaped pattern

# One way to deal with this is to use White robust standard errors
omega = diag( diag( e %*% t(e) ) ) # independent, distinct variances
meat = t(x) %*% omega %*% x # X' Omega X

# Lets pullout the robust variance-covariance matrix
robvarcovBeta = solve(xTx) %*% meat %*% solve(xTx)
robvarcovBeta

# Lets pull out the standard errors with a degrees of freedom adjustment
dfc = n/(df)
rseBeta = sqrt(dfc * diag( robvarcovBeta ) )
rseBeta

# This is equivalent to:
coeftest(lm1, vcov=function(x) vcovHC(x, method="white1", type="HC1"))

# Coefficient p values...same as before
rtStat = beta/rseBeta
rtStat
rPval = 2 * pt(abs(rtStat), df = df, lower.tail = FALSE)
rPval
###################

###################
# R squared
ssReg = sum( (yhat - mean(y))^2 )
ssTot = sum( (y - mean(y))^2 )
Rsq = ssReg/ssTot
Rsq

# Adj r squared
adjRsq = 1-(1-Rsq) * ( (n-1) / (df) )
adjRsq
###################

###################
# F statistic
muSqReg = sum(yhat^2)/p
muSqRes = sum(e^2)/( df )
Fstat = muSqReg/muSqRes
Fstat
pf(Fstat, p, df, lower.tail=FALSE)
###################

# Combine info
coefRes=cbind(beta, seBeta, tStat, pval)
coefRes
# Compare to 
summary(lm1)$'coefficients'

# How to incorporate other info?
###################################################

###################################################
# Missing data

# Lets introduce some random missingness to our dataset
set.seed(6886)
missMatrix = cbind(rbinom(n, 1, .9),rbinom(n, 1, .85),rbinom(n, 1, .83))
missMatrix[missMatrix==0]=NA
apsrMiss = apsr[,c(dv, ivs)] * missMatrix

# Check out new dimensions
dim(na.omit(apsrMiss))

# Simply throwing away missing data is increasingly
## being recognized in our field as bad practice

# Lets first rerun our regression without the missing values
lm1ListDel = lm(form, data=apsrMiss)
round(summary(lm1ListDel)$'coefficients',3)

# Compare with earlier results
round(summary(lm1)$'coefficients',3)

# Most commonly used package in our field for missing
## data is Amelia
apsrAmelia = amelia(x=apsrMiss, m=1)
names(apsrAmelia$imp)

# lets just use the second imputed dataset for now
## and reestimate our model
lm1Amelia = lm(form, data=apsrAmelia$imp$imp1)
round(summary(lm1Amelia)$'coefficients',3)

# Another tool for imputing missing data is sbgcop
# sbgcop is a bayesian model and thus we sample our way to 
## posterior values. the nsamp parameter controls how many
## samples we draw, here I just set it to 5000, but you might 
## need to set it higher, so that you can be sure the results have converged
apsrSbgcop = sbgcop.mcmc(Y=apsrMiss, nsamp=5000, seed=6886)

# apsrSbgcop is a list comprised of a number of objects, the most relevant
names(apsrSbgcop)
## for us is Y.impute, in this case it has dimensions: 46 x 3 x 1000
## the first two dims correspond to the size of our original dataset
## and the last dim, 1000, represents the iterations from our sampler
dim(apsrSbgcop$Y.impute)

# We need to account for the fact that it took some time for our sampler
## to converge, so we need to burn some of our initial draws. below i burn
## the first 50% of the draws we pulled from our bayesian model
toKeep = (dim(apsrSbgcop$Y.impute)[3]/2 + 1):dim(apsrSbgcop$Y.impute)[3]
apsrSbgcopPost = apsrSbgcop$Y.impute[,,toKeep]

# Next we average across these results so that we can have just one 
## consolidated imputed dataset for our analysis
apsrSbgcopPostAvg = apply(apsrSbgcopPost, c(1,2), mean)

# Now lets add colnames back in and turn this back into a dataframe
colnames(apsrSbgcopPostAvg) = names(apsrMiss)
apsrSbgcopPostAvg = data.frame(apsrSbgcopPostAvg)

# Last lets rerun our model using the imputed data from sbgcop
lm1Sbgcop = lm(form, data=apsrSbgcopPostAvg)
round(summary(lm1Sbgcop)$'coefficients',3)

# Compare three sets of results: 

# True model, i.e., no randomly excluded data
round(summary(lm1)$'coefficients',3)

# Listwise deletion
round(summary(lm1ListDel)$'coefficients',3)

# Amelia
round(summary(lm1Amelia)$'coefficients',3)

# Sbgcop
round(summary(lm1Sbgcop)$'coefficients',3)
###################################################