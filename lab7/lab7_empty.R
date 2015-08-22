# Lab 7
# SM

###################################################
# Set up workspace
rm(list=ls()) 
setwd("~/Dropbox/Duke/Spring 2015/PS 733/lab7")

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
A * (1/5)

# Matrix multiplication
A %*% B

# Matrix multiplication is not always commutative
B %*% A

# Invering a matrix
solve(A)

# Matrix Multiplication with Identity Matrix
round(solve(A) %*% A)

###################################################

###################################################
# Real world example
apsr = read.dta("apsrfinaldata.dta")

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


# First get out residuals


# calculating e'e (sum of squared residuals): assuming homoskedasticity


# Adjust by degrees of freedom


# Pull out the standard errors for the coefficient estimates


# Calculate t values


# Calculate p-values


# Compare against lm
summary(lm1)$'coefficients'

# Also...the variance-covariance matrix of beta

vcov(lm1) # Compare with lm output
###################

###################
# Robust standard errors
# se(Beta) = sqrt( (X'X)^-1 ( X' Omega X ) (X'X)^-1 )

# Note that when we calculated eTe above we are essentially assuming
## that we have constant error variance, however, what if our model
## suffers from heteroskedasticity

# Lets check
 # There does appear to be some fan-shaped pattern

# One way to deal with this is to use White robust standard errors
 # independent, distinct variances
 # X' Omega X

# Lets pullout the robust variance-covariance matrix



# Lets pull out the standard errors with a degrees of freedom adjustment



# This is equivalent to:
coeftest(lm1, vcov=function(x) vcovHC(x, method="white1", type="HC1"))

# Coefficient p values...same as before

###################

###################
# R squared


# Adj r squared

###################

###################
# F statistic

###################

# Combine info


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


# lets just use the second imputed dataset for now
## and reestimate our model


# Regression on imputed data is closer to actual result
###################################################