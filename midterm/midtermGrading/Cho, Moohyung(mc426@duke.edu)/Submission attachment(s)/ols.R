# Set up workspace
# rm(list=ls()) 
# setwd("/Users/moohyungcho/Desktop/Spring 2015/MLE/Lab")

# # Function to load packages
# loadPkg=function(toLoad){
#   for(lib in toLoad){
#     if(! lib %in% installed.packages()[,1])
#     { install.packages(lib, repos='http://cran.rstudio.com/') }
#     suppressMessages( library(lib, character.only=TRUE) ) }
# }

# # Load libraries
# packs=c('foreign', 'lmtest', 'sandwich', 'Amelia', 'sbgcop')
# loadPkg(packs)

##################

# OLS function

ols <- function(formula, data, impute=FALSE){
  # impute missing data if impute=TRUE
  if(impute){
    set.seed(6886)
    data=amelia(x=data, m=1) # m is the number of imputed dataset
    data=data$imp$imp1
  }  
  else{data=na.omit(data)}
  
  # Obtain datasets from inputs
  
  # Retrieve vars from formula input
  dv = all.vars(form)[1]
  ivs = all.vars(form)[ 2:length(all.vars(form)) ]
  
  # Create matrix with column for intercept and
  ## data from independent variables
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  colnames(x)[colnames(x)=="1"] <- "(Intercept)"
  
  # General parameters
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df = n-p-1 # degrees of freedom
  
  # coefficient estimates = Beta 
  xTx = t(x) %*% x
  xTy = t(x) %*% y
  beta = solve(xTx) %*% xTy
  beta
  
  # standard errors
  yhat = x %*% beta
  e = y - yhat
  eTe = sum(e^2) 
  sigma2 = eTe/df  ## standard error = (sum(e^2)/n-p-1), df = n-p-1
  
  varcov = sigma2 * solve(xTx) # variance-covariance matrix
  serrors = sqrt(diag(varcov)) # standard errors of each parameter
  
  # t-statistics and p-value
  tstats = beta/serrors
  pval = 2*pt(abs(tstats), df, lower.tail=FALSE)
  
  # 95% confidence intervals
  upr = beta + qt(.975, df)*serrors
  lwr = beta - qt(.975, df)*serrors
  
  # R squared
  Rsq = sum((yhat-mean(y))^2)/sum((y-mean(y))^2)
  
  # F statistic
  MSreg = sum(yhat^2)/p
  MSres = sum(e^2)/(n-p-1)
  Fstat = round(MSreg/MSres, 3)
  Fstat_pval = round(pf(Fstat, p, df, lower.tail=FALSE), 3)
  Fstat2 = paste("F-statistic:", Fstat, "on", p, "and", df, "DF, p-value:", Fstat_pval)
  
  # creating coefficient matrix
  coefficients = data.frame(cbind(beta, serrors, tstats, pval, lwr, upr))
  names(coefficients)=c('Estimate', 'Std. Error', 'T-Statistic', 'P-Value', 'Lower 95% CI', 'Upper 95% CI')
  
  list(coefficients=coefficients, varcov=varcov, Rsq=Rsq, Fstat=Fstat2)
}

###############

# # Applying to the actual data

# # load data
# load("/Users/moohyungcho/Desktop/Spring 2015/MLE/Lab/midTermData.rda")

# # Run the various models
# form = formula(gini_net_std ~ ELF_ethnic + polity2)

# # 1. Original model
# model = ols(formula=form, data=data)
# model

# # 2. Model with listwise deletion
# modelListDel = ols(formula=form, data=dataMiss)
# modelListDel

# # 3. Model with imputation (Amelia)
# modelAmelia = ols(formula=form, data=dataMiss, impute=TRUE)
# modelAmelia

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

expect_that(olsSM(form, data)$coefficients, equals(ols(form, data=data)$coefficients))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))