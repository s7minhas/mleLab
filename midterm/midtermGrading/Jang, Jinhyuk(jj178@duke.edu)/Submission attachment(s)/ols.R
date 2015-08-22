# Lab 7
# MLE Midterm: Jinhyuk Jang
# March 17, 2015

###################################################
# Start with a clean workspace
rm(list=ls())

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

# # Set a theme for gg
# theme_set(theme_bw())

# # Functions that I use frequently
# char = function(x){ as.character(x) }
# num = function(x){ as.numeric(char(x)) }

# # Relevant paths
# labPath=file.path('C:','Users', 'Jinhyuk Jang', 'Documents', 'MLE Lab')
# ###################################################

# ###################################################
# # Set path
# setwd(labPath)

# # Load Data
# load("midTermData.rda")

###################################################

###################################################

ols <- function(formula, data, impute=FALSE){
  # Create a condition for imputed data
  if(impute==TRUE){
    # Create imputed dataset
    dataAmelia = amelia(x=data, m=1) # use m=1 
    names(dataAmelia$imp)
    data=data=dataAmelia$imp$imp1
    # Retrieve vars from formula input
    dv = all.vars(formula)[1]
    ivs = all.vars(formula)[ 2:length(all.vars(formula)) ]
    # Create matrix with column for intercept and
    ## data from independent variables
    y = data[,dv]
    x = data.matrix(cbind(1, data[,ivs]))
    # General parameters
    n = nrow(x) # Number of observations
    p = length(ivs) # Number of parameters
    df = n-p-1 # degrees of freedom
    # Calculate coefficient estimates
    xTx = t(x) %*% x
    xTy = t(x) %*% y
    beta = solve(xTx) %*% xTy
    # Calculate standard errors
    yhat = x %*% beta
    e = y-yhat
    ssr=sum(e^2)
    df=n-p-1
    sigma2=ssr/df
    varcov = sigma2 * solve(xTx)
    serrors = sqrt(diag(varcov))
    # Calculate tstats
    tstats = beta/serrors
    # Calculate pvalue
    pvalue = 2*pt(-abs(tstats), df)
    # Calculate 95% confidence interval
    upper = beta + 1.96*serrors
    lower = beta - 1.96*serrors
    # Calculate R2
    ssReg = sum((yhat-mean(y))^2)
    ssTot = sum((y-mean(y))^2)
    R2 = ssReg/ssTot
    # F statistic
    MsReg = sum(yhat^2)/p
    MsRes = sum(e^2)/df
    Fstat = MsReg/MsRes
    Fstat = round(Fstat, digits = 3)
    # Degree of freedom
    P = p
    Df = df
    # Pvalue in F states
    Pvalue = pf(Fstat, p, df, lower.tail=FALSE)
    Pvalue = round(Pvalue, digits = 3)
    coefficients=matrix(cbind(beta, serrors, tstats, pvalue, lower, upper), nrow=3, ncol=6)
    rownames(coefficients) <- c('(Intercept)','ELF_ethnic','polity2')
    colnames(coefficients) <- c('Estimate','Std. Error','T-Statistic','P-Value', 'Lower 95% CI', 'Upper 95% CI')
    varcov=matrix(varcov, nrow=3, ncol=3)
    rownames(varcov) <- c('(Intercept)','ELF_ethnic','polity2')
    colnames(varcov) <- c('(Intercept)','ELF_ethnic','polity2')
    Rsq = as.numeric(R2)
    Fstat=paste("F-statistic:",Fstat,"on", P,"and", Df,"DF,","p-value:", Pvalue)
    list(coefficients=coefficients, varcov=varcov, Rsq=Rsq, Fstat=Fstat)
  } else {
    # Omit Missing Values
    data=na.omit(data)
    # Retrieve vars from formula input
    dv = all.vars(formula)[1]
    ivs = all.vars(formula)[ 2:length(all.vars(formula)) ]
    # Create matrix with column for intercept and
    ## data from independent variables
    y = data[,dv]
    x = data.matrix(cbind(1, data[,ivs]))
    # General parameters
    n = nrow(x) # Number of observations
    p = length(ivs) # Number of parameters
    df = n-p-1 # degrees of freedom
    # Calculate coefficient estimates
    xTx = t(x) %*% x
    xTy = t(x) %*% y
    beta = solve(xTx) %*% xTy
    # Calculate standard errors
    yhat = x %*% beta
    e = y-yhat
    ssr=sum(e^2)
    df=n-p-1
    sigma2=ssr/df
    varcov = sigma2 * solve(xTx)
    serrors = sqrt(diag(varcov))
    # Calculate tstats
    tstats = beta/serrors
    # Calculate pvalue
    pvalue = 2*pt(-abs(tstats), df)
    # Calculate 95% confidence interval
    upper = beta + qt(.975, df)*serrors
    lower = beta - qt(.975, df)*serrors
    # Calculate R2
    ssReg = sum((yhat-mean(y))^2)
    ssTot = sum((y-mean(y))^2)
    R2 = ssReg/ssTot
    # F statistic
    MsReg = sum(yhat^2)/p
    MsRes = sum(e^2)/df
    Fstat = MsReg/MsRes
    Fstat = round(Fstat, digits = 3)
    # Degree of freedom
    P = p
    Df = df
    # Pvalue in F states
    Pvalue = pf(Fstat, p, df, lower.tail=FALSE)
    Pvalue = round(Pvalue, digits = 3)
    coefficients=matrix(cbind(beta, serrors, tstats, pvalue, lower, upper), nrow=3, ncol=6)
    rownames(coefficients) <- c('(Intercept)','ELF_ethnic','polity2')
    colnames(coefficients) <- c('Estimate','Std. Error','T-Statistic','P-Value', 'Lower 95% CI', 'Upper 95% CI')
    varcov=matrix(varcov, nrow=3, ncol=3)
    rownames(varcov) <- c('(Intercept)','ELF_ethnic','polity2')
    colnames(varcov) <- c('(Intercept)','ELF_ethnic','polity2')
    Rsq = as.numeric(R2)
    Fstat=paste("F-statistic:",Fstat,"on", P,"and", Df,"DF,","p-value:", Pvalue)
    list(coefficients=coefficients, varcov=varcov, Rsq=Rsq, Fstat=Fstat)
  }
}

###################################################

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

rownames(olsSM(form, data)[[1]])[2]==rownames(ols(form, data=data)[[1]])[2]

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))

# ###################################################

# # First set a seed
# set.seed(6886)

# # Define formula
# form = formula(gini_net_std ~ ELF_ethnic + polity2)

# # Run the various models
# model = ols(formula=form, data=data)
# modelListDel = ols(formula=form, data=dataMiss)
# modelAmelia = ols(formula=form, data=dataMiss, impute=TRUE)
