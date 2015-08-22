#Jeremy Spater OLS function
#POLSCI 733 MLE
#3/10/15
rm(list=ls())
ols <- function(formula, data, impute=FALSE){
  # Retrieve vars from formula input
  
  if(impute==TRUE){
    print("Imputing missing data with Amelia")
    AmeliaOutput=amelia(x=data,m=1)
    data=AmeliaOutput$imp$imp1
  }
  
  #listwise deletion of missing data (won't be any if impute==true)
  data=na.omit(data)
  
  formula=form
  dv = all.vars(formula)[1]
  ivs = all.vars(formula)[ 2:length(all.vars(formula)) ]
  
  # Create matrix with column for intercept and
  ## data from independent variables
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))#include 1
  
  # General parameters
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df = n-p-1 # degrees of freedom
  
  # Derive OLS estimator of beta
  xTx = t(x) %*% x 
  xTxinv=solve(xTx)
  xTy = t(x) %*% y
  beta = solve(xTx) %*% xTy
  yhat=x %*% beta
  
  #calculate residuals
  e=y-yhat
  ssr = t(e) %*% e
  ssr2=sum(e^2)
  
  #calculate standard error
  stder=ssr/df
  
  #calculate estimated covariance matrix of betahat
  Sigma2=as.numeric(stder)*xTxinv #square root of these diagonals give std errors as reported in lm1; confirmed p620 of S&W
  sigma=sqrt(diag(Sigma2))
  
  #t values
  t=beta/sigma
  
  #p values
  pvals=2*pt(abs(t),df,lower.tail=FALSE)#2 accounts for prob of getting more-extreme value, on other side of distribution
  
  #confidence intervals
  numstd=qt(0.975,43,lower.tail=TRUE)
  conf=matrix(data=NA,nrow=p+1,ncol=2)
  conf=cbind(beta-numstd*sigma,beta+numstd*sigma)
  
  #make matrix
  coefficients=cbind(beta,sigma,t,pvals,conf)
  colnames(coefficients)=c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
  
  #rename varcov matrix
  varcov=Sigma2
  
  #calculate r2
  ymean=sum(y)/n
  tss=sum((y-ymean)^2)
  ess=sum((yhat-ymean)^2)
  rss=ssr
  Rsq=1-rss/tss
  
  #F statistic
  
  MSreg = sum(yhat^2)/p
  MSres = sum(e^2)/(n-p-1)
  fstat=MSreg/MSres
  pvalf=1-pf(fstat,p,n-p-1)
  
  Fstat =paste("F-statistic: ", as.character(round(fstat,3)), " on ", as.character(p), " and ", as.character(n-p-1), " DF,  p-value: ", format(pvalf,digits=4,scientific=TRUE))
  
  output=list(coefficients, varcov, Rsq, Fstat)
  names(output)=c("coefficients", "varcov", "Rsq", "Fstat")
  
  return(output)
  
}

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

olsSM(form, data=data)
ols(form, data=data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))