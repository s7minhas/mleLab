rm(list=ls())
ols = function(formula, data, impute=FALSE){
  #first impute data if missing is true
  if(impute==TRUE) {
    #print(data)
    library(Amelia)
    set.seed(6886)
    imputed = amelia(x=data, m=1)
    data = imputed$imp$imp1
    #print(data)
  }
  if(impute==FALSE) {
    data = na.omit(data)
    
  }
  #print(data)
  library(stats)
  
  #print(data)
  # Retrieve vars from formula input
  dv = all.vars(formula)[1]
  #print(dv)
  ivs = all.vars(formula)[ 2:length(all.vars(formula)) ]
  #print(ivs)
  
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  # General parameters
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df = n-p-1 # degrees of freedom
  #print(n)
  #print(p)
  
  #calculating beta
  xTx = t(x) %*%x
  xTy = t(x) %*%y
  beta= solve(xTx) %*% xTy
  rownames(beta)[1] = "(Intercept)"
  #print(beta)
  yhat = x %*% beta #yhats
  resid = y - yhat
  #print(resid)
  sigmasq = (sum(resid^2) )/df
  varcov = sigmasq * solve(xTx) #produces vcov matrix
  rownames(varcov)[1] = "(Intercept)"
  colnames(varcov)[1] = "(Intercept)"
  #print(varcov)
  stderror = sqrt(diag(varcov))
  #print(stderror)
  
  #calculating conf intervals
  upper = beta + qt(.975,df)*stderror
  lower = beta - qt(.975,df)*stderror
  #print(upper)
  
  #calculating t-statistic
  tstats = beta/stderror
  
  pval = 2*pt(abs(tstats), df=df, lower.tail=FALSE)
  #print(tstats)
  #print(pval)
  
  #calculating R squared (3rd object)
  Rsq = 1 - ( sum(resid^2)/(sum( (y-mean(y) )^2) ) )
  #print(Rsq)
  
  #calculating F statistic (4th object)
  meansqreg = ( sum(yhat^2))/p
  meansqresid = sum(resid^2)/df
  Fstatistic = meansqreg/meansqresid
  Fstatistic = round(Fstatistic, digits = 3)
  pvalF = pf(Fstatistic, p, df, lower.tail=FALSE)
  pvalF = round(pvalF, digits = 3)
  Fstat = paste0('F-statistic: ', Fstatistic, ' on ', p, ' and ', df, ' DF, p-value: ', pvalF)
  #print(Fstat)
  
  #compiling eveyrthing into desired output
  coefficients = cbind(beta, stderror, tstats, pval, lower, upper)
  colnames(coefficients) = c('Estimate', 'Std. Error', 'T-Statistic', 'P-Value', 'Lower 95% CI', 'Upper 95% CI' )
  rownames(coefficients)[1] = "(Intercept)"
  
  #print(coefficients)
  output = list(coefficients, varcov, Rsq, Fstat)
  names(output)[1] = "coefficients"
  names(output)[2] = "varcov"
  names(output)[3] = "Rsq"
  names(output)[4] = "Fstat"
  return(output)
  
}

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))