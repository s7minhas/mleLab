rm(list=ls())
# the function
ols = function(formula, data, impute=FALSE){
  # load var names
  dv = all.vars(formula)[1]
  ivs = all.vars(formula)[ 2:length(all.vars(form)) ]
  # impute missing if impute==TRUE
  if(impute==TRUE){
    data = amelia(x=data[,c(dv, ivs)], m=1)
    data = data$imp$imp1
  }
  # listwise delete missing if impute==FALSE  
  if(impute==FALSE){
    data = na.omit(data)
  }
  # obtain data and general parameters from inputs
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  n = nrow(x)
  k = length(ivs)
  df = n-k-1
  # estimate model parameters
  # coefficients
  xTx = t(x) %*% x
  xTy = t(x) %*% y
  beta = solve(xTx) %*% xTy
  # standard errors
  yHat = x %*% beta
  sigma2 = sum((y - (yHat))^2)/df
  varcov = (sigma2 * solve(xTx))
  serrors = sqrt(diag(varcov))
  # confidence intervals  
  upper95 = beta + (qt(.975, df)*serrors)
  lower95 = beta - (qt(.975, df)*serrors)
  # t-statistics and p-values
  tstats = beta/serrors
  pvals = 1.96*pt(abs(tstats), df, lower.tail=FALSE)
  # model performance metrics
  R2 = 1 - sum((yHat - mean(y))^2)/sum((y - mean(y))^2)
  MSreg = sum(yHat^2)/(k-1)
  MSres = sum((y - (yHat))^2)/(n-k)
  Fstat = (MSreg)/(MSres)
  FF = pf(Fstat, k-1, n-k, lower.tail=FALSE)
  # create output
  # coefficients
  coefficients = cbind(beta, serrors, tstats, pvals, lower95, upper95)
  coefficients = as.matrix(coefficients)
  rownames(coefficients)[1] = "(Intercept)"
  colnames(coefficients) = c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
  # varcov
  varcov = as.matrix(varcov)
  rownames(varcov)[1] = "(Intercept)"
  colnames(varcov)[1] = "(Intercept)"
  # Rsq
  Rsq = as.numeric(R2)
  # Fstat
  Fstatistic = as.character(round(Fstat, digits = 3))
  Fpv = as.character(round(FF, digits = 3))
  df1 = as.character(k-1)
  df2 = as.character(n-k)
  Fstat = paste0(c("F-statistic: ", Fstatistic, " on ", df1, " and ",
                   df2, " DF, p-value: ", Fpv), collapse="")
  # Keep them
  output = list(coefficients, varcov, Rsq, Fstat)
  names(output) = c("coefficients", "varcov", "Rsq", "Fstat")
  return(output)
}

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

ols(form, data=data)
olsSM(form, data=data)


expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))