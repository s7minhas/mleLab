#############################################################################
# INPUT: 
# 1) formula = formula object specifying model to be run
# 2) data    = dataframe object containing data (DV and IVs)
# 3) impute  = logical parameter indicating whether missing data should be
#              imputed via Amelia (m=1)
#      ^ default  = FALSE
#############################################################################
# OUTPUT: list object
# 1) coefficients = (p+1)x6 matrix object summarizing coefficient estimates
#      ^ rownames = '(Intercept)' and original IV names
#      ^ colnames = 'Estimate', 'Std. Error', 'T-Statistic', 'P-Value', 
#                   'Lower 95% CI', 'Upper 95% CI'
# 2) varcov       = matrix object containing variance-covariance matrix of 
#                   coefficient estimates
#      ^ rownames = rownames(coefficients)
#      ^ colnames = rownames(coefficients)
# 3) Rsq          = numeric object containing value of model's R^2 statistic
# 4) Fstat        = character object providing F-statistic, relevant degrees 
#                   of freedom, and p-value (rounded to 3 decimal places) 
#      ^ format   = 'F-statistic: 0.204 on 1 and 8 DF,  p-value: 0.663'
#############################################################################

ols = function(formula, data, impute=FALSE) {
  dv     = all.vars(formula)[1]
  ivs    = all.vars(formula)[2:length(all.vars(formula))]
  if(impute==T) {
    if(! 'Amelia' %in% installed.packages()[,1]) { 
      install.packages('Amelia', repos='http://cran.rstudio.com/') 
    }
    suppressMessages(library('Amelia', character.only=T))
    # set.seed(6886)
    set.seed(6886)
    found = amelia(data, m=1, p2s=0)
    Adata = found$imputations$imp1
    y     = Adata[,dv]
    x     = data.matrix(cbind(1, Adata[,ivs]))
  } else {
    data  = na.omit(data[,c(dv, ivs)])
    y     = data[,dv]
    x     = data.matrix(cbind(1, data[,ivs]))
  }
  n      = nrow(x)
  p      = length(ivs)
  df     = n - p - 1
  xTx    = t(x) %*% x
  xTy    = t(x) %*% y
  betas  = solve(xTx) %*% xTy
  yHats  = x %*% betas
  resids = y - yHats
  msRes  = sum(resids^2) / df
  msReg  = sum(yHats^2) / p
  fstat  = round((msReg / msRes), 3)
  fpval  = round(pf(fstat, p, df, lower.tail=F), 3)
  Fstat  = paste('F-statistic:', fstat, 'on', p, 'and', df, 'DF, p-value:', fpval)
  ssRes  = sum(t(resids) %*% resids)
  ssReg  = sum((mean(y) - yHats)^2)
  ssTot  = ssRes + ssReg
  Rsq    = ssReg / ssTot
  sigma2 = ssRes / df
  varcov = solve(xTx) * sigma2
  ses    = sqrt(diag(varcov))
  tstats = betas / ses
  pvals  = 2 * pt(abs(tstats), df, lower.tail=F)
  band   = qt(.975, df) * ses
  lower  = betas - band
  upper  = betas + band
  coefficients              = cbind(betas, ses, tstats, pvals, lower, upper)
  colnames(coefficients)    = c('Estimate', 'Std. Error', 'T-Statistic', 
                                'P-Value', 'Lower 95% CI', 'Upper 95% CI')
  rownames(coefficients)[1] = '(Intercept)'
  rownames(varcov) = colnames(varcov) = rownames(coefficients)
  output = list(coefficients, varcov, Rsq, Fstat)
  names(output) = c('coefficients', 'varcov', 'Rsq', 'Fstat')
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