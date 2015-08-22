# ##---------------------------------##
# ## MLE Midterm: OLS function       ##
# ## Jochen Rehmert                  ##
# ## 07.03.2015                      ##
# ##---------------------------------##


rm(list=ls())

# load("C:/Users/jochen/Dropbox/4_Fachsemester/MLE/midterm/midTermData.rda")

# form = formula(gini_net_std ~ polity2 + ELF_ethnic)


# #----- problems 
# # 1) output for p-val of Fstat does not round as wanted




# set.seed(6886)

ols <- function(formula, data, impute=FALSE){
  #-------------------------------------------- handling missing data  
  if((!require(Amelia)) && (impute==TRUE)){print("The 'Amelia' package is not installed. Please use the 'install.packages()' command to install 'Amelia'.")}
  if(impute == TRUE){      
  set.seed(6886)       
    data <-  Amelia::amelia(data, m = 1)[['imputations']][['imp1']]
    data <- na.omit(data)       
  } else {
    data <- na.omit(data)            
  } 
  mf <- model.frame(formula=formula, data=data)
  #-------------------------------------------- model matrix % parameter     
  X  <- model.matrix(attr(mf, "terms"), data=mf)
  y  <- model.response(mf)
  n  <- nrow(mf)
  p  <- length(mf)
  df <- n-p
  #-------------------------------------------- estimation
  # beta: (X'X)^-1 %*% X'Y
  xTx = t(X) %*% X # X'X   dimensions: (p+1) * (p+1)
  xTy = t(X) %*% y # X'Y
  beta = solve(xTx) %*% xTy # getting beta
  # standard errors: sqrt(sigma^2 * (X'X)^-1)
  yhat = X %*% beta # predicted values
  e = y - yhat # residuals
  sigma2 = (sum(e^2))/(n-p) # e'e = sum of squared residuals
  varcov = sigma2 * solve(xTx) # variance-covariance matrix
  serror = sqrt(diag(varcov)) # standard error 
  # t-stats 
  tstat = beta/serror
  # p-value
  pval = 2 * (pt(abs(tstat), df, lower.tail=FALSE)) 
  #------------------------------------------- diagnostics & fit
  # R-squared 
  SSreg = sum((yhat - mean(y))^2) # regression sum of squares
  SStot = sum((y - mean(y))^2) # total sum of squares
  rsq = SSreg/SStot 
  # Adjusted R-squared
  adjr2 = 1-(1-rsq)*((n-1)/df)
  # F statistic 
  MSReg = sum(yhat^2)/(p-1)
  MSRes = sum(e^2)/df
  fstat = MSReg/MSRes # f statistic
  fpval = pf(fstat, p, df, lower.tail=FALSE) # p-value for f statistic
  # root mean squared error
  rmse = sqrt(mean(e^2))
  # confidence intervals
  cilow95 = beta - qt(.975,df) * serror
  cihigh95 = beta + qt(.975,df) * serror
  #------------------------------------------- output essentials & cosmetics
  coefMat <- cbind(
    'Coefficient' = beta,
    'Std. Error' = serror,
    'T Stat' = tstat,
    'P-Value' = pval,
    'Low CI' = cilow95,
    'High CI' = cihigh95
  )
  colnames(coefMat) <- c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
  
  ftest = paste0("F-statistic: ", round(fstat,3)," on ", (p-1), " and " ,df, " DF, p-value: ", round(fpval,3))
  
  return(list(
    coefficients= coefMat, 
    varcov=varcov,
    Rsq=rsq,
    Fstat=ftest
      ))
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

# model = ols(formula=form, data=data)

# modelListDel = ols(formula = form, data=dataMiss)

# modelAmelia = ols(formula = form, data=dataMiss, impute=TRUE)



