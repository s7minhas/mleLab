#This function has three inputs: 
#1. formula = A formula class object that specifies the model to be run (e.g., y ~ x1 + x2)
#2. data = A dataframe object that contains the data for the DV and covariates
#3. impute = A logical parameter (is TRUE or FALSE) to indicate whether missing data should be imputed using Amelia (set number of imputed datasets to 1). By default this parameter should be set to FALSE.

#The output of this function should be a list that contains four named objects:
#1. coefficients = A matrix object that contains rows for every coefficient estimate (including the Intercept). The row with the intercept should be labeled as "(Intercept)", the remaining rows should simply take whatever the name is of the independent variables. This matrix should contain six columns with the following names: "Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI".
#2. varcov = A matrix object that contains the variance covariance matrix of the estimated coefficients. The rows and column names should match the rownames of the coefficients matrix.
#3. Rsq =  A numeric object that is just the value of the R2 statistic for the model.
#4. Fstat = A character object that gives the F statistic, relevant degrees of freedom, and p-value. It should be in the form: "F-statistic: 0.204 on 1 and 8 DF,  p-value: 0.663". The F-statistic and p-value should be rounded to three decimal places. 

#############################################################

ols <- function(formula, data, impute=FALSE){
  dv = all.vars(formula)[1] #retrieve dependent vars from formula
  ivs = all.vars(formula)[ 2:length(all.vars(formula))] #retrieve independent vars from formula
 
#impute missing if impute=TRUE, listwise deletion if impute=FALSE
if (!impute) {  
  fulldata = na.omit(data)
  y = fulldata[,dv] #retrieve vector of dependent var from data
  x = data.matrix(cbind(1, fulldata[,ivs]))#retrieve matrix of independent var from data
} else { 
  library(Amelia)
  dataAmelia = amelia(x=data, m=2)
  imputedata = dataAmelia$imp$imp1
  y = imputedata[,dv] #retrieve vector of dependent var from data
  x = data.matrix(cbind(1, imputedata[,ivs]))#retrieve matrix of independent var from data
}

  n = nrow(x) #number of observations
  p = length(ivs) #number of parameters
  df = n-p-1 #degrees of freedom
  xTx = t(x) %*% x #calculating X'X
  xTy = t(x) %*% y #calculating X'Y
  beta = solve(xTx) %*% xTy #calculating beta
  yhat = x %*% beta #estimating dependent vars
  error = y-yhat #error term
  sigma2=sum(error^2)/df #sum of squared residuals
  varcov = sigma2 * solve(xTx) #variance-covariance matric
  rownames(varcov) <- c('(Intercept)', ivs) #change row names
  colnames(varcov) <- c('(Intercept)', ivs) #change column names
  std.error = as.matrix(sqrt(diag(varcov))) #calculate standard errors
  tstat = beta/std.error #calculate t-statistics
  pval = pt(abs(tstat), df=df, lower.tail=FALSE)*2 #calculate p-values
  lowerq = beta + qt(0.025, df=df)*std.error/sqrt(n) #calculate lower confidence interval
  upperq = beta + qt(0.975, df=df)*std.error/sqrt(n) #calculate upper confidence interval
  coefs = cbind(beta, std.error, tstat, pval, lowerq, upperq) #combine into matrix object
  rownames(coefs) <- c('(Intercept)', ivs) #change row names
  colnames(coefs) <- c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI") #change column names
  ssreg = sum((yhat-mean(y))^2)
  sstotal = sum((y-mean(y))^2)
  r2 = ssreg/sstotal #calculate R-squared
  msreg = sum(yhat^2)/p
  msres = sum(error^2)/df
  fstat = round(msreg/msres,3) #calculate F-statistic 
  probf = round(pf(fstat,p,df,lower.tail=FALSE),3) #calculate p-value
  sentence <- paste0("F-statistic: ",fstat, " on 1 and ",df, " DF,  p-value: ", probf, ".") #create character object
  
  output <- list(coefficients=round(coefs,4), varcov=round(varcov,6), Rsq=r2, Fstat=sentence)
  print(output)
}


library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

olsSM(form, data)
ols(form, data=data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))