ols <- function(formula, data, impute=FALSE){
  
  #Step 1: do an IF function to run Amelia
  if (impute==TRUE){
    a.out <- amelia(x=data,m=1)
    data <-a.out$imputation[[1]]
    #data <- amelia(x=data, m=1)
  }
  
  #Step 2: tells the function to omit missing data... if ran amelia, nothing will be omitted
  data <-na.omit(data)
  
  #Step 3: create value (DV & IV) that come from the elements of the TBD formula
  dv = all.vars(formula)[1] #for ease of coding, this labels the first variable in formula as the DV
  ivs = all.vars(formula)[ 2:length(all.vars(formula)) ] #same... but label the others vars as IVs
  
  #Step 4: create matrix with columns for intercept and IDVs: (n=# observation; p=# parameters)
  y = data[,dv]  #this creates a n*1 vector of DV values by subsetting data matrix 
  x = data.matrix(cbind(1, data[,ivs])) #creates n*(p+1) Design Matrix of intercept and coefficients.  
  
  colnames(x)[1] <- c('(Intercept)') #this labels the 1st column of Design Matrix as Intercept 
  
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df = n-p-1 # degrees of freedom is # observatioms - # coefficients -1
  
  #Step 5: determine the estimated Coefficients from OLS "Beta Hat"
  xtx <-  t(x) %*%x #matrix multiply transposed Design Matrix against itself(sizde of matrix?)
  xty <- t(x) %*% y #matrix multiply transposed Design Matrix againt vector of DVs 
  betahat <-  solve(xtx) %*% xty #p+1 vector of OLS estimated coefficients from xTx and xTy
  
  #Step 6: determine Standard errors
  yhat <- x %*% betahat #n*1 vector of Predicted values of y (Design matrix * Coeff vector)
  resid <- y- yhat #n*1 vector of Residuals or Errors (difference between real and estimated)
  error2 <- sum(resid^2) #numeric (1*1 vector) Residual Sum of Squares
  sigmasq <- error2/df #numeric (1*1 vector) of RSS/degress of freedom
  varcov <- sigmasq*solve(xtx) #p+1 x p+1 variance/covariance matrix (xTX * sigma2)
  stderrors <-  sqrt(diag(varcov)) #p+1 vector of Standard Deviations (diagonal of varcov matrix
  
  #Step 7: determine the T-stat
  tstats <- betahat/stderrors #p+1 vector that gives a T-stat for each parameter and constant
  
  #Step 8: determine the p-value (probablity you would find this answer assuming Null Hypoth true)
  pval <-  2 * pt(abs(tstats), df=n-p-1, lower.tail=F)
  pval <- round(pval, digits=3)
  
  #Step 9: determine confidence intervals    
  upper95 = betahat + qnorm(.975)*stderrors #n*1 vector of upper 95% CI for each y-hat value
  lower95 = betahat - qnorm(.975)*stderrors #n*1 vector of lower 95% CI for each y-hat value 
  
  #Step 10: creates the matrix object for output
  coefficients <- (cbind(betahat, stderrors,
                         round(tstats, digits =3), round(pval,digits = 3),
                         lower95, upper95)) #this binds all of the column names into a single value
  colnames(coefficients) <- c("Estimate", "Std. Error", "T-Statistic", "P-Value","Lower 95% CI", "Upper 95% CI") #this labels the columns 
  
  #Step 11: creates the Variance-Covariance Matrix for the OLS
  varcov <- sigmasq*solve(xtx) #p+1 x p+1 variance/covariance matrix (xTX * sigma2)
  
  #Step 12: create the R^2
  SSreg <-  sum((yhat-mean(y))^2) #sum of the square of difference btw yhat and avg y (numeric)
  SStot <- sum((y-mean(y))^2) #sum of the square of difference btw actual y and avg y (numeric)
  Rsq <-  SSreg/SStot #R-squared (numeric)
  names(Rsq) <- "R-squared" #labels the value (char)
  
  #Step 13: calculate the F-statistic..output is a sentence statement with variables in it
  MSreg <-  sum(yhat)^2/ p
  MSresid <-  sum(resid^2)/df
  Fstat <-  round(MSreg/ MSresid, digits =3)
  pval2 <- pf(Fstat, p, df, lower.tail=FALSE) #p=#parameters, df+degree freedom n-p-1
  FSoutput <- paste0("F-statistic: ", round(Fstat,digits= 3)," on ",p," and ",df," DF,"," p-value: ",round(pval2,digits=3))
  
  #Step 14: putting together the output
  outList <- list(coefficients, varcov, Rsq, FSoutput) #Steps:Coeff=5-10; varcov=11, Rsq=12, Fstat=13
  names(outList) <- c('coefficients', 'varcov', 'Rsq', "Fstat")
  return(outList)
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