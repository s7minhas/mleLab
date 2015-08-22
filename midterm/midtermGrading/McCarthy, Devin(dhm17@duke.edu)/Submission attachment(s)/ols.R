rm(list=ls()) 
# setwd("c:/Users/devin_000/Documents/R/MLE")

# loadPkg=function(toLoad){
#   for(lib in toLoad){
#     if(! lib %in% installed.packages()[,1])
#     { install.packages(lib, repos='http://cran.rstudio.com/') }
#     suppressMessages( library(lib, character.only=TRUE) ) }
# }

# packs=c('foreign', 'lmtest', 'sandwich', 'Amelia', 'sbgcop')
# loadPkg(packs)

# midterm = load("midTermData.rda")
# set.seed(6886)

# datamisslistdel <- subset(dataMiss, !is.na(dataMiss$polity2))
# form = formula(gini_net_std ~ ELF_ethnic + polity2)


ols <- function(formula, data, impute=FALSE ){
  # form = formula(gini_net_std ~ ELF_ethnic + polity2)
  
  dv = all.vars(form)[1]
  ivs = all.vars(form)[ 2:length(all.vars(form)) ]
  
  if(impute==TRUE)  {
    set.seed(6886)
    dataimp <- amelia(x=data, m=1)
    data <- dataimp$imputations[[1]]
                   }
  
  data = na.omit(data)
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  
  n = nrow(x) 
  p = length(ivs)
  df = n-p-1
  
  
  xTx = t(x) %*% x
  
  xTy = t(x) %*% y
  
  beta = solve(xTx) %*% xTy
  
  yhat = x %*% beta
  e = y-yhat
  
  sigma2 = sum(e^2)/df 
  varcov= sigma2 * solve(xTx)
  serrors = sqrt (diag(varcov))
  tstats = beta/serrors
  pvals = 2*pt(abs(tstats), df=df, lower.tail=F)
  ssReg=sum((yhat-mean(y))^2)
  ssTot=sum((y-mean(y))^2) 
  
  Rsq = ssReg/ssTot
  adjr2 = 1- (1-Rsq)*((n-1)/df)
  
  msReg = sum(yhat^2)/p
  msRes = sum(e^2)/df
  
  Fstat = msReg/msRes
  
  Fpdf = pf(Fstat, p, df, lower.tail=F)
  
  beta.hat = cbind(beta, serrors, tstats, pvals)
  
  uppers = beta + serrors*qt(.975,df)
  lowers = beta - serrors*qt(.975,df)
  
  coefficients = cbind(beta, serrors, tstats, pvals, lowers, uppers)
  colnames(coefficients) = c( "Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
  rownames(coefficients)[1] = c("(Intercept)")
  colnames(varcov)[1] = c("(Intercept)")
  rownames(varcov)[1] = c("(Intercept)")
  
  Fstat <- paste0("F-statistic: ", round(Fstat, 3), " on ", p, " and ", df, " DF, p-value: ", round(Fpdf, 3))
  # Fstat <- noquote(paste(Fstat, collapse=""))
  
  list(coefficients=coefficients, varcov=varcov, Rsq=Rsq, Fstat=Fstat)
}

library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

# olsSM(form, data)
# ols(form, data=data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))

expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))

set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))

# model = ols(formula=form, data)
# model

# modelListDel =ols(formula=form, na.omit(dataMiss))
# modelListDel

# modelAmelia = ols(formula=form, dataMiss, impute=TRUE)
# modelAmelia