###################################################
# MIDTERM ASSIGNMENT
# Soomin Oh
# March 17, 2015

####################### setting up the environment ###################
# Set up workspace
rm(list=ls()) 
# setwd("~/Dropbox/Grad School/Spring 2015/MLE/lab/lab7")

# # Function to load packages
# loadPkg=function(toLoad){
#   for(lib in toLoad){
#     if(! lib %in% installed.packages()[,1])
#     { install.packages(lib, repos='http://cran.rstudio.com/') }
#     suppressMessages( library(lib, character.only=TRUE) ) }
# }

# # Load libraries
# packs=c('foreign', 'lmtest', 'sandwich', 'Amelia', 'sbgcop', 'ggplot2','arm')
# loadPkg(packs)

# midterm <- load("midTermData.rda")

# ################## OLS FUNCTION FROM SCRATCH ##########################

# set.seed(6886)
# form = formula(gini_net_std ~ ELF_ethnic + polity2)

ols <- function(formula, data, impute=FALSE){
  # for listwise deletion and imputation with Amelia
  if (impute == F){
    data <- data[which(!is.na(apply(data, 1, function(x) sum(x)))), ]
  }
  else {
    coln <- colnames(data)
    set.seed(6886)
    mtAmelia <- amelia(x=data, m=1)
    data <- as.data.frame(mtAmelia$imputations)
    colnames(data) <- coln
    rm(coln)
  }
  

  dv <- all.vars(form)[1]
  ivs <- all.vars(form)[2:length(all.vars(form))]

  y <- data[,dv]
  x <- data.matrix(cbind(1,data[,ivs]))

  n = nrow(x)
  p = length(ivs)
  df = n-p-1

## Coefficient estimates
# X'X
  xTx = t(x) %*% x

# X'Y
  xTy = t(x) %*% y

#calculating beta
  beta = solve(xTx) %*% xTy

# Calculating standard errors
  yhat = x %*% beta
  e = y - yhat
  
  sst <- sum((y-mean(y))^2)
  sigma2 = sum(e^2)
  adj = sum(e^2)/(n-p-1)

  varcov = adj*solve(xTx)
  serrors = sqrt(diag(varcov))

# t-values
  tstats = beta/serrors

# p-values
  pval = 2*pt(abs(tstats),df=n-p-1,lower.tail=F)

# lower 95% CI
error <- qnorm(0.975)*sigma2/sqrt(n)
lower <- beta - error

# upper 95% Ci
upper <- beta+error

# R-squared
 Rsq <- 1-(sigma2/sst)

# F-Stat
fstat = ((sum(yhat^2)/p)/(sum(e^2)/df))
fstatp = pf(fstat,p,df,lower.tail=FALSE)

coefficients <- cbind(beta, serrors, tstats, pval, lower, upper)
colnames(coefficients) <- c('Estimate', 'Std. Error', 'T-Statistic', 'P-Value', 'lower 95% CI', 'Upper 95% CI')
rownames(coefficients) <- c('(Intercept)', 'ELF_ethnic', 'polity2')

colnames(varcov) <- c('(Intercept)', 'ELF_ethnic', 'polity2')
rownames(varcov) <- c('(Intercept)', 'ELF_ethnic', 'polity2')

fstat3dp <- round(fstat, digits = 3)
fstatp3dp <- round(fstatp, digits = 6)

fstatistic <- sprintf("F-statistic: %s on %s DF, p-value: %s", fstat3dp, df, fstatp3dp)

print(list(coefficients=coefficients, varcov=varcov, Rsq=Rsq, fstatistic=fstatistic))
}


library(testthat)
setwd('~/Dropbox/Duke/Spring 2015/PS 733/midterm')
source('midterm.R')
load('midTermData.rda')
library(Amelia)
set.seed(6886)
form = formula(gini_net_std ~ polity2 + ELF_ethnic)

olsSM(form, data)
ols(form, data)

expect_that(olsSM(form, data), equals(ols(form, data=data)))
expect_that(olsSM(form, dataMiss), equals(ols(form, data=dataMiss)))
set.seed(6886); expect_that(
  olsSM(form, dataMiss,impute=T), 
  equals(ols(form, data=dataMiss,impute=T)))

# ##################End of Function####################

# # Running the function with different models
# form = formula(gini_net_std ~ ELF_ethnic + polity2)

# model = ols(formula=form, data=data)

# modelListDel = ols(formula=form, data=dataMiss)

# modelAmelia = ols(formula = form, data=dataMiss, impute=TRUE)

# # Extracting coefficients, varcov, Rsq and fstat from the models

# ## Original model
# model$coefficients
# model$varcov
# model$Rsq
# model$fstat

# ## Listwise deleted model
# modelListDel$coefficients
# modelListDel$varcov
# modelListDel$Rsq
# modelListDel$fstat

# ## Imputed model
# modelAmelia$coefficients
# modelAmelia$varcov
# modelAmelia$Rsq
# modelAmelia$fstat

