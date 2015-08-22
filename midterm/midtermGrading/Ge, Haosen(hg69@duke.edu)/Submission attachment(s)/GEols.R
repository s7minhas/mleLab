# loadPkg=function(toLoad){
#   for(lib in toLoad){
#     if(! lib %in% installed.packages()[,1])
#     { install.packages(lib, repos='http://cran.rstudio.com/') }
#     suppressMessages( library(lib, character.only=TRUE) ) }
# }
# setwd("C:/Users/Richard/Desktop/MLE/HW7")
# load("midTermData.rda")
# packs=c('foreign', 'lmtest', 'sandwich', 'Amelia', 'sbgcop','stargazer','ggplot2','gridExtra')
# loadPkg(packs)


# OLS function
  ols = function(formula, data, impute = FALSE){
  if(impute==TRUE){
  set.seed(6886)
  data=amelia(x=data,m=1)$imp$imp1}
  data=data[complete.cases(data),]
  dv = all.vars(formula)[1]
  ivs = all.vars(formula)[2:length(all.vars(formula))]
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df = n-p-1 # degrees of freedom
  ## Get the coefficients
  xTx = t(x) %*% x
  xTy = t(x) %*% y
  beta = solve(xTx) %*% xTy
  ## Get the se
  yhat=x %*% beta
  e= y - yhat
  sigma2=sum(e^2)/df
  vcov=sigma2 * solve(xTx)
  se=sqrt(diag(vcov))
  ## Get the tstat, p-val and CI
  tstat = beta/se
  pval1=2*pt(abs(tstat),df,lower.tail=FALSE)
  up95=beta+qt(0.975,df)*se
  lo95=beta-qt(0.975,df)*se
  ## Get R-sq
  SSreg=sum((yhat-mean(y))^2)
  SStot=sum((y-mean(y))^2)
  R2=SSreg/SStot
  ## Get F stat
  Fstat=round((sum((yhat-mean(y))^2)/p)/(sum(e^2)/df),3)
  pval2=round(pf(Fstat,p,df,lower.tail=FALSE),3)
  ## Prepare the result
  coefficients = cbind(beta,se,tstat,pval1,lo95,up95)
  colnames(coefficients)=c("Estimate", "Std. Error", "T-Statistic", "P-Value", "Lower 95% CI", "Upper 95% CI")
  rownames(coefficients)=c('(Intercept)',ivs)
  varcov=vcov
  rownames(varcov)=colnames(varcov)=rownames(coefficients)
  Rsq=R2
  Fstat=paste('F-statistic:',Fstat, 'on' ,p, 'and' ,df, 'DF,',  'p-value:', pval2,sep=' ')
  result = list(coefficients=coefficients,varcov=varcov,Rsq=Rsq,Fstat=Fstat)
  return(result)
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