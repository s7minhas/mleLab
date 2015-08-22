rm(list=ls())
# ---
#   title: "MLE : Mid-term"
# author: "Jinhee Kwak"
# date: "March 17, 2015"
# header-includes:
# - \usepackage{multirow}
# - \usepackage{dcolumn}
# output: 
#   pdf_document:
#   fig_caption: yes
# ---


ols=function(formula,data,impute=FALSE){
  
  if(impute==TRUE){
    set.seed(6886)
    data=amelia(x=data,m=1)
    data=data$imp$imp1}
  data=data[complete.cases(data),]
  
  # Retrieve vars from formula input
  dv = all.vars(form)[1]
  ivs = all.vars(form)[ 2:length(all.vars(form)) ]
  # Create matrix with column for intercept and
  ## data from independent variables
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  # General parameters
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df=n-p-1# degrees of freedom
  
  #coefficient
  
  # Beta = (X'X)^-1 %*% X'Y
  # calculating X'X
  xTx = t(x) %*% x
  
  # calculating X'y
  xTy = t(x) %*% y
  
  # calculating Beta
  beta = solve(xTx) %*% xTy
  
  #standard errors
  #se(Beta) = sqrt( sigma^2 * (X'X)^-1 )
  # First lets calculate our yhat
  yhat=x%*%beta
  # First get out residuals
  e=y-yhat
  # calculating e'e (sum of squared residuals): assuming homoskedasticity
  # Adjust by degrees of freedom
  sigma2 = sum(e^2)/df
  varcov = sigma2*solve(xTx)
  # Pull out the standard errors for the coefficient estimates
  se = sqrt(diag(varcov))
  # Calculate t values
  tval = beta/se
  # Calculate p-values
  pval = 2*pt(abs(tval),df, lower.tail=FALSE)
  up95=beta+qt(0.975,df)*se
  lo95=beta-qt(0.975,df)*se
  # R squared
  ssReg = sum((yhat-mean(y))^2)
  ssTot = sum((y-mean(y))^2)
  R2 = ssReg/ssTot
  # F statistic
  msReg = sum((yhat-mean(y))^2)/p
  msRes = sum(e^2)/df
  Fstat = round(msReg/msRes,3)
  Fpval = round(pf(Fstat,p,df,lower.tail=F),3)
  resul=paste("F-statistic:",Fstat,"on",p,"and",df,"DF,","p-value:",Fpval,sep=' ')
  # creating matrix
  coef=cbind(beta, se, tval, pval,lo95, up95)
  colnames(coef)=c("Estimate","Std. Error","T-Statistic","P-Value",
                   "Lower 95% CI","Upper 95% CI")
  rownames(coef)[1]=c("(Intercept)")
  colnames(varcov)=rownames(coef)
  rownames(varcov)=rownames(coef)
  a<-list("coefficients"=coef, "varcov"=varcov, "Rsq"=R2, "Fstat"=resul)
  return(a)
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