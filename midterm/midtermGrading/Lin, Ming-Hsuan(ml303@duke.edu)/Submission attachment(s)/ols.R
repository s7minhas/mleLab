#Mid-term ols function
#Ming-Hsuan Lin

##################
ols=function(formula,data,impute=FALSE){
  #missing if miss=TURE
  if(impute==T){
    set.seed(6886)
    data=amelia(x=data, m=1)$imp$imp1}
  data=data[complete.cases(data),] 
  #(new)y,x,n,p,df
  dv = all.vars(formula)[1]
  ivs = all.vars(formula)[ 2:length(all.vars(formula)) ]
  y = data[,dv]
  x = data.matrix(cbind(1, data[,ivs]))
  n = nrow(x) # Number of observations
  p = length(ivs) # Number of parameters
  df = n-p-1 # degrees of freedom
  #coefficients
  xTx = t(x)%*%x
  xTy = t(x)%*%y
  beta = solve(xTx)%*%xTy
  #serrors
  yhat = x%*%beta
  e = y-yhat
  sigma2 = sum(e^2)/df
  varcov = sigma2*solve(xTx)
  serror = sqrt(diag(varcov)) 
  #t test, pvalue and CI
  tstate = beta/serror
  pvalue = round(2*pt(abs(tstate),df=n-p-1, lower.tail=F))
  up95=beta+qnorm(0.975)*serror
  lo95=beta-qnorm(0.975)*serror  
  #R2
  ssReg = sum((yhat-mean(y))^2)
  ssTot = sum((y-mean(y))^2)
  R2 = ssReg/ssTot
  #Fstat and pvalue
  msReg = sum((yhat-mean(y))^2)/p
  msRes = sum(e^2)/df
  Fstat = round(msReg/msRes,3) 
  Fstatp = round(pf(Fstat,p,df,lower.tail=F),3)
  xx=paste("F-statistic:",Fstat,"on",p,"and",df,"DF,","p-value:",Fstatp,sep=' ') 
  #creat matrix and names
  coese=cbind(beta, serror, tstate,pvalue,lo95,up95)
  colnames(coese)=c("Estimate","Std.Error","T-Statistic","P-Value",
                    "Lower 95% CI","Upper 95% CI")
  rownames(coese)=c("intercept", "Elf-ethnic","polity2")
  colnames(varcov)=c("intercept","Elf-ethnic","polity2")
  rownames(varcov)=c("intercept","Elf-ethnic","polity2")
  #list
  list1=list("coefficients"=coese, "varcov"=varcov, "Rsq"=R2, "Fstat"=xx)
  return(list1)   
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